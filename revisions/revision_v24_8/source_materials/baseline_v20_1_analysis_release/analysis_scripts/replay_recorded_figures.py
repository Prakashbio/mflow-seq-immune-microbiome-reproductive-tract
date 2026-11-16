"""Replay recorded plotting scripts in an isolated layout; never modify inputs.

Requires Python Pillow and R packages recorded in the source environment.
No experimental validation or upstream FASTQ reprocessing is performed.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import sys
from PIL import Image, ImageChops, ImageStat


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def inventory(root):
    return {str(p.relative_to(root)): digest(p) for p in sorted(root.rglob('*')) if p.is_file()}


def write_csv(path, rows):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--candidate', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--rscript', default='Rscript')
    args = ap.parse_args()
    source = args.candidate.resolve()
    out = args.output.resolve()
    if source == out or source in out.parents:
        ap.error('Output must be outside the read-only candidate.')
    if out.exists():
        ap.error('Output must be a new directory; previous runs are preserved.')
    assert (source / 'documentation/Current_figure_exports.csv').is_file()
    before = inventory(source)
    out.mkdir(parents=True)
    stage = out / 'stage'
    # All numeric inputs, saved values and recorded source scripts are copied
    # exactly. The historical arguments then point into this isolated layout.
    shutil.copytree(source / 'bundle', stage)
    shutil.copytree(source / 'analysis_scripts/source_archive', out / 'recorded_scripts')
    (out / 'logs').mkdir()
    (out / 'selected_replayed_figures').mkdir()
    scripts = out / 'recorded_scripts'
    jobs = [
        ('PRE context', 'verify_and_plot_PRE_v7_2026_10_29.R', [stage / 'pre_context']),
        ('FACS', 'plot_facs_publication_v12_2026_10_30.R', [stage]),
        ('Diversity', 'plot_diversity_integration_2026_10_28.R', [stage, stage / 'statistics']),
        ('Taxon representation', 'plot_taxon_integration_v3_2026_10_28.R', [stage, stage / 'taxon_analysis', stage / 'statistics']),
        ('Legacy supplementary descriptions', 'verify_plot_legacy_supplement_v6_2026_10_29.R', [stage]),
        ('Analytical feature selection', 'plot_verify_source_filter_v5_2026_10_28.R', [stage]),
        ('Storage descriptions', 'plot_verify_storage_review_v4_2026_10_28.R', [stage]),
        ('Confirmed storage-pair notes', 'update_storage_figure_notes_v14_2026_11_03.R', [stage]),
    ]
    log_rows = []
    for n, (component, name, argv) in enumerate(jobs, 1):
        print(f'{n}/{len(jobs)} {component}', flush=True)
        command = [args.rscript, '--vanilla', str(scripts / name), *map(str, argv)]
        run = subprocess.run(command, cwd=stage, capture_output=True, text=True)
        logfile = out / 'logs' / f'{n:02d}_{name}.txt'
        logfile.write_text(run.stdout + '\nSTDERR:\n' + run.stderr)
        log_rows.append(dict(component=component, script=name, source_sha256=digest(source / 'analysis_scripts/source_archive' / name), returncode=run.returncode, log=str(logfile.relative_to(out))))
        write_csv(out / 'Script_execution.csv', log_rows)
        if run.returncode:
            raise RuntimeError(f'{component} failed; see {logfile}')

    # Current figure numbers differ from historical plot-output filenames.
    # S9 is selected by exact raster agreement among the two recorded storage
    # displays, avoiding an assumption from an old crosswalk description.
    mapping = {
        'Figure 2': ['pre_context/results/PRE_baseline_context_REVIEW.png'],
        'Figure 3': ['figure/Main_Figure_3_v12.png'],
        'Figure 4': ['figures/Figure2_Shannon_preparation_pairs.png'],
        'Figure 5': ['figures/Figure3_PCoA_participant_pairs.png'],
        'Figure 6': ['figures/Figure5_matched_site_scores.png'],
        'Figure S4': ['figures/FigureS14_filtering_reconstruction.png'],
        'Figure S5': ['figures/FigureS6_PRE_descriptive_rank_correlations.png'],
        'Figure S6': ['figures/FigureS8_PRE_participant_accumulation.png'],
        'Figure S7': ['figures/FigureS11_PRE_and_retained_count_descriptions.png'],
        'Figure S8': ['figures/FigureS2_storage_alpha_descriptive.png'],
        'Figure S9': ['figures/Figure7_descriptive_storage.png', 'figures/FigureS7_storage_group_means.png'],
        'Figure S10': ['figure_note_update_v14/Supplementary_Figure_S10.png'],
        'Figure S11': ['figure_note_update_v14/Supplementary_Figure_S11.png'],
        'Figure S12': ['figures/FigureS3_Shannon_woman_site_pairs.png'],
        'Figure S13': ['figures/FigureS4_matched_own_PRE.png'],
        'Figure S14': ['figures/FigureS9_matched_PRE_contrast_intervals.png'],
        'Figure S15': ['figures/FigureS12_Lactobacillus_paired_values.png'],
        'Figure S16': ['figures/FigureS10_genus_representation.png'],
        'Figure S17': ['figures/Figure6_PRE_group_descriptives.png'],
    }
    registry = list(csv.DictReader((source / 'documentation/Current_figure_exports.csv').open()))
    results = []
    for row in registry:
        label = row['current_label']
        if label not in mapping:
            results.append(dict(figure=label, status='PRESERVED_NONSTATISTICAL_ASSET', recorded_output='', published_sha256=row['sha256'], replayed_sha256='', size_equal='', pixels_equal='', mean_absolute_channel_difference='', note='Workflow diagram or original native experimental panel; not regenerated.'))
            continue
        original = source / 'bundle' / row['asset']
        assert digest(original) == row['sha256']
        best = None
        for name in mapping[label]:
            replay = stage / name
            with Image.open(original) as oa, Image.open(replay) as ra:
                a, b = oa.convert('RGBA'), ra.convert('RGBA')
                equal_size = a.size == b.size
                equal = equal_size and a.tobytes() == b.tobytes()
                mad = sum(ImageStat.Stat(ImageChops.difference(a, b)).mean) / 4 if equal_size else None
            candidate = dict(figure=label, status='PASS_IDENTICAL_PIXELS' if equal else 'REVIEW_RASTER_DIFFERENCE', recorded_output=name, published_sha256=digest(original), replayed_sha256=digest(replay), size_equal=equal_size, pixels_equal=equal, mean_absolute_channel_difference=mad, note='Exact RGBA comparison; PDF metadata/bytes are not treated as plot identity.')
            if best is None or equal:
                best = candidate
            if equal:
                break
        results.append(best)
        if best['pixels_equal']:
            filename = Path(row['asset']).name
            shutil.copyfile(stage / best['recorded_output'], out / 'selected_replayed_figures' / filename)
            for ext in ('.pdf', '.svg'):
                vector = (stage / best['recorded_output']).with_suffix(ext)
                if vector.exists():
                    shutil.copyfile(vector, (out / 'selected_replayed_figures' / filename).with_suffix(ext))
    write_csv(out / 'Figure_replay_comparison.csv', results)
    after = inventory(source)
    source_unchanged = before == after
    bad = [r['figure'] for r in results if r['status'] == 'REVIEW_RASTER_DIFFERENCE']
    summary = dict(status='PASS' if source_unchanged and not bad else 'REVIEW_REQUIRED', recorded_plot_scripts_run=len(jobs), figures_replayed=len(mapping), identical_raster_figures=sum(r['status'] == 'PASS_IDENTICAL_PIXELS' for r in results), preserved_nonstatistical_assets=4, source_candidate_files=len(before), source_candidate_unchanged=source_unchanged, figures_requiring_review=bad, raw_FASTQ_pipeline_reproduced=False, experimental_validation_added=False, manuscript_or_sealed_package_changed=False, remote_actions=False, submission_ready=False)
    (out / 'Verification.json').write_text(json.dumps(summary, indent=2) + '\n')
    hashes = {str(p.relative_to(out)): digest(p) for p in sorted(out.rglob('*')) if p.is_file() and p.name != 'Replay_SHA256SUMS.json'}
    (out / 'Replay_SHA256SUMS.json').write_text(json.dumps(hashes, indent=2) + '\n')
    print(json.dumps(summary, indent=2))
    return 0 if summary['status'] == 'PASS' else 1


if __name__ == '__main__':
    sys.exit(main())
