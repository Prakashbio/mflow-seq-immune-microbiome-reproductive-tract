#!/usr/bin/env python3
"""
run_all.py — Execute the complete mFLOW-Seq bioinformatics pipeline
====================================================================
Runs all analysis scripts in order, generating all publication figures
in the new_figures/ directory.

Usage:
    python run_all.py              # Run everything
    python run_all.py --fig 1 4 8  # Run only specified figures (1, 4, 8)
    python run_all.py --list       # List all available scripts
    python run_all.py --help       # Show this help message

Author: mFLOW-Seq Analysis Pipeline
"""

import os
import sys
import subprocess
import argparse
import time
from pathlib import Path

# Get the scripts directory
SCRIPT_DIR = Path(__file__).parent / 'scripts'
FIG_DIR = Path(__file__).parent / 'new_figures'

# Define all figure scripts in logical execution order
FIGURE_SCRIPTS = [
    # Main Figures
    ('2', 'fig2_alpha_diversity.py', 'Figure2.png'),
    ('3', 'fig3_beta_diversity.py', 'Figure3.png'),
    ('4', 'fig4_facs_binding.py', 'Figure4.png'),
    ('5', 'fig5_differential_abundance.py', 'Figure5.png'),
    ('6', 'fig6_healthy_vs_dysbiotic.py', 'Figure6.png'),
    ('7', 'Fig7_fresh_vs_frozen.py', 'Figure7.png'),
    ('8', 'fig8_anatomical_site.py', 'Figure8.png'),
    ('9', 'fig9_multivariate_analysis.py', 'Figure9.png'),

    # Supplementary Figures (ordered by number)
    ('S1', 'sfig1_phylum_composition.py', 'SupplementaryFigure1.png'),
    ('S2', 'sfig2_alpha_diversity_storage.py', 'SupplementaryFigure2.png'),
    ('S3', 'sfig3_per_location_pcoa.py', 'SupplementaryFigure3.png'),
    ('S4', 'sfig4_vs_presorted.py', 'SupplementaryFigure4.png'),
    ('S5', 'sfig5_condition_volcano.py', 'SupplementaryFigure5.png'),
    ('S6', 'sfig6_cooccurrence_matrix.py', 'SupplementaryFigure6.png'),
    ('S7', 'sfig7_fresh_frozen_correlation.py', 'SupplementaryFigure7.png'),
    ('S8', 'sfig8_taxon_accumulation.py', 'SupplementaryFigure8.png'),
    ('S9', 'sfig9_supplementary_igg.py', 'SupplementaryFigure9.png'),
    ('S10', 'sfig10_iga_coating_score.py', 'SupplementaryFigure10.png'),
    ('S11', 'sfig11_patient_paired_readdepth.py', 'SupplementaryFigure11.png'),
    ('S12', 'sfig12_lactobacillus_gradient.py', 'SupplementaryFigure12.png'),
    ('S13', 'sfig13_composition.py', 'SupplementaryFigure13.png'),
    ('S14', 'sfig14_prevalence_network.py', 'SupplementaryFigure14.png'),
    ('S15', 'sfig15_fresh_frozen.py', 'SupplementaryFigure15.png'),
    ('S16', 'sfig16_contamination_analysis.py', 'SupplementaryFigure16.png'),
]


def list_scripts():
    """Print a list of all available scripts."""
    print("\n" + "="*70)
    print("Available Figure Scripts")
    print("="*70)
    print(f"{'Fig #':<8} {'Script Name':<40} {'Output'}")
    print("-"*70)

    for fig_num, script, output in FIGURE_SCRIPTS:
        print(f"{fig_num:<8} {script:<40} {output}")

    print("="*70)
    print(f"\nTotal: {len(FIGURE_SCRIPTS)} scripts")
    print(f"Output directory: {FIG_DIR.absolute()}\n")


def run_script(fig_num, script_name, output_name, verbose=True):
    """
    Run a single figure script.

    Parameters
    ----------
    fig_num : str
        Figure number (e.g., '1', 'S1')
    script_name : str
        Name of the Python script
    output_name : str
        Expected output file name
    verbose : bool
        Whether to print verbose output

    Returns
    -------
    bool
        True if successful, False otherwise
    """
    script_path = SCRIPT_DIR / script_name

    if not script_path.exists():
        print(f"  ✗ ERROR: Script not found: {script_path}")
        return False

    if verbose:
        print(f"\n{'─'*70}")
        print(f"Running Figure {fig_num}: {script_name}")
        print(f"{'─'*70}")

    start_time = time.time()

    try:
        # Run the script
        result = subprocess.run(
            [sys.executable, str(script_path)],
            cwd=SCRIPT_DIR,
            capture_output=True,
            text=True,
            timeout=300  # 5 minute timeout
        )

        elapsed = time.time() - start_time

        if result.returncode == 0:
            # Check if output file was created
            output_path = FIG_DIR / output_name
            if output_path.exists():
                if verbose:
                    print(f"  ✓ SUCCESS: {output_name} created in {elapsed:.1f}s")
                    if result.stdout.strip():
                        print(f"  Output: {result.stdout.strip()}")
                return True
            else:
                print(f"  ✗ WARNING: Script completed but {output_name} not found")
                if verbose and result.stdout.strip():
                    print(f"  Output: {result.stdout.strip()}")
                if verbose and result.stderr.strip():
                    print(f"  Errors: {result.stderr.strip()}")
                return False
        else:
            print(f"  ✗ FAILED: Script exited with code {result.returncode}")
            if result.stderr.strip():
                print(f"  Error output:")
                for line in result.stderr.strip().split('\n')[:20]:  # Show first 20 lines
                    print(f"    {line}")
            return False

    except subprocess.TimeoutExpired:
        print(f"  ✗ TIMEOUT: Script took longer than 5 minutes")
        return False
    except Exception as e:
        print(f"  ✗ ERROR: {str(e)}")
        return False


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description='Run mFLOW-Seq figure generation scripts',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_all.py              # Run all scripts
  python run_all.py --fig 1 4 8  # Run figures 1, 4, and 8
  python run_all.py --fig S1 S2  # Run supplementary figures 1 and 2
  python run_all.py --list       # List all available scripts
        """
    )

    parser.add_argument(
        '--fig', '--figure',
        nargs='+',
        metavar='NUM',
        help='Run only specified figure(s). Use numbers (1-8) or S# for supplementary'
    )

    parser.add_argument(
        '--list',
        action='store_true',
        help='List all available scripts and exit'
    )

    parser.add_argument(
        '-q', '--quiet',
        action='store_true',
        help='Suppress detailed output'
    )

    args = parser.parse_args()

    # Handle --list option
    if args.list:
        list_scripts()
        return 0

    # Ensure output directory exists
    FIG_DIR.mkdir(exist_ok=True)

    # Filter scripts if specific figures requested
    if args.fig:
        requested = set(args.fig)
        scripts_to_run = [
            (num, script, output)
            for num, script, output in FIGURE_SCRIPTS
            if num in requested
        ]

        # Check if any requested figures were not found
        found_nums = set(num for num, _, _ in scripts_to_run)
        not_found = requested - found_nums
        if not_found:
            print(f"\nWarning: The following figures were not found: {', '.join(sorted(not_found))}")
            print("Use --list to see all available figures.\n")

        if not scripts_to_run:
            print("Error: No valid figures specified.")
            return 1
    else:
        scripts_to_run = FIGURE_SCRIPTS

    # Print header
    print("\n" + "="*70)
    print("mFLOW-Seq Figure Generation Pipeline")
    print("="*70)
    print(f"Running {len(scripts_to_run)} script(s)")
    print(f"Output directory: {FIG_DIR.absolute()}")
    print("="*70)

    # Run scripts
    start_time = time.time()
    results = []

    for fig_num, script_name, output_name in scripts_to_run:
        success = run_script(fig_num, script_name, output_name, verbose=not args.quiet)
        results.append((fig_num, script_name, output_name, success))

    # Print summary
    elapsed = time.time() - start_time
    successful = sum(1 for _, _, _, success in results if success)
    failed = len(results) - successful

    print("\n" + "="*70)
    print("Execution Summary")
    print("="*70)
    print(f"Total scripts:  {len(results)}")
    print(f"Successful:     {successful} ✓")
    print(f"Failed:         {failed} ✗")
    print(f"Total time:     {elapsed:.1f}s ({elapsed/60:.1f} min)")
    print("="*70)

    # Print detailed results
    if failed > 0:
        print("\nFailed scripts:")
        for fig_num, script_name, _, success in results:
            if not success:
                print(f"  ✗ Figure {fig_num}: {script_name}")

    print(f"\nOutput directory: {FIG_DIR.absolute()}")
    print(f"Generated figures: {successful}\n")

    return 0 if failed == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
