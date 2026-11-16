> Historical v20 guide. Use [REPRODUCE_v20_1.md](REPRODUCE_v20_1.md) for the consolidated release and corrected read-only verifier.

# Reproduce the v20 downstream revision

Use a writable directory adjacent to this read-only candidate for outputs.

1. Check all payload hashes and selected downstream arithmetic:

   `python3 verify_bundle.py`

2. Refit the recorded Shannon, focused beta and dispersion models:

   `Rscript --vanilla analysis_scripts/refit_recorded_models.R .`

   This writes to `reproduced/fresh_model_check/`; retain prior runs separately.

3. Regenerate all 19 statistical PNG figures using the eight recorded plotting
   scripts and compare them with the selected current figures:

   `python3 analysis_scripts/replay_recorded_figures.py --candidate . --output ../v20_figure_replay --rscript Rscript`

   Python Pillow and the R packages listed in the recorded sessions are needed.
   Four nonstatistical original assets are preserved rather than regenerated.

4. Recompute the joint factorial model, all 9,999 permutations and 1,000 cluster
   bootstrap draws per arm, and run the independent NumPy verification:

   `python3 analysis_scripts/run_joint_factorial_q15.py --candidate . --output ../v20_joint_replay --rscript Rscript --verify`

5. Independently verify the bundled joint results without rerunning R:

   `python3 analysis_scripts/verify_Q15_joint_factorial_2026_11_05.py bundle/joint_factorial_q15`

R 4.4.1, vegan 2.7-2 and permute 0.9-8 were used for the new joint analysis.
The independent checker requires NumPy and pandas; it does not need SciPy.
The environment directory records other model/figure dependencies and versions.
These are tested local replay commands, not a claim of a fresh-machine package
installation test or of FASTQ-to-table reproduction. The source archive retains
historical path conventions; use the supplied runners for the supported replays.

The manuscript's direct quantitative DAPI and upstream sequencing limitations
remain. Successful statistical replay cannot supply missing experimental controls.
