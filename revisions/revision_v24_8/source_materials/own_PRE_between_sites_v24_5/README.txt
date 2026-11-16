Direct between-site own-PRE follow-up for Reviewer 2 Major comment 2 (Q15).

Run: python own_PRE_between_sites_v24_5.py --inputs inputs --output replay
Dependencies: Python, NumPy, pandas. Original source counts/metadata are copied unchanged. Primary full-cohort 39-assignment nine-test family is displayed in Table S2R.11; related 41-assignment/proliferative sensitivities are in results. All pairs, exact sign assignments and 5000 bootstrap indices are retained. Matched means weight women equally; HUT099 is not counted twice. No original data or microscopy was edited. See Analysis_specification.json for all assumptions and limitations.

Independent audit: python verify_own_PRE_between_sites_v24_5.py --workdir .
