# Optional plotms patch

- `calibration.py` no longer imports or invokes `plotms` in the active path.
- `selfcal.py` no longer imports or invokes `plotms` in the active path.
- Native, headless QA plots from `qa_products.py` remain enabled.
- Missing or broken ilifu plotting components cannot abort calibration or self-calibration.
