# Reproduction environment

The final Python rerun used Python 3.11.5 with the exact packages in
`requirements.txt`. The replicate reconstruction used R 4.4.0; the complete R
and Bioconductor package versions are frozen in
`results/reconstruction_session_info.txt`. The manuscript was compiled with
pdfTeX/MiKTeX 25.12.

Install Python dependencies in a clean environment with:

```powershell
python -m venv .venv
.\.venv\Scripts\python -m pip install -r requirements.txt
```

The R reconstruction additionally requires `GenomicRanges` and `data.table`.
The archived source objects listed in `input_provenance.csv` are deliberately
not bundled while redistribution permission remains unresolved.
