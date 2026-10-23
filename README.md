# Residual-weighted ensemble Kalman–particle filtering for behavioral COVID-19 modeling and retrospective forecasting

Publication figures and tables from saved numerical exports for Korea, Japan,
and Pakistan. Repository: https://github.com/Wasim5620/enkpf

## Contents and use

- `EnKPF_Manuscript_Supplementary.ipynb` is self-contained: eight Python helpers
  are embedded in isolated namespaces. No separate project Python files are needed.
- `Data/` contains observation workbooks, fitted-result workbooks, reduced
  endpoint ensemble arrays, native forecast CSVs, saved recovery/policy summaries,
  and caption/provenance metadata. MAT files contain numerical data only.
- `Results/Tables/` contains four main and twenty supplementary numbered CSVs,
  the 25-sheet workbook (Index plus 24 tables), and output/source records.
- `Results/manuscript_tables/` contains the corresponding 24 TeX fragments.
- `Figures/PNG/` contains seven main and seven supplementary PNG/PDF figure sets;
  `Figures/TIFF/` contains the same fourteen RGB LZW-compressed TIFFs.
- `Results/Forecast/` records the figure clipping audit, captions, provenance,
  and both completed-origin and convergence-subset forecast summaries.

Use Python 3.11 or newer. In this folder:

```sh
python -m venv .venv
# Activate the environment using the command for your operating system.
python -m pip install -r requirements.txt
python -m jupyter lab EnKPF_Manuscript_Supplementary.ipynb
```

Select that environment's kernel, then Restart Kernel and Run All. The notebook
finds this folder from the working directory or its parents; set
`PACKAGE_ROOT_OVERRIDE` if needed. All input/output paths are package-relative.
Default cells read saved data, export current tables, display reviewed figures,
and redraw the two native forecast views from CSVs. `REBUILD_EXISTING_FIGURES`
also redraws the other fitting/parameter/diagnostic plots from saved arrays.
Sensitivity and synthetic recovery results are read from their saved summaries;
no ODE sensitivity, fitting, filtering, resampling, or forecast simulation runs.
No MATLAB code or installation is included or required. This package reproduces
publication outputs from exports, not the original MATLAB estimation experiments.

## Scope and limitations

The primary retrospective comparison includes all three completed peak origins,
six original Q/R outcomes, and all five computational seeds per country.
EnKPF has lower original RMSE in six of six comparisons and lower 14-day RMSE in
21 of 30 matched comparisons. The separately reported Korea/Japan converged-
initializer robustness subset is 16 of 20. Pakistan's initializer reaches its
evaluation budget (flag 0; 1,501 evaluations); its completed, finite forecasts
are retained with that limitation disclosed. These scopes were defined during
revision and are not independent prospective validation. Repeated seeds describe
computational variability, not independent epidemiological replications.

Korean original active-case R-squared remains negative. Latent interval coverage
is imperfect, longer-horizon aggregate errors exceed deterministic references,
and a completed alternative benchmark did not establish a general hybrid
advantage. Incomplete extensions and unresolved profiles are not publication
results in this release.

Training means in Figures 2/S7 are solid lines without uncertainty shading.
Central empirical 95% latent bands begin at the first held-out observation;
S7 member extrema are held-out only. Figure 2 uses disclosed observed-scale
clipping, with all full ranges available in S7. Future measurement noise is
excluded. Figure 1 separately uses full-wave mean plus/minus 1.96 ensemble SD.
Behavioral scenarios compare equal +5% changes at 10% versus 30% elapsed duration;
they are conditional simulations, not causal policy-effect estimates.

The notebook and cached figure assets follow the approved manuscript title above.
See `Data/README.md`, `Data/Source_Manifest.json`, and `Package_Contents.json`
for data attribution, source hashes, extraction details and release checks.
No files have been uploaded or submitted by this build.
