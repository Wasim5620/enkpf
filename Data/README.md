# Data provenance

The observation workbooks contain public COVID-19 series compiled by
[Coronaboard](https://coronaboard.kr/en/) from its documented upstream providers.
Those include Korean public-health reporting, WHO, Worldometer, NHK and Japan's
Ministry of Health, Labour and Welfare. Follow the portal's reference links.
Source documentation was consulted on 16 October 2026; this is not the original
data-download date, which was not recorded. Observation dates describe the
included snapshots, not their retrieval date. No new third-party data license
is assigned by this package; applicable terms remain with the providers.

Source_Manifest.json and Source_Manifest.csv record project-relative original
names and hashes, packaged names and hashes, and numerical-array extraction.
The observation and fitted-result workbooks are copied unchanged. NativeForecast
contains all nine unchanged cache files required for the three original runs,
fifteen computational repetitions, and source-reference checks. All countries,
methods, outcomes and declared seeds remain included.

Reduced MAT files contain final EnKF endpoint ensembles, params1 and N_pop for
optional Gaussianity figure regeneration, with numerical values unchanged.
They contain data only. Initialization JSON records saved-run metadata; its
original script is identified by name/hash but no MATLAB script is included.
Validation contains saved synthetic recovery, local-information and policy
summaries. The original fitting, filtering, forecasting, sensitivity simulations
and synthetic estimation computations are not reproduced by this release.
