# Sample data

The recordings in this folder are third-party data. They are **not** covered by
the GPL-3.0 license of the code; each one keeps its original license, listed below.

| File | Source | Rate | Columns | License |
|---|---|---|---|---|
| `bidmc_01_Signals.csv` | BIDMC PPG and Respiration Dataset, subject 01 | 125 Hz | `Time [s]`, `RESP`, `PLETH`, `V`, `AVR`, `II` | ODC-By 1.0 |
| `Sample1.CSV`, `Sample2.CSV` | Pulse Transit Time PPG Dataset, 1.1.0 (excerpts; IR, distal phalanx, seated) | 500 Hz | `Time (s)`, `Pleth` (`;` separator, decimal comma) | ODbL 1.0 |

## Notes

- **BIDMC.** Pimentel, M. A. F. et al. (2017). *Toward a robust estimation of respiratory rate from pulse
  oximeters.* IEEE Trans. Biomed. Eng. 64(8), 1914–1923. PhysioNet: <https://physionet.org/content/bidmc/>.
  The dataset is based on the MIMIC-II matched waveform database.
- **Pulse Transit Time PPG Dataset.** Mehrgardt, P., Khushi, M., Poon, S., & Withana, A. (2022).
  *Pulse Transit Time PPG Dataset* (version 1.1.0). PhysioNet. RRID:SCR_007345.
  <https://doi.org/10.13026/jpan-6n92>. Released under the Open Database License v1.0
  (<https://opendatacommons.org/licenses/odbl/1-0/>). `Sample1.CSV` and `Sample2.CSV` are extracts of the
  dataset. Any adapted database made from them must be shared under the same license.
  - Acquisition: infrared (IR) light at the distal phalanx, subject seated at rest.
  - The `Pleth` column holds raw MAX30101 counts. Higher counts mean less absorbed light, so the pulse is
    **inverted**. The pipeline detects this and flips the signal (`preprocessing.polarity = "auto"`).
  - Around 23 % of consecutive samples repeat the previous value exactly. The sensor runs at
    1000 Hz and the stored data are synchronized at 500 Hz, so the repeats probably come from that
    resynchronization. They do not affect the analysis once the signal is resampled to the processing
    rate.
- **PhysioNet citation.** Goldberger, A. L. et al. (2000). *PhysioBank, PhysioToolkit, and PhysioNet.*
  Circulation 101(23), e215–e220.
- Only recordings whose license allows redistribution belong in this folder. Other local files can be
  kept here for your own analyses, but add them to `.gitignore` so that they are never committed.
