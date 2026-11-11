# Analyses for "Scaled Polyspherical PCA"

Run from the repository root (requires the `sppca` package sources in
`sppca/`, plus `devtools`, `polykde`, `shapes`, `rotasym`):

| Script | Output |
|---|---|
| `simulation_study.R` | Simulation study (`SPPCA_CORES` workers, default 4) |
| `sim_figures.R` | Simulation figures and summary table |
| `hippocampus.R`, `hippo_figs.R` | Hippocampus s-reps analysis and figures |
| `rna.R` | RNA torsion angles on T^7 (data in `st-pca/`) |
| `rho_profile.R` | Stress profile in the radius |
| `landmark_study.R` | Stress and time vs number of landmarks |
| `pipeline_figure.R` | Pipeline and OODA mode figures |
| `benchmark_old.R` | Original implementation (`R/`) vs `sppca` (slow) |
| `make_tables.R` | LaTeX tables |

Figures go to `analysis/figures/`, cached results to `analysis/results/`.
