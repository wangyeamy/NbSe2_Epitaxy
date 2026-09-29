# Interfacial forces and mechanical work in confined NbSe2

Analysis and figure-reproduction code accompanying **Resolving interfacial mechanics of monolayer NbSe2 grown by confined epitaxy** (Ye Wang).

Version: 1.0.0-review, 30 September 2026. This is a local review package; no repository URL or DOI has yet been assigned.

## Run

Python 3.12 was used for verification. From this directory:

```sh
python -m venv .venv
# Activate the virtual environment using the command appropriate to your OS.
python -m pip install -r requirements.txt
python reproduce.py
```

The last command needs neither network access nor Quantum ESPRESSO. It reads the supplied files and writes plots (PNG, SVG, PDF), numerical tables (CSV) and audit results (JSON) to `results/`. Installing the dependencies normally requires internet access. No API key, institutional account, server connection or AI service is required. Typical local runtime is under one minute; hardware will affect this.

Plots use DejaVu Sans for portability. This may change typography relative to the manuscript, whose figures use Arial; it does not change data. Generated reference outputs are included for comparison.

## Manuscript mapping

| Item | Script | Input | Output |
| --- | --- | --- | --- |
| Figure 1 | Conceptual workflow/structure artwork; no numerical analysis | Not part of this numerical code package | — |
| Figure 2 | `scripts/plot_figure2.py` | Archived low-moment atom-resolved force CSV | `results/Figure2.*` |
| Figure 3, latest three-panel revision | `scripts/plot_figure3.py` | Accepted high-moment endpoints and actual reference atomic coordinates | `results/Figure3.*` |
| Figure 4, DOS | `scripts/plot_figure4.py` | Three original DOS output tables | `results/Figure4.*`, `dos_audit.json` |
| Figure S1 | `scripts/plot_figureS1.py` | Plane means recovered by `analyze.py` from archived mode projections | `results/FigureS1.*` |
| Figure S2 | `scripts/plot_figureS2.py` | Archived high-moment four-state summary and recomputed five-mode work | `results/FigureS2.*` |
| Tables S1–S3 and force/work checks | `scripts/analyze.py` | QE output files, benchmark CSV, endpoint CSV, mode summary | `results/TableS1.csv`–`TableS3.csv`, `analysis.json` |

An earlier four-panel Figure 3 can be generated separately with `python scripts/plot_figure3_legacy.py`. That version contains the work-comparison panel and differs from the author's latest three-panel figure. Its intermediate curve is interpolation, not additional calculations. Both versions use the same accepted endpoint data. The default reproduction uses the latest three-panel revision. These scripts do not modify any manuscript document.

## What is reproduced

1. The low-moment matched confined-minus-unconfined force differences are independently parsed from the final complete total-force blocks of six supplied QE outputs (two meshes, each with confined/unconfined/isolated-hBN cases). The results are checked against the atom-resolved CSV underlying Figure 2. Verbose QE output also prints component-force blocks; these are not substituted for total forces.
2. The high-moment plane means are recovered by solving the three linear relations between the normalized Nb–Se, rigid-z and Se-breathing projections and the three plane sums. This recovers averages, not individual atomic forces.
3. The work is computed by the endpoint trapezoidal rule. With `f_int = F_confined - F_unconfined`, `W = integral(f_int dq)` along q = 0 to -0.02 angstrom. Thus the comparison is `deltaI` against `-W`. The results are +6.712504750751146 and +6.709269341832361 meV/cell; their difference is 0.0032354089187851542 meV/cell. It is not an uncertainty bound.
4. Figure 3's coordinate inputs are checked: the common atoms have identical reference positions, and the supplied negative displacements match the normalized mode. Nb moves by -0.02/sqrt(24) angstrom, each Se by +0.01/sqrt(24) angstrom; all other atoms remain fixed. The plotted arrows show collective displacement direction, not force. Six nearest Se neighbors per Nb are verified with periodic translations and a 2.95 angstrom cutoff. The displayed cluster is a crop of a periodic crystal, not an isolated molecule or its stoichiometric formula.
5. The DOS is summed over spin channels and referenced to each calculation's own Fermi energy. The value at zero energy is linearly interpolated from actual samples; spectra are not smoothed or put on a common absolute energy reference.
6. The four-state five-mode work diagnostic is recomputed from mode projections and geometric displacements. It is separate from the single-coordinate work comparison.

## Data provenance and coverage

`data/forces/atom_forces.csv` is the archived matched-coordinate low-moment export. It retains the original `capped`/`uncapped` terminology in source-field values where present; throughout the paper these mean confined/unconfined. Source JSON keys are retained for traceability.

`data/qe_outputs/` contains six archived low-moment outputs with numerical content preserved and local filesystem paths redacted. `data/dos/` contains original numerical DOS outputs. `data/high_moment/four_state.json` is the archived high-moment summary, not a newly recomputed raw-output analysis. `data/high_moment/displacement.csv` records only the two accepted endpoints. The high-moment per-atom output files and saved electronic checkpoints are not included in this local package. Consequently, the package reproduces the high-moment analyses from their recorded inputs/summaries; it does not independently regenerate those electronic solutions.

`data/benchmark/` holds the archived primitive-cell CSV and its inputs; the individual benchmark output blocks were not recovered. `historical/high_spin_four_state.py` preserves the original four-state analysis program with site paths redacted. It is included for method inspection, is not part of the default runner, and needs the original full output/checkpoint-marker collection to run.

QE inputs in `data/inputs/` are reference records. Site-specific paths are redacted, pseudopotential binaries and restart checkpoints are not bundled, and these are not advertised as ready-to-submit cluster jobs. No unconverged point has been relabeled as converged. No incorporation path or reaction barrier is generated by this code.

`data_manifest.json` records original and distributed hashes for copied inputs. `SHA256SUMS.txt` covers all package files other than itself. Source fingerprints refer to the archived originals; redacted copies have their own fingerprints.

## Verification

The default runner was executed successfully locally with the versions in `results/environment.json`. Numerical checks include the CSV-versus-raw-force comparison, matched reference geometry, prescribed endpoint displacement, reconstruction of force projections, and the endpoint force/energy relationship. These checks validate this analysis package; they are not a DFT convergence study.

## Authorship and availability

Author/contact: Ye Wang, Eindhoven University of Technology, y.wang19@tue.nl.

AI tools assisted analysis and code preparation under human scientific direction, as disclosed in the manuscript. The review scripts are deterministic and do not call an AI model.

This package is prepared for inspection and execution during peer review. A general software reuse license has not yet been selected by the author. GitHub/Zenodo access links and a DOI should be added only after deposit. No link in this package implies that publication or deposit has already occurred.
