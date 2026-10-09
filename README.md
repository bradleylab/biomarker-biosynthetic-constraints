# Biosynthetic constraints on molecular biomarker interpretation

Analysis code for Hinkston, M. and Bradley, A. S., *Quantifying biosynthetic constraints on molecular biomarker interpretation*. The analyses examine the documented biochemistry of 86 lipids and biosynthetic intermediates related to molecular fossils, using reaction networks built from MetaCyc 29.5 and curated literature records. Projected paths in a simplified compound graph are measured (path multiplicity ψ, normalized branch depth λ and fraction shared σ) and compared with molecular descriptors. The reactants of each producing reaction are traced upstream, and alternative reaction sets are enumerated within a search limit. Pairs of biomarkers are compared by the enzymes and pathways recorded for their producing reactions, and the two pairs shown in Figure 3B also by their projected paths. The code also finds which combinations of source categories can account for the sterane and hopanoid assemblages of 156 geological samples reported by Brocks et al. (2023).

## Requirements

Python 3.11 and the packages pinned in `requirements.lock`:

```bash
python3.11 -m venv .venv
source .venv/bin/activate
pip install -r requirements.lock
```

The analyses read molecular descriptors from `data/descriptors.csv`. Recalculating them is optional and needs further software. `analyses/02_projected_paths/extract_structures.py` reads compound structures through Pathway Tools, with the MetaCyc 29.5 database and the pythoncyc package. `analyses/02_projected_paths/compute_descriptors.py` uses the packages in `requirements-descriptors.txt` and the Böttcher scorer `bottchscore3.py` (https://github.com/forlilab/bottchscore, revision 6f1368cd), and checks structural eligibility with RDKit 2025.09.4 installed in a separate environment. Commands are in each script's docstring.

## Input data

Two inputs are not included:

- **MetaCyc 29.5** flat files `classes.dat`, `compounds.dat`, `enzrxns.dat`, `genes.dat`, `pathways.dat`, `proteins.dat` and `reactions.dat`, in one folder. MetaCyc is available from SRI International under license (https://metacyc.org).
- The two Supplementary Information workbooks of **Brocks et al. (2023)**, saved in one folder as `brocks2023_SI1_MOESM4.xlsx` and `brocks2023_SI2_MOESM5.xlsx`. Brocks, J. J. et al. Lost world of complex life and the late rise of the eukaryotic crown. *Nature* 618, 767–773 (2023). https://doi.org/10.1038/s41586-023-06170-w

Everything else is in `data/`:

| File | Contents |
|---|---|
| `biomarker_identities.jsonl` | The 86 biomarkers and the compound identities considered for each |
| `literature_reactions.json` | The 27 reaction records compiled from the literature (Table S8) and the literature-defined compounds they connect |
| `base_literature_reactions.json` | Literature reactions added to the base network; the final networks replace them with the records of `literature_reactions.json` and keep their compounds |
| `metacyc_corrections.json` | Decisions on 14 MetaCyc reactions, with reasons and sources |
| `compound_aliases.csv` | Literature compound identifiers and the MetaCyc compound each one denotes |
| `starting_compounds.csv` | Compounds at which tracing stops in the matched controls |
| `representative_targets.csv` | The reaction target that represents each precursor class in the geological analysis |
| `comparison_set.csv` | The 10,473 compounds of the projected-path comparison |
| `descriptors.csv` | Molecular descriptors, Böttcher scores and structural eligibility of the 7,380 compounds with descriptor records, with the path counts and cumulative branching scores of the reference graph that select the descriptor cohort |
| `literature_structures.csv` | Structures of nine literature-defined compounds from PubChem, used only by `compute_descriptors.py` |
| `geology/` | Workbook layout, sample matching between the two workbooks, fossil-to-precursor and precursor-to-source assignments, source categories and assignment criteria, source profiles, citations and analysis settings |

## Running the analyses

```bash
make all METACYC=/path/to/metacyc-29.5 WORKBOOKS=/path/to/brocks2023
```

`METACYC` is the folder of MetaCyc flat files and `WORKBOOKS` the folder of the two workbooks. With the default of four worker processes (`WORKERS`), a full run takes about 20 minutes and writes about 1.3 GB to the output folder (`RESULTS`, default `results`).

Steps can also be run one at a time (for example, `make tracing`) once the steps before them have finished; `make sources` needs only `WORKBOOKS`.

| Step | Command | Output folder |
|---|---|---|
| Reaction networks | `make networks` | `results/networks/` |
| Enzyme and pathway comparisons | `make pairs` | `results/pairs/` |
| Projected paths and descriptors | `make projected_paths` | `results/projected_paths/`, `results/descriptors/` |
| Upstream tracing | `make tracing` | `results/tracing/` |
| Source compatibility | `make sources` | `results/sources/` |
| Literature records | `make literature` | `results/literature/` |
| Matched controls | `make controls` | `results/controls/` |
| Figures | `make figures` | `results/figures/` |
| Supplementary Data | `make supplementary_data` | `results/supplementary_data/` |
| Numbers in the text | `make numbers` | printed |

The `Makefile` lists the scripts behind each target, and a script can be run on its own from the repository root with the environment variables set at the top of the `Makefile`. Inputs and outputs are described in each script's docstring and by `--help`.

## Outputs

The Supplementary Data tables are collected in `results/supplementary_data/`:

| Table | Files |
|---|---|
| Table 1 | `Table1_projection_class_means.csv`, `Table1.md` |
| Table S1 | `TableS1_*.csv` (five files) |
| Table S2 | `TableS2_*.csv` (four files) |
| Table S3 | `TableS3_pair_comparisons.csv` |
| Table S4 | `TableS4_*.csv` (four files) |
| Table S5 | `TableS5_source_observations.csv` |
| Table S6 | `TableS6_*.csv` (six files) |
| Table S7 | `TableS7_*.csv` (nine files) |
| Table S8 | `TableS8_*.csv` (three files) |
| Table S9 | `TableS9_*.csv` (three files) |

The figures are written to `results/figures/`, each with a JSON file of the plotted values:

| Figure | Files |
|---|---|
| Figure 2B | `Figure2_worked_example.pdf` |
| Figure 2C | `Figure2_diploptene.pdf`, `Figure2_cholesterol.pdf` |
| Figure 3 | `Figure3.pdf`, `Figure3.png` |
| Figure 4 | `Figure4.pdf`, `Figure4.png` |
| Supplementary Figures 1–6 | `SupplementaryFigure1` to `SupplementaryFigure6` (PDF and PNG), and the detail pages of Supplementary Figure 2 (`SupplementaryFigure2_DETAIL.pdf`) |

Supplementary Figures 1, 5 and 6 are also drawn for the graph with conditional connections (`*_with_conditional`).

The following data records are archived with the code on Zenodo:

| File | Contents |
|---|---|
| `results/sources/detected_class_sets.jsonl` | All compatible combinations of source categories and all sensitivity tests for each set of detected fossil classes |
| `results/sources/sample_scenarios.jsonl` | Every sample scenario (detection threshold, low-signal treatment and assignment criteria) and its set of detected classes |
| `results/sources/samples.jsonl` | Workbook rows, labels and reported ages of each sample |
| `results/sources/delegation_summary.csv` | For each set of detected classes, the smallest number of classes attributed to sources outside the five categories that changes each conclusion |
| `results/sources/delegated_class_sets.csv` | Every smallest such attribution behind `TableS7_delegation_by_sample.csv` |
| `results/descriptors/TableS4_diagnostic_sensitivity.csv` | Descriptor statistics under each treatment of the Böttcher score, summarized in `TableS4_score_treatment_summary.csv` |
| `results/projected_paths/Projection_graph_and_comparison_scope.csv` | Graph and comparison-set counts given in the Supplementary Methods |

## Tests

```bash
make test
```

## License

The code is released under the MIT license (`LICENSE`). MetaCyc and the Brocks et al. (2023) workbooks are not included and remain under their own terms.

## Citation

If you use this code, please cite the article and the archived code:

- Hinkston, M. and Bradley, A. S. Quantifying biosynthetic constraints on molecular biomarker interpretation.
- Hinkston, M. and Bradley, A. S. (2026). Code and data for "Quantifying biosynthetic constraints on molecular biomarker interpretation", version 1.0.0. Zenodo. https://doi.org/10.5281/zenodo.23253971

Citation metadata are in `CITATION.cff`.
