PYTHON ?= python
RESULTS ?= results
WORKERS ?= 4
METACYC ?=
WORKBOOKS ?=

export PYTHONPATH := $(CURDIR)
export PYTHONHASHSEED := 2
export OMP_NUM_THREADS := 1
export OPENBLAS_NUM_THREADS := 1
export MKL_NUM_THREADS := 1
export NUMEXPR_NUM_THREADS := 1

NETWORKS = $(RESULTS)/networks
FIGURES = $(RESULTS)/figures

.PHONY: all networks pairs projected_paths tracing sources literature controls \
	supplementary_data figures numbers test check-inputs

all: networks pairs projected_paths tracing sources literature controls figures \
	supplementary_data numbers

check-inputs:
	@test -n "$(METACYC)" || (echo "Set METACYC to the MetaCyc 29.5 flat-file folder" && exit 1)

networks: check-inputs
	$(PYTHON) analyses/01_networks/build_networks.py --metacyc "$(METACYC)" --out "$(NETWORKS)"

pairs: check-inputs
	$(PYTHON) analyses/04_enzyme_comparisons/compare_pairs.py --metacyc "$(METACYC)" \
		--networks "$(NETWORKS)" --out "$(RESULTS)/pairs"

projected_paths: check-inputs
	$(PYTHON) analyses/02_projected_paths/projected_paths.py --metacyc "$(METACYC)" \
		--networks "$(NETWORKS)" --out "$(RESULTS)/projected_paths" --workers $(WORKERS)
	$(PYTHON) analyses/02_projected_paths/zero_path_causes.py \
		--paths "$(RESULTS)/projected_paths" --out "$(RESULTS)/projected_paths" --workers $(WORKERS)
	$(PYTHON) analyses/02_projected_paths/descriptor_cohort.py \
		--paths "$(RESULTS)/projected_paths" --out "$(RESULTS)/descriptors"
	$(PYTHON) analyses/02_projected_paths/descriptor_associations.py \
		--cohort "$(RESULTS)/descriptors" --out "$(RESULTS)/descriptors"
	$(PYTHON) analyses/09_figures/table1.py --results "$(RESULTS)" --out "$(RESULTS)/projected_paths"

tracing:
	$(PYTHON) analyses/03_upstream_tracing/trace_biomarkers.py --networks "$(NETWORKS)" \
		--pairs "$(RESULTS)/pairs" --out "$(RESULTS)/tracing" --workers $(WORKERS)
	$(PYTHON) analyses/03_upstream_tracing/formylhopane_reductase_test.py \
		--networks "$(NETWORKS)" --out "$(RESULTS)/tracing"

sources:
	@test -n "$(WORKBOOKS)" || (echo "Set WORKBOOKS to the folder of the Brocks et al. (2023) workbooks" && exit 1)
	$(PYTHON) analyses/05_source_compatibility/read_workbooks.py --workbooks "$(WORKBOOKS)" \
		--out "$(RESULTS)/sources"
	$(PYTHON) analyses/05_source_compatibility/source_tables.py --out "$(RESULTS)/sources"
	$(PYTHON) analyses/05_source_compatibility/compatible_combinations.py \
		--observations "$(RESULTS)/sources" --out "$(RESULTS)/sources"
	$(PYTHON) analyses/05_source_compatibility/sensitivity_tests.py \
		--results "$(RESULTS)/sources" --out "$(RESULTS)/sources"
	$(PYTHON) analyses/05_source_compatibility/select_examples.py --results "$(RESULTS)/sources"

literature: check-inputs
	$(PYTHON) analyses/06_literature_records/literature_reactions.py --networks "$(NETWORKS)" \
		--out "$(RESULTS)/literature"
	$(PYTHON) analyses/06_literature_records/metacyc_corrections.py --metacyc "$(METACYC)" \
		--out "$(RESULTS)/literature"

controls: check-inputs
	$(PYTHON) analyses/07_matched_controls/matched_controls.py --metacyc "$(METACYC)" \
		--networks "$(NETWORKS)" --out "$(RESULTS)/controls" --workers $(WORKERS)

figures: check-inputs
	$(PYTHON) analyses/09_figures/figure2.py --results "$(RESULTS)" --out "$(FIGURES)"
	$(PYTHON) analyses/09_figures/figure3.py --results "$(RESULTS)" --metacyc "$(METACYC)" --out "$(FIGURES)"
	$(PYTHON) analyses/09_figures/figure4.py --results "$(RESULTS)" --out "$(FIGURES)"
	$(PYTHON) analyses/09_figures/supplementary_figure1.py --results "$(RESULTS)" --out "$(FIGURES)"
	$(PYTHON) analyses/09_figures/supplementary_figure2.py --results "$(RESULTS)" --out "$(FIGURES)"
	$(PYTHON) analyses/09_figures/supplementary_figure3.py --results "$(RESULTS)" --out "$(FIGURES)"
	$(PYTHON) analyses/09_figures/supplementary_figure4.py --results "$(RESULTS)" --out "$(FIGURES)"
	$(PYTHON) analyses/09_figures/supplementary_figure5.py --results "$(RESULTS)" --out "$(FIGURES)"
	$(PYTHON) analyses/09_figures/supplementary_figure6.py --results "$(RESULTS)" --out "$(FIGURES)"

supplementary_data:
	$(PYTHON) analyses/08_supplementary_data/assemble.py --results "$(RESULTS)" \
		--out "$(RESULTS)/supplementary_data"

numbers: check-inputs
	$(PYTHON) analyses/02_projected_paths/text_numbers.py --results "$(RESULTS)"
	$(PYTHON) analyses/03_upstream_tracing/text_numbers.py --results "$(RESULTS)"
	$(PYTHON) analyses/04_enzyme_comparisons/text_numbers.py --results "$(RESULTS)" --metacyc "$(METACYC)"
	$(PYTHON) analyses/05_source_compatibility/text_numbers.py --results "$(RESULTS)/sources"
	$(PYTHON) analyses/06_literature_records/text_numbers.py --results "$(RESULTS)"
	$(PYTHON) analyses/07_matched_controls/text_numbers.py --results "$(RESULTS)"

test:
	$(PYTHON) -m pytest -q
