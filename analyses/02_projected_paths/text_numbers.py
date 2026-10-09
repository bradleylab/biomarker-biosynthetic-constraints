"""Print the projected-path and descriptor numbers quoted in the article and supplement.

Reads the outputs of projected_paths.py, zero_path_causes.py, descriptor_cohort.py
and descriptor_associations.py, and prints each number in a plain sentence, with
the precision used in the text.

  python analyses/02_projected_paths/text_numbers.py --results results
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from biosynthesis import descriptors as desc  # noqa: E402
from biosynthesis import simplified_graph as graph  # noqa: E402
from biosynthesis import statistics as stats  # noqa: E402
from biosynthesis.paths import (  # noqa: E402
    FIGURE_2B_EXAMPLE,
    SHARED_FRACTION,
    compute_dp_metrics,
    node_measurements,
    psi_bits,
    target_graph,
)
from biosynthesis.tables import read_csv  # noqa: E402

PRIMARY, CONDITIONAL = graph.NETWORK_VERSIONS
TREATMENT_NAMES = {
    "unscreened": "Unscreened input scores",
    "larger_search_values": "Expanded symmetry-search values only",
    "withhold_failures": "Withhold memory-limited and atomic-error scores",
    "withhold_all_warnings": "Withhold all diagnostic scores",
    "exclude_flagged_compounds": "Exclude all flagged compounds from every descriptor",
    "main_analysis": "Main analysis: expanded symmetry-search values; withhold unresolved scores",
}


def n(value: int) -> str:
    """An integer with thousands separators."""
    return f"{int(value):,}"


def f(value: float, digits: int) -> str:
    """A number at the given precision, with a true minus sign."""
    text = f"{value:.{digits}f}"
    return text.replace("-", "−")


def rows_of(path: Path, version: str | None = None) -> list[dict]:
    rows = read_csv(path)
    return [r for r in rows if version is None or r["network_version"] == version]


def path_counts(rows: list[dict]) -> Counter:
    return Counter(min(int(r["num_paths"]), 2) for r in rows)


def projection_numbers(results: Path) -> list[str]:
    paths_dir = results / "projected_paths"
    s2 = {
        v: rows_of(paths_dir / "TableS2_projection.csv", v)
        for v in graph.NETWORK_VERSIONS
    }
    s1 = {
        v: rows_of(paths_dir / "TableS1_projection.csv", v)
        for v in graph.NETWORK_VERSIONS
    }
    summaries = {
        (r["network_version"], r["population"]): r
        for r in rows_of(paths_dir / "TableS2_population_summaries.csv")
    }
    classes = {
        r["class"]: r
        for r in rows_of(paths_dir / "Table1_projection_class_means.csv", PRIMARY)
    }
    scope = {
        r["network_version"]: r
        for r in rows_of(paths_dir / "Projection_graph_and_comparison_scope.csv")
    }
    zero = json.loads((paths_dir / "zero_path_summary.json").read_text())
    causes = zero["zero_path_compounds_by_cause"]
    absent = zero["queried_compounds_not_in_graph_by_cause"]
    counts = path_counts(s2[PRIMARY])
    total = len(s2[PRIMARY])
    everything = summaries[(PRIMARY, "all")]
    limited = summaries[(PRIMARY, "score_at_most_1e6")]
    conditional_all = summaries[(CONDITIONAL, "all")]
    conditional_limited = summaries[(CONDITIONAL, "score_at_most_1e6")]
    by_id = {r["compound_id"]: r for r in s1[PRIMARY]}
    with_path = {v: sum(int(r["num_paths"]) > 0 for r in s1[v]) for v in s1}
    archaeal = [
        r
        for r in s1[PRIMARY]
        if r["class"] == "Archaeal lipid" and int(r["num_paths"]) > 0
    ]
    example = compute_dp_metrics(FIGURE_2B_EXAMPLE, "B")
    nodes = node_measurements(target_graph(FIGURE_2B_EXAMPLE, "B"), "B")["nodes"]
    shared = sum(node["q"] >= SHARED_FRACTION for node in nodes)
    single = [r for r in s2[PRIMARY] if int(r["num_paths"]) == 1]
    single_values = {(r["psi_bits"], r["lambda"], r["sigma"]) for r in single}
    main_classes = [c for c in classes.values() if int(c["n_with_path"]) >= 5]
    lines = [
        "Comparison set and path counts (primary graph):",
        f"  In the comparison set of {n(total)} compounds, {n(counts[0])} have no projected path, "
        f"{n(counts[1])} have one, and {n(counts[2])} have several.",
        f"  {n(causes['every_producing_record_lacks_an_eligible_reactant'])} of {n(counts[0])} "
        "compounds without a path lack one because every reactant of every producing record was "
        "excluded as a hub or carrier.",
        f"  Of the remaining {n(counts[0] - causes['every_producing_record_lacks_an_eligible_reactant'])}, "
        f"{causes['not_a_product']} were not recorded as products, {causes['no_producing_records']} had "
        f"no producing records and {causes['only_self_links_remain']} retained only self-links.",
        f"  For the {n(everything['multi_path_n'])} compounds with several paths, the medians of psi, "
        f"lambda and sigma are {f(float(everything['multi_path_psi_median']), 2)} bits "
        f"(about {2 ** float(everything['multi_path_psi_median']):.0f} paths), "
        f"{f(float(everything['multi_path_lambda_median']), 2)} and "
        f"{f(float(everything['multi_path_sigma_median']), 2)}.",
        f"  Pearson r = {f(float(everything['pearson_psi_lambda']), 2)} (psi, lambda), "
        f"{f(float(everything['pearson_psi_sigma']), 2)} (psi, sigma), "
        f"{f(float(everything['pearson_lambda_sigma']), 2)} (lambda, sigma); after excluding "
        f"cumulative scores above 10⁶, Pearson r (psi, lambda) = {f(float(limited['pearson_psi_lambda']), 2)}.",
        f"  Diploptene has {by_id['HOP-2229-ENE']['num_paths']} paths and cholesterol "
        f"{by_id['CHOLESTEROL']['num_paths']}.",
        f"  Every one of the {n(len(single))} single-path compounds has (psi, lambda, sigma) = "
        f"{', '.join(sorted(single_values)[0]) if len(single_values) == 1 else single_values}.",
        f"  After the cumulative-score filter (at most 10⁶), {n(limited['positive_path_n'])} compounds "
        f"with at least one path remained, including {n(limited['multi_path_n'])} with several paths, "
        f"whose medians were psi = {f(float(limited['multi_path_psi_median']), 2)} bits, "
        f"lambda = {f(float(limited['multi_path_lambda_median']), 3)} and "
        f"sigma = {f(float(limited['multi_path_sigma_median']), 2)}.",
        f"  With conditional connections, {n(conditional_all['positive_path_n'])} compounds had at least "
        f"one path and {n(conditional_all['multi_path_n'])} had several; after filtering, "
        f"{n(conditional_limited['positive_path_n'])} and {n(conditional_limited['multi_path_n'])} remained.",
        f"  Few paths branching far upstream (psi ≤ {stats.SCREEN_PSI_MAX}, lambda > "
        f"{stats.SCREEN_LAMBDA_MIN}, sigma > {stats.SCREEN_SIGMA_MIN}): {everything['few_paths_far_branching_n']} "
        "compounds with several paths (primary graph, all).",
        "",
        "Graph scope:",
        f"  {total - int(scope[PRIMARY]['queried_compounds_n'])} compounds of the comparison set share "
        f"a query with another, giving {n(scope[PRIMARY]['queried_compounds_n'])} distinct queries for the "
        f"{n(scope[PRIMARY]['comparison_compounds_n'])} compounds.",
        f"  The primary graph contains {n(scope[PRIMARY]['graph_nodes_n'])} nodes and the graph with "
        f"conditional connections {n(scope[CONDITIONAL]['graph_nodes_n'])}.",
        f"  Of the {n(total)} compounds, {scope[PRIMARY]['queried_compounds_not_in_graph_n']} do not appear "
        f"in the primary graph (by query): "
        f"{absent.get('every_producing_record_lacks_an_eligible_reactant', 0)} are produced only by records "
        f"whose reactants are all excluded, and "
        f"{absent.get('not_a_product', 0) + absent.get('no_producing_records', 0)} have no producing "
        f"reaction ({absent.get('not_a_product', 0)} not recorded as products, "
        f"{absent.get('no_producing_records', 0)} recorded only as their own reactant).",
        "",
        "Hub exemptions:",
        f"  Without the two exemptions, {n(zero['comparison_compounds_fewer_paths'])} of the {n(total)} "
        f"compounds had fewer paths, {zero['comparison_compounds_more_paths']} had more, and "
        f"{zero['comparison_compounds_losing_every_path']} lost every path.",
        f"  Among the {zero['biomarkers']} biomarkers, {zero['biomarkers_fewer_paths']} had fewer paths, "
        f"including {zero['biomarkers_fewer_paths_by_class'].get('Sterol', 0)} of the "
        f"{zero['biomarker_class_sizes']['Sterol']} sterols.",
        f"  {', '.join(zero['biomarkers_losing_every_path'])} each had "
        f"{'/'.join(map(str, zero['biomarkers_losing_every_path_counts']))} paths and lost them all "
        f"({len(zero['biomarkers_losing_every_path'])} biomarkers).",
        "",
        "Biomarkers (primary graph; Table 1):",
        f"  Of the {len(s1[PRIMARY])} biomarkers, {with_path[PRIMARY]} had at least one projected path; "
        f"with conditional connections, {with_path[CONDITIONAL]}.",
        f"  Classes with at least five members with a path: {len(main_classes)} "
        f"({', '.join(c['class'] for c in main_classes)}).",
    ]
    for c in main_classes:
        lines.append(
            f"  {c['class']}: mean psi = {f(float(c['mean_psi']), 2)} bits, mean lambda = "
            f"{f(float(c['mean_lambda']), 2)}, mean sigma = {f(float(c['mean_sigma']), 2)}."
        )
    for name in ("Hopanoid", "Archaeal lipid", "Algal hydrocarbon", "Fatty acid"):
        lines.append(
            f"  {name}: {classes[name]['n_with_path']} of {classes[name]['n']} have a path."
        )
    lines += [
        f"  The archaeal lipids with a path have at most {max(int(r['num_paths']) for r in archaeal)} path.",
        "",
        "Measurement definitions and the worked example of Figure 2B:",
        f"  Hub limit {graph.HUB_LIMIT} product-specific records; {len(graph.HUB_EXEMPTIONS)} exemptions "
        f"({', '.join(sorted(graph.HUB_EXEMPTIONS))}); {len(graph.TRANSFER_AND_REDOX_CARRIERS)} transfer "
        "and redox carriers.",
        f"  Shared nodes are on at least {SHARED_FRACTION:.0%} of paths ({SHARED_FRACTION:.2f}).",
        f"  Worked example: {example['num_paths']:.0f} paths, psi = {psi_bits(example['num_paths']):.0f} bit, "
        f"lambda = {example['mean_branching_depth_normalized']:.1f}, {shared} nodes on both paths and "
        f"{len(nodes) - shared} on one, sigma = {shared}/({shared} + 0.5) = "
        f"{example['fraction_shared_80']:.2f}.",
        f"  Population filter: cumulative branching score at most 10^{round(math.log10(stats.SCORE_LIMIT))}; "
        "class SDs use the n − 1 denominator and need at least two biomarkers with a path.",
    ]
    return lines


def table1_lines(results: Path) -> list[str]:
    table = results / "projected_paths/Table1.md"
    if not table.exists():
        return []
    return ["", "Table 1 (Table1.md):"] + [
        "  " + line
        for line in table.read_text(encoding="utf-8").splitlines()
        if line.startswith("| ")
    ]


def descriptor_numbers(results: Path) -> list[str]:
    folder = results / "descriptors"
    selection = read_csv(folder / "TableS4_descriptor_selection.csv")
    cohort = {
        v: rows_of(folder / "TableS4_descriptor_cohort.csv", v)
        for v in graph.NETWORK_VERSIONS
    }
    associations = rows_of(folder / "TableS4_associations.csv")
    treatments = rows_of(folder / "TableS4_score_treatment_summary.csv")

    def fit(version: str, cohort_name: str) -> dict:
        return next(
            r
            for r in associations
            if r["network_version"] == version
            and r["cohort"] == cohort_name
            and r["analysis_type"] == "ols_standardized"
        )

    def strongest(version: str, cohort_name: str) -> dict:
        rows = [
            {**r, "spearman_rho": float(r["spearman_rho"])}
            for r in associations
            if r["network_version"] == version
            and r["cohort"] == cohort_name
            and r["analysis_type"] == "spearman"
        ]
        return stats.strongest_psi_correlation(rows)

    def panel_n(version: str, descriptor: str) -> str:
        return next(
            r["n"]
            for r in associations
            if r["network_version"] == version
            and r["cohort"] == "all"
            and r["target"] == "psi_bits"
            and r["descriptor"] == descriptor
        )

    primary = cohort[PRIMARY]
    status = Counter(r["bottcher_status"] for r in primary)
    memory = [
        r for r in primary if "automorphism_memory_limit" in r["bottcher_warnings"]
    ]
    flagged = sum(r["bottcher_warnings"] != "" for r in primary)
    fixed, eligible = fit(PRIMARY, "all"), fit(PRIMARY, "structurally_eligible")
    cfixed, celigible = (
        fit(CONDITIONAL, "all"),
        fit(CONDITIONAL, "structurally_eligible"),
    )
    top, top_eligible, ctop = (
        strongest(PRIMARY, "all"),
        strongest(PRIMARY, "structurally_eligible"),
        strongest(CONDITIONAL, "all"),
    )
    other_n = {panel_n(PRIMARY, d) for d in stats.FEATURES if d != "bottchscore"}
    cother_n = {panel_n(CONDITIONAL, d) for d in stats.FEATURES if d != "bottchscore"}
    lines = [
        "",
        "Descriptor cohort (Table S4):",
        f"  {len(stats.FEATURES)} descriptors; {len(stats.FEATURES) - 1} calculated with RDKit "
        f"{desc.RDKIT_DESCRIPTOR_VERSION}; structural eligibility with RDKit {desc.RDKIT_ELIGIBILITY_VERSION}.",
        f"  {n(len(selection))} compounds with descriptor records; "
        f"{sum(r['has_path_reference_graph'] == 'false' for r in selection)} without a path in the reference "
        f"graph; {sum(float(r['cumulative_score_reference_graph']) > desc.REFERENCE_SCORE_LIMIT for r in selection)} "
        f"above the cumulative-score limit of 10⁶; "
        f"{sum(r['compound_id'] in desc.CURRENCY_COMPOUNDS for r in selection)} currency compounds; "
        f"{n(sum(r['selected'] == 'true' for r in selection))} selected.",
        f"  Without a path: {sum(int(r['num_paths']) == 0 for r in primary)} cohort compounds in the primary "
        f"graph and {sum(int(r['num_paths']) == 0 for r in cohort[CONDITIONAL])} with conditional connections.",
        f"  Strongest psi correlation: {top['descriptor']} (rho = {f(top['spearman_rho'], 3)}); "
        f"R² = {float(fixed['r_squared']):.1%} = {f(float(fixed['r_squared']), 3)} (adjusted "
        f"{f(float(fixed['adjusted_r_squared']), 3)}) for {n(fixed['complete_case_n'])} complete cases.",
        f"  Structurally eligible: {len(primary) - int(eligible['analysis_rows'])} removed, leaving "
        f"{n(eligible['analysis_rows'])} compounds and {n(eligible['complete_case_n'])} complete cases; "
        f"largest correlation {f(abs(top_eligible['spearman_rho']), 3)}; R² {f(float(eligible['r_squared']), 3)} "
        f"(adjusted {f(float(eligible['adjusted_r_squared']), 3)}).",
        f"  With conditional connections: strongest correlation {ctop['descriptor']} (rho = "
        f"{f(ctop['spearman_rho'], 3)}); R² = {f(float(cfixed['r_squared']), 3)} and "
        f"{f(float(celigible['r_squared']), 3)} (adjusted {f(float(cfixed['adjusted_r_squared']), 3)} and "
        f"{f(float(celigible['adjusted_r_squared']), 3)}) for {n(cfixed['complete_case_n'])} and "
        f"{n(celigible['complete_case_n'])} complete cases.",
        f"  Supplementary Figure 1 denominators: Böttcher n = {n(panel_n(PRIMARY, 'bottchscore'))}, other "
        f"panels n = {', '.join(n(x) for x in sorted(other_n))}; with conditional connections "
        f"{n(panel_n(CONDITIONAL, 'bottchscore'))} and {', '.join(n(x) for x in sorted(cother_n))}.",
        "",
        "Böttcher scores:",
        f"  Settings: disable_mesomer={desc.BOTTCHER_SETTINGS['disable_mesomer']}, first symmetry-search "
        f"limit {n(desc.BOTTCHER_SETTINGS['automorp_memory_maxsize'])}, larger limit "
        f"{n(desc.LARGER_SYMMETRY_SEARCH)}; openbabel-wheel {desc.OPENBABEL_WHEEL_VERSION} (library "
        f"{desc.OPENBABEL_LIBRARY_VERSION}).",
        f"  {len(memory)} compounds repeated with the larger limit; {status[desc.RERUN]} had no warning and "
        f"{sum(r['bottcher_status'] == desc.WITHHELD for r in memory)} still had a stereochemical warning.",
        f"  {status[desc.WITHHELD]} scores withheld; {status[desc.NO_SCORE]} compounds without an input "
        f"score; {flagged} flagged compounds.",
        "",
        "Böttcher score treatments (primary graph, all cohort compounds):",
    ]
    for r in treatments:
        if r["network_version"] == PRIMARY and r["cohort"] == "all":
            lines.append(
                f"  {TREATMENT_NAMES[r['score_treatment']]}: {n(r['regression_complete_case_n'])} complete "
                f"cases, R² = {f(float(r['r_squared']), 4)}, largest absolute psi correlation "
                f"{f(float(r['largest_absolute_psi_rho']), 4)} ({r['strongest_psi_descriptor']})."
            )
    return lines


def main(results: Path) -> None:
    print(
        "\n".join(
            projection_numbers(results)
            + table1_lines(results)
            + descriptor_numbers(results)
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--results", type=Path, required=True)
    args = parser.parse_args()
    main(args.results)
