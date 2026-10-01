"""
Check which "extra" analysis capabilities are actually supported by each
MalariaGEN vector dataset, as opposed to merely being present as an inherited
method.

Why this is needed
-------------------
Ag3, Af1, As1, Amin1, Adir1 and Adar1 all inherit from the same
`AnophelesDataResource` class (malariagen_data/anopheles.py), which mixes in
every analysis module: CNV (cnv_data.py / cnv_frq.py), ancestry-informative
markers / species calls (aim_data.py), phased haplotype data and everything
built on it (hap_data.py, hap_frq.py, hapclust.py, dipclust.py, h12.py,
h1x.py, xpehh.py), and insecticide-resistance phenotypes (phenotypes.py).

That means every dataset class exposes methods like `cnv_hmm()`,
`aim_variants()`, `haplotypes()`, `h12_gwss()`, `phenotype_data()`, etc,
*regardless* of whether the underlying data actually exists for that
dataset. Comparing the per-dataset "data downloads" pages on
https://malariagen.github.io/vector-data/ against the analysis IDs each
dataset class is actually configured with (malariagen_data/<dataset>.py)
turned up these cases worth checking manually:

  dataset   CNV data                         AIM / species calls   Haplotype / phased data
  -------   -------------------------------  --------------------  ------------------------
  Ag3       documented + configured          documented+configured documented + configured
  Af1       documented + configured          NOT configured        documented + configured
  As1       NOT configured                   NOT configured        NOT configured
  Amin1     configured with a placeholder    NOT configured        NOT configured
            analysis id "minimus_noneyet" -
            looks like a stub, not real data
  Adir1     configured in code, but NOT      NOT configured        configured in code, but NOT
            mentioned on the download page                         mentioned on the download page
  Adar1     configured in code, but NOT      NOT configured        configured in code, but NOT
            mentioned on the download page                         mentioned on the download page

None of the six datasets' download pages mention insecticide-resistance
phenotype data at all (phenotypes.py), so that is worth checking for every
dataset.

Affected methods, grouped by capability
----------------------------------------
CNV (malariagen_data/anoph/cnv_data.py, cnv_frq.py):
    cnv_hmm, cnv_coverage_calls, cnv_discordant_read_calls, gene_cnv,
    gene_cnv_frequencies, gene_cnv_frequencies_advanced,
    plot_cnv_hmm_coverage_track, plot_cnv_hmm_coverage,
    plot_cnv_hmm_heatmap_track, plot_cnv_hmm_heatmap

AIM / species calls (malariagen_data/anoph/aim_data.py):
    aim_variants, aim_calls, plot_aim_heatmap

Phased haplotype data and everything built on it
(hap_data.py, hap_frq.py, hapclust.py, dipclust.py, h12.py, h1x.py, xpehh.py):
    haplotypes, haplotype_sites, haplotypes_frequencies,
    haplotypes_frequencies_advanced, plot_haplotype_clustering,
    haplotype_pairwise_distances, plot_haplotype_clustering_advanced,
    transcript_haplotypes, cut_dist_tree, plot_haplotype_sharing_arc,
    plot_haplotype_sharing_chord, plot_diplotype_clustering,
    diplotype_pairwise_distances, plot_diplotype_clustering_advanced,
    h12_calibration, plot_h12_calibration, h12_gwss, plot_h12_gwss_track,
    plot_h12_gwss, plot_h12_gwss_multi_overlay_track,
    plot_h12_gwss_multi_overlay, plot_h12_gwss_multi_panel, h1x_gwss,
    plot_h1x_gwss_track, plot_h1x_gwss, xpehh_gwss, plot_xpehh_gwss_track,
    plot_xpehh_gwss, plot_haplotype_network (median-joining network)

Insecticide resistance phenotypes (malariagen_data/anoph/phenotypes.py):
    phenotype_data, phenotypes_with_snps, phenotypes_with_haplotypes,
    phenotype_sample_sets, phenotype_binary

NOT included above (SNP-based, works off data every dataset does have):
    pca, average_fst, pairwise_average_fst, g123_gwss,
    biallelic_snp_calls_ld_pruned - these only need SNP genotype data,
    which the download pages confirm is available for all six datasets, so
    they are not expected to be dataset-specific failures (worth spot
    checking anyway, but lower priority).

What this script does
----------------------
For each dataset, two checks are run:

1. Config check (no network access): inspect the private attributes each
   capability's methods rely on internally (`aim_ids`,
   `_default_phasing_analysis`, `_default_coverage_calls_analysis`,
   `_discordant_read_calls_analysis`) to see whether the capability is wired
   up at all. This is fast and gives a first-pass answer, but note some of
   these attributes (e.g. `_discordant_read_calls_analysis`) are only
   resolved from a config file loaded from GCS, so this step still needs a
   dataset instance (a small file read), just not any of the bulk data.

2. Live check (network access, reads a small amount of real data): calls
   the cheapest representative method for each capability and reports
   whether it returns data, returns "no data found", or raises an error.
   These are deliberately the lightest possible calls (opening a zarr
   store rather than loading/concatenating a whole dataset), but they do
   still talk to GCS, so this part is slower and requires connectivity.

Run this directly (`python check_dataset_capabilities.py`), or import
`check_dataset` and call it for a single dataset while you're poking around
interactively, e.g.:

    from check_dataset_capabilities import check_dataset
    import malariagen_data
    check_dataset("Amin1", malariagen_data.Amin1())
"""

from __future__ import annotations

import traceback
from dataclasses import dataclass, field
from typing import Any, Callable, Optional

import malariagen_data

# Same GCS access pattern already used in test_methods_for_docs.py. Adjust
# cache_storage per dataset if you want to keep the caches separate.
DATASET_FACTORIES: dict[str, Callable[[], Any]] = {
    "Ag3": lambda: malariagen_data.Ag3(
        "simplecache::gs://vo_agam_release_master_us_central1",
        simplecache=dict(cache_storage="gcs_cache_ag3"),
    ),
    "Af1": lambda: malariagen_data.Af1(
        "simplecache::gs://vo_afun_release_master_us_central1",
        simplecache=dict(cache_storage="gcs_cache_af1"),
    ),
    "As1": lambda: malariagen_data.As1(
        "simplecache::gs://vo_aste_release_master_us_central1",
        simplecache=dict(cache_storage="gcs_cache_as1"),
    ),
    "Amin1": lambda: malariagen_data.Amin1(
        "simplecache::gs://vo_amin_release_master_us_central1",
        simplecache=dict(cache_storage="gcs_cache_amin1"),
    ),
    "Adir1": lambda: malariagen_data.Adir1(
        "simplecache::gs://vo_adir_release_master_us_central1",
        simplecache=dict(cache_storage="gcs_cache_adir1"),
    ),
    "Adar1": lambda: malariagen_data.Adar1(
        "simplecache::gs://vo_adar_release_master_us_central1",
        simplecache=dict(cache_storage="gcs_cache_adar1"),
    ),
}


@dataclass
class CheckResult:
    capability: str
    status: str  # "configured" / "not configured" / "data found" / "no data found" / "error"
    detail: str
    affected_methods: tuple[str, ...] = field(default_factory=tuple)


CNV_METHODS = (
    "cnv_hmm",
    "cnv_coverage_calls",
    "cnv_discordant_read_calls",
    "gene_cnv",
    "gene_cnv_frequencies",
    "gene_cnv_frequencies_advanced",
    "plot_cnv_hmm_coverage_track",
    "plot_cnv_hmm_coverage",
    "plot_cnv_hmm_heatmap_track",
    "plot_cnv_hmm_heatmap",
)

AIM_METHODS = ("aim_variants", "aim_calls", "plot_aim_heatmap")

HAPLOTYPE_METHODS = (
    "haplotypes",
    "haplotype_sites",
    "haplotypes_frequencies",
    "haplotypes_frequencies_advanced",
    "plot_haplotype_clustering",
    "haplotype_pairwise_distances",
    "plot_haplotype_clustering_advanced",
    "transcript_haplotypes",
    "cut_dist_tree",
    "plot_haplotype_sharing_arc",
    "plot_haplotype_sharing_chord",
    "plot_diplotype_clustering",
    "diplotype_pairwise_distances",
    "plot_diplotype_clustering_advanced",
    "h12_calibration",
    "plot_h12_calibration",
    "h12_gwss",
    "plot_h12_gwss_track",
    "plot_h12_gwss",
    "plot_h12_gwss_multi_overlay_track",
    "plot_h12_gwss_multi_overlay",
    "plot_h12_gwss_multi_panel",
    "h1x_gwss",
    "plot_h1x_gwss_track",
    "plot_h1x_gwss",
    "xpehh_gwss",
    "plot_xpehh_gwss_track",
    "plot_xpehh_gwss",
    "plot_haplotype_network",
)

PHENOTYPE_METHODS = (
    "phenotype_data",
    "phenotypes_with_snps",
    "phenotypes_with_haplotypes",
    "phenotype_sample_sets",
    "phenotype_binary",
)


def _config_check(ds: Any) -> list[CheckResult]:
    results = []

    aim_ids = ds.aim_ids
    results.append(
        CheckResult(
            capability="AIM / species calls",
            status="configured" if aim_ids else "not configured",
            detail=f"aim_ids={aim_ids!r}",
            affected_methods=AIM_METHODS,
        )
    )

    phasing = getattr(ds, "_default_phasing_analysis", None)
    results.append(
        CheckResult(
            capability="Phased haplotype data",
            status="configured" if phasing else "not configured",
            detail=f"default_phasing_analysis={phasing!r}",
            affected_methods=HAPLOTYPE_METHODS,
        )
    )

    coverage = getattr(ds, "_default_coverage_calls_analysis", None)
    discordant = getattr(ds, "_discordant_read_calls_analysis", None)
    status = "configured" if (coverage or discordant) else "not configured"
    if coverage and "noneyet" in coverage:
        status = "configured (looks like a placeholder, not real data)"
    results.append(
        CheckResult(
            capability="CNV data",
            status=status,
            detail=(
                f"default_coverage_calls_analysis={coverage!r}, "
                f"discordant_read_calls_analysis={discordant!r}"
            ),
            affected_methods=CNV_METHODS,
        )
    )

    return results


def _live_check(ds: Any) -> list[CheckResult]:
    results = []

    # --- AIM / species calls -------------------------------------------
    aim_ids = ds.aim_ids
    if not aim_ids:
        results.append(
            CheckResult(
                "AIM / species calls",
                "not configured",
                "ds.aim_ids is empty, skipping live call",
                AIM_METHODS,
            )
        )
    else:
        try:
            aim_ds = ds.aim_variants(aims=aim_ids[0])
            n = aim_ds.sizes.get("variants", 0)
            results.append(
                CheckResult(
                    "AIM / species calls",
                    "data found" if n else "no data found",
                    f"aim_variants(aims={aim_ids[0]!r}) -> {n} variants",
                    AIM_METHODS,
                )
            )
        except Exception as exc:  # noqa: BLE001 - intentionally broad, this is a probe
            results.append(
                CheckResult(
                    "AIM / species calls",
                    "error",
                    f"aim_variants(aims={aim_ids[0]!r}) raised "
                    f"{exc.__class__.__name__}: {exc}",
                    AIM_METHODS,
                )
            )

    # --- Phased haplotype data ------------------------------------------
    try:
        root = ds.open_haplotype_sites()
        results.append(
            CheckResult(
                "Phased haplotype data",
                "data found",
                f"open_haplotype_sites() -> zarr group with keys {list(root.keys())[:5]!r}",
                HAPLOTYPE_METHODS,
            )
        )
    except Exception as exc:  # noqa: BLE001
        results.append(
            CheckResult(
                "Phased haplotype data",
                "error",
                f"open_haplotype_sites() raised {exc.__class__.__name__}: {exc}",
                HAPLOTYPE_METHODS,
            )
        )

    # --- CNV data ---------------------------------------------------------
    try:
        sample_sets_df = ds.sample_sets()
        first_sample_set = sample_sets_df["sample_set"].iloc[0]
        root = ds.open_cnv_hmm(sample_set=first_sample_set)
        if root is None:
            results.append(
                CheckResult(
                    "CNV data",
                    "no data found",
                    f"open_cnv_hmm(sample_set={first_sample_set!r}) -> None "
                    "(no CNV HMM zarr for this sample set)",
                    CNV_METHODS,
                )
            )
        else:
            results.append(
                CheckResult(
                    "CNV data",
                    "data found",
                    f"open_cnv_hmm(sample_set={first_sample_set!r}) -> zarr group found",
                    CNV_METHODS,
                )
            )
    except Exception as exc:  # noqa: BLE001
        results.append(
            CheckResult(
                "CNV data",
                "error",
                f"open_cnv_hmm(...) raised {exc.__class__.__name__}: {exc}",
                CNV_METHODS,
            )
        )

    # --- Insecticide resistance phenotypes --------------------------------
    try:
        sample_sets_with_phenotypes = ds.phenotype_sample_sets()
        results.append(
            CheckResult(
                "Insecticide resistance phenotypes",
                "data found" if sample_sets_with_phenotypes else "no data found",
                f"phenotype_sample_sets() -> {sample_sets_with_phenotypes!r}",
                PHENOTYPE_METHODS,
            )
        )
    except Exception as exc:  # noqa: BLE001
        results.append(
            CheckResult(
                "Insecticide resistance phenotypes",
                "error",
                f"phenotype_sample_sets() raised {exc.__class__.__name__}: {exc}",
                PHENOTYPE_METHODS,
            )
        )

    return results


def check_dataset(
    name: str, ds: Any, *, live: bool = True, verbose: bool = True
) -> dict[str, list[CheckResult]]:
    """
    Run the config check (and, unless live=False, the live check) for a
    single already-instantiated dataset, e.g.:

        import malariagen_data
        check_dataset("Amin1", malariagen_data.Amin1())
    """
    out = {"config": _config_check(ds)}
    if live:
        out["live"] = _live_check(ds)

    if verbose:
        _print_dataset_report(name, out)

    return out


def _print_dataset_report(name: str, report: dict[str, list[CheckResult]]) -> None:
    print(f"\n=== {name} ===")
    for stage, results in report.items():
        print(f"  [{stage}]")
        for r in results:
            print(f"    {r.capability:35s} {r.status:45s} {r.detail}")


def main(datasets: Optional[list[str]] = None, live: bool = True) -> dict:
    """
    Run the checks for all datasets (or a subset, e.g. main(["Amin1", "Adir1"])).
    Returns the full set of results in case you want to inspect them further
    or turn them into a table.
    """
    names = datasets or list(DATASET_FACTORIES.keys())
    all_results = {}

    for name in names:
        print(f"\nInstantiating {name}...")
        try:
            ds = DATASET_FACTORIES[name]()
        except Exception as exc:  # noqa: BLE001
            print(f"  FAILED to instantiate {name}: {exc.__class__.__name__}: {exc}")
            traceback.print_exc()
            continue

        all_results[name] = check_dataset(name, ds, live=live)

    return all_results


if __name__ == "__main__":
    main()
