# check output of to_vcf
import malariagen_data

def check_wgs_data_catalog(dataset)-> tuple:
    sample_sets = dataset.sample_sets()["sample_set"].to_list() 
    accessions = {}
    samples_no_data = {}
    for sample_set in sample_sets:
        try:
            result = dataset.wgs_data_catalog(sample_set)
            accessions[sample_set] = result
        except ValueError:
            samples_no_data.append(sample_set)
    return accessions, samples_no_data

ag3 = malariagen_data.Ag3(
    "simplecache::gs://vo_agam_release_master_us_central1",
    simplecache=dict(cache_storage="/Users/katie.barr/malariagen-data-python/gcs_cache_ag3"),
    results_cache="/Users/katie.barr/malariagen-data-python/results_cache",
)

res = check_wgs_data_catalog(ag3)
sample_set_data = res[0]['bergey-2019']

sample_id = sample_set_data["sample_id"].iloc[0]
vcf_fileneme = sample_set_data["snp_genotypes_vcf"].iloc[0]

# bgzip -f eva_variants.vcf
# tabix -p vcf eva_variants.vcf.gz
# output_path: vcf_params.vcf_output_path,
#         region: base_params.regions,
#         sample_sets: Optional[base_params.sample_sets] = None,
#         sample_query: Optional[base_params.sample_query] = None,
#         sample_query_options: Optional[base_params.sample_query_options] = None,
#         sample_indices: Optional[base_params.sample_indices] = None,
#         site_mask: Optional[base_params.site_mask] = base_params.DEFAULT,
#         inline_array: base_params.inline_array = base_params.inline_array_default,
#         chunks: base_params.chunks = base_params.native_chunks,
#         overwrite: plink_params.overwrite = False,
#         fields: vcf_params.vcf_fields = ("GT",),


# this fails eventually, traceback saved for investigation
ag3.snp_calls_to_vcf(output_path="/Users/katie.barr/malariagen_api_utils/eva/test_vcf.vcf", 
                     sample_sets = ['bergey-2019'],
                     sample_query = f"sample_id == '{sample_id}'",
                     region=ag3.contigs,
                     fields = {"GT", "GQ", "AD", "MQ"})


# retry with 1 chrom
ag3.snp_calls_to_vcf(output_path="/Users/katie.barr/malariagen_api_utils/eva/test_vcf_chrom_2R.vcf", 
                     sample_sets = ['bergey-2019'],
                     sample_query = f"sample_id == '{sample_id}'",
                     region="2R",
                     fields = {"GT", "GQ", "AD", "MQ"})

# runs in a couple of minutes

# test claude changes re batching
import time

start = time.time()
ag3.snp_calls_to_vcf(output_path="/Users/katie.barr/malariagen_api_utils/eva/test_vcf_chrom_2R_claude_changes_numpy.vcf", 
                     sample_sets = ['bergey-2019'],
                     sample_query = f"sample_id == '{sample_id}'",
                     region="2R",
                     fields = {"GT", "GQ", "AD", "MQ"},
                     overwrite = True)
end = time.time()
total_time = end-start
print(f"Total time taken is: {total_time}") 
# Total time taken is: 103.31019186973572
# after numpy edit
# Total time taken is: 93.28258991241455


# claude re numpy 
# Fixed. The diff is small and surgical: `sample_fields = np.empty(n_samples, dtype=object)` + indexed assignment (`sample_fields[k] = ...`) is now `sample_fields = []` + `sample_fields.append(...)`. `"\t".join(sample_fields)` afterward is unaffected — it works identically on a list or a numpy object array, so nothing downstream needed to change.

# Why this was the right call: that array held short Python strings, built one at a time via ordinary string concatenation — there's no vectorized/numeric operation happening on it at all, so numpy bought nothing but its own allocation overhead, paid on every single variant row (up to ~150M times for a whole-genome export). A quick microbenchmark confirms the per-call overhead: `np.empty(1, dtype=object)` + assignment vs a plain list append is about **6.6x slower** per iteration for this exact pattern. That's not the dominant cost by itself at 143M rows (order of a few tens of seconds, not hours), but it's real, free to fix, and exactly the kind of thing that adds up alongside the I/O batching fix from earlier.

# 20/20 tests pass, ruff and mypy clean. The two performance issues I flagged from the original 13-hour run are now both addressed; the third thing I mentioned — whether a 13-hour whole-genome single-shot export is the right shape at all for an IGV-feeding exporter — is still an open design question whenever you want to revisit it.
## output with None for sample_sets and sample_query looks literally nothing like a VCF if sample sets and query are none, don't think they should be able to be none, makes no sense for VCF

# to compare
# bcftools isec -p output_dir SAMN12920115.vcf.gz test_vcf.vcf.gz
# bcftools stats SAMN12920115.vcf.gz test_vcf.vcf.gz > comparison_stats.txt

# No quality on output VCF
# format string varies- not sure that matters

# so far (16.49) been writing for nearly 40 minutes

# Claude
# 2. Suspected AD / Number=R mismatch — likely, worth confirming against real data.

# variant_allele always has 4 slots (REF + 3 ALT, padded with empty string), and the code correctly drops empty ALT slots when building the ALT column (if s: alt_alleles.append(s)). But call_AD is also always 4 slots, and the AD-writing code (",".join(... for x in ad_vals), to_vcf.py:254-258) emits all 4 values regardless of how many ALT alleles were actually written. So for an ordinary biallelic SNP, you'd get ALT=T but AD=12,8,-1,-1 — 4 values against a header declaring Number=R (which for 1 REF + 1 ALT should be exactly 2). That's a real spec violation if it plays out the way I think, and something like bcftools view or vcf-validator should catch it immediately. I haven't run it against real data to confirm, so treat this as a strong hypothesis, not a confirmed bug.

# 3. Silent partial failure on optional fields.

# gq_chunk/ad_chunk/mq_chunk are wrapped in try/except (FileNotFoundError, KeyError): pass per chunk (to_vcf.py:169-183). If a field fails to load for one chunk but not others (patchy data, a transient GCS hiccup, a genuinely missing sub-array for part of a region), that chunk silently gets . for every sample in every requested extra field — no warning, no error. The output VCF looks complete and well-formed but has quietly dropped real data in a way a casual read-through wouldn't catch. Given you're specifically worried about ENA/EVA reachability, this pattern is worth a second look — it's exactly the kind of thing that turns "URL flaky" into "wrong answer with no error," which is worse than a loud failure.

# 4. Cosmetic but real fidelity gaps vs. a genuine VCF.

# QUAL, FILTER, ID, and INFO are unconditionally .. QUAL/ID/INFO absence is expected (this data isn't in the Zarr arrays), so that's a documented limitation rather than a bug. FILTER=. is more debatable: if site_mask has already been applied as a hard subset filter upstream (removing failing sites before they ever reach the writer), then every row here technically passed, and a stricter tool would write PASS, not . — . means "filter not evaluated," which isn't quite true. Not wrong, but worth deciding deliberately rather than by default.

# Core GT logic itself looks sound. Missing-genotype detection ((a0<0)|(a1<0) → ./.), unphased "/" joining, and the vectorized string-building all look internally consistent, and the allele-index-to-REF/ALT mapping follows the same convention used elsewhere in the codebase (index 0 = REF, GT integers index into variant_allele). I'd put this at "probably fine" rather than "verified," since it's exactly the kind of thing a real ground-truth diff would either confirm or immediately blow up.

# How I'd validate against ground truth:

# The cleanest ground truth you have access to is the actual per-sample VCF this data was derived from — the same one referenced by wgs_data_catalog()["snp_genotypes_vcf"], which igv.py already points IGV at. That means you're comparing the exporter against the file it's meant to be a fallback for, which is exactly the right comparison.

# Pick one real sample and a small but non-trivial region (a few kb, ideally one that contains at least one multiallelic site, to specifically stress-test the AD/ALT-count question above).
# Pull the ground truth slice: bcftools view -r <region> -s <sample> <snp_genotypes_vcf_url> > truth.vcf.
# Generate the same slice from malariagen_data: api.snp_calls_to_vcf(output_path=..., region=<region>, sample_sets=..., sample_query=f"sample_id == '{sample}'", site_mask=None, fields=("GT","GQ","AD","MQ")) — use site_mask=None for this first pass so you're comparing raw calls, not filtered ones; validate the site-mask filtering as a separate step afterward.
# Don't diff the raw text — normalize first (bcftools norm, or just load both into scikit-allel/pandas and compare arrays). Compare, per site: POS/REF/ALT set equality, then GT per sample, then GQ/AD/MQ per sample if included.
# Specifically assert len(AD_values) == 1 + len(ALT_alleles) per row — that directly tests point 2 above.
# Compare row counts over the whole region too, not just the rows that do align — a silently-dropped chunk (point 3) would show up as missing rows rather than wrong values in a naive row-by-row diff, so check that separately.
# Once GT/AD/GQ/MQ match, layer site_mask back in and confirm the included/excluded site set matches the corresponding variant_filter_pass_* array from snp_calls() directly, independent of the VCF round-trip.
# I'd fix #1 first (it blocks everything else from even running), then use that harness to settle #2 empirically before trusting the exporter for anything beyond the single-sample case it's now scoped to.

# I stopped the sample_Sets and sample_query from being optional because VCF as a format makes no sense if you don't specify what you are calling variants based on. Leaving those as "None" resulted in nonesense VCF output
# Show more
# That reasoning holds up — a VCF without a defined cohort is genuinely ambiguous (which samples were the ALT alleles even called against?), so requiring sample_sets makes sense. Two things are still worth separating out from that design decision, though:

# The test suite is now uniformly broken regardless of whether the decision was right. tests/anoph/test_vcf_exporter.py predates both commits and was never updated to pass sample_query, so right now there's zero passing test coverage on this exporter — not "some edge cases fail," literally all 10 tests raise TypeError before touching data. That needs fixing either way, independent of whether required args was the correct call.

# The comment says something narrower than "you must specify a cohort." It says sample_query "must be of form f\"sample_id == '{sample_id}'\"" — i.e. single-sample only, not just non-None. That's a different and stronger constraint than what you just described. Is the intent that this exporter is now scoped to single-sample export only (which would make sense if it's meant to feed per-sample IGV tracks), or is any valid multi-sample query still fine and the comment is just documenting the one pattern you'd actually tested? That distinction matters for the AD-padding concern I raised — it's relevant either way, but a hard single-sample restriction would also make sense of "nonsense output": if you'd tried a multi-sample export and hit something broken (rather than just conceptually undefined), that's a different, more specific bug than "None was accepted" — worth telling me what the actual nonsense looked like if you still have it, since that might point at something more concrete than the ambiguity argument alone.

# crashes eventually

# check output of to_vcf
import malariagen_data

def check_wgs_data_catalog(dataset)-> tuple:
    sample_sets = dataset.sample_sets()["sample_set"].to_list() 
    accessions = {}
    samples_no_data = {}
    for sample_set in sample_sets:
        try:
            result = dataset.wgs_data_catalog(sample_set)
            accessions[sample_set] = result
        except ValueError:
            samples_no_data.append(sample_set)
    return accessions, samples_no_data

ag3 = malariagen_data.Ag3(
    "simplecache::gs://vo_agam_release_master_us_central1",
    simplecache=dict(cache_storage="/Users/katie.barr/malariagen-data-python/gcs_cache_ag3"),
    results_cache="/Users/katie.barr/malariagen-data-python/results_cache",
)

res = check_wgs_data_catalog(ag3)
sample_set_data = res[0]['bergey-2019']

sample_id = sample_set_data["sample_id"].iloc[0]
vcf_fileneme = sample_set_data["snp_genotypes_vcf"].iloc[0]

import time

start = time.time()
ag3.snp_calls_to_vcf(output_path="/Users/katie.barr/malariagen_api_utils/eva/test_vcf_chrom_2R_no_pysam.vcf.gz", 
                     sample_sets = ['bergey-2019'],
                     sample_query = f"sample_id == '{sample_id}'",
                     region="2R",
                     fields = {"GT", "GQ", "AD", "MQ"},
                     overwrite = True,
                     site_mask=None)
end = time.time()
total_time = end-start
print(f"Total time taken is: {total_time}") 

start = time.time()
ag3.snp_calls_to_vcf(output_path="/Users/katie.barr/malariagen_api_utils/eva/test_vcf_chrom_2L_pysam.vcf.gz", 
                     sample_sets = ['bergey-2019'],
                     sample_query = f"sample_id == '{sample_id}'",
                     region="2L",
                     fields = {"GT", "GQ", "AD", "MQ"},
                     overwrite = True,
                     site_mask=None)
end = time.time()
total_time = end-start
print(f"Total time taken is: {total_time}") 
# pre edit: Total time taken is: 140.784765958786


ag3.snp_calls_to_vcf(output_path="/Users/katie.barr/malariagen_api_utils/eva/test_vcf_chrom_2R_restarting_exclude_0_0.vcf.gz", 
                     sample_sets = ['bergey-2019'],
                     sample_query = f"sample_id == '{sample_id}'",
                     region="2R",
                     site_mask = None,
                     fields = {"GT", "GQ", "AD", "MQ"},
                     non_ref_only = True,
                     overwrite = True)


start = time.time()
ag3.snp_calls_to_vcf(output_path="/Users/katie.barr/malariagen_api_utils/eva/test_vcf_all_pysam_refactored.vcf.gz", 
                     sample_sets = ['bergey-2019'],
                     sample_query = f"sample_id == '{sample_id}'",
                     region=ag3.contigs,
                     fields = {"GT", "GQ", "AD", "MQ"},
                     overwrite = True,
                     site_mask=None)
end = time.time()
total_time = end-start
print(f"Total time taken is: {total_time}") 


from malariagen_data.anoph import base_params
import numpy as np 

sno_calls = ag3.snp_calls(
            region="2R",
            sample_sets=['bergey-2019'],
            sample_query=f"sample_id == '{sample_id}'",
            site_mask=None,
            inline_array=base_params.inline_array_default,
            chunks=base_params.native_chunks,
        )

gt_data = sno_calls["call_genotype"].data
pos_data = sno_calls["variant_position"].data
contig_data = sno_calls["variant_contig"].data
allele_data = sno_calls["variant_allele"].data

chunk_sizes = gt_data.chunks[0] # tuple of size of each chunk
offsets = np.cumsum((0,) + chunk_sizes) # np array of amount to iterate through each time

from dataclasses import dataclass
from dask.array.core import Array
from collections import OrderedDict
import dask

@dataclass
class VariantChunkData:
    gq_chunk: Array
    ad_chunk: Array
    mq_chunk: Array
    gt_chunk: Array
    pos_chunk: Array
    contig_chunk: Array
    allele_chunk: Array

gq_data = sno_calls["call_GQ"].data

ad_data = sno_calls["call_AD"].data
mq_data = sno_calls["call_MQ"].data

optional_arrays = OrderedDict()
optional_arrays["GQ"] = gq_data
optional_arrays["AD"] = ad_data
optional_arrays["MQ"] = mq_data

ci=5
start = offsets[ci]
stop = offsets[ci + 1]
gt_chunk, pos_chunk, contig_chunk, allele_chunk = dask.compute(
    gt_data[start:stop],
    pos_data[start:stop],
    contig_data[start:stop],
    allele_data[start:stop],
)


computed = dask.compute(
    *(arr[start:stop] for arr in optional_arrays.values())
)
optional_chunks = dict(zip(optional_arrays.keys(), computed))
gq_chunk = optional_chunks.get("GQ")
ad_chunk = optional_chunks.get("AD")
mq_chunk = optional_chunks.get("MQ")

gt_chunk_2d = gt_chunk.reshape( # for single sample this makes it the same shape
    gt_chunk.shape[0],
    gt_chunk.shape[1],
    2,
)
a0 = gt_chunk_2d[:, :, 0]  # (n_variants, n_samples)
a1 = gt_chunk_2d[:, :, 1]  # (n_variants, n_samples)
missing = (a0 < 0) | (a1 < 0) # array of loci not present

# Build formatted GT strings using NumPy vectorization
gt_formatted = np.empty(
    (gt_chunk.shape[0], 1), dtype=object
)
gt_formatted[missing] = "./."
present_idx = ~missing
if np.any(present_idx):
    a0_str = a0[present_idx].astype(str)
    a1_str = a1[present_idx].astype(str)
    gt_formatted[present_idx] = np.char.add(
        np.char.add(a0_str, "/"), a1_str
    )

# then this is done for each position in chunk
contigs = sno_calls.attrs.get("contigs", ag3.contigs)

j = 12
chrom = contigs[contig_chunk[j]]
pos = str(pos_chunk[j])
alleles = allele_chunk[j]
ref = (
    alleles[0].decode()
    if hasattr(alleles[0], "decode")
    else str(alleles[0])
)

alt_alleles = []
for a in alleles[1:]:
    s = a.decode() if hasattr(a, "decode") else str(a)
    if s:
        alt_alleles.append(s)
alt = ",".join(alt_alleles) if alt_alleles else "."

parts = [gt_formatted[j, 0]]

gt = gt_formatted[:, 0]

not_0_0 = [i for i in list(gt) if i != '0/0']

# format string
format_str = "GT:GQ:AD:MQ"

import pysam
import pysam.bcftools
pysam.bcftools.index("--csi", "ex2.vcf.gz")


