# need to compare vcf output to vcf downloaed based on accession
bgzip SAMN12920115.vcf -o SAMN12920115.vcf.gz

bcftools index SAMN12920115.vcf.gz
bcftools view -r "2R" SAMN12920115.vcf.gz -u -o SAMN12920115_2R.vcf
bcftools view -r "2R" SAMN12920115.vcf.gz -o SAMN12920115_2R_no_u_passed_to_bcftools_view.vcf


bgzip SAMN12920115_2R.vcf -o SAMN12920115_2R.vcf.gz

bcftools index SAMN12920115_2R.vcf.gz


bcftools isec -p compare_vcfs SAMN12920115.vcf.gz test_vcf.vcf.gz

# https://davetang.org/muse/2019/09/02/comparing-vcf-files/
bedtools intersect -u -a SAMN12920115.vcf.gz -b test_vcf.vcf.gz | wc -l

 
# calculate Jaccard index
bedtools jaccard -a SAMN12920115.vcf.gz -b test_vcf.vcf.gz

SnpSift concordance -v SAMN12920115.vcf test_vcf.vcf > compare_vcfs/snp_concordance.txt

# now have 3 vcfs, one based on original code, 2 based on claude changes
md5sum test_vcf_chrom_2R_claude_changes.vcf
#2bdc280d45a3c54f684afccbc027d8e2  test_vcf_chrom_2R_claude_changes.vcf
md5sum test_vcf_chrom_2R_claude_changes_numpy.vcf
#3d6e0c3a3e06c58073f0a89807cc0110  test_vcf_chrom_2R_claude_changes_numpy.vcf
md5sum test_vcf_chrom_2R.vcf 
#9a196524aa8fca0db076ce7dea8ebef1  test_vcf_chrom_2R.vcf

wc -l test_vcf_chrom_2R_claude_changes.vcf
# 40561680 test_vcf_chrom_2R_claude_changes.vcf
wc -l test_vcf_chrom_2R_claude_changes_numpy.vcf
# 40561680 test_vcf_chrom_2R_claude_changes_numpy.vcf
wc -l test_vcf_chrom_2R.vcf
# 40561680 test_vcf_chrom_2R.vcf

wc -l SAMN12920115_2R.vcf
# 5102135 SAMN12920115_2R.vcf
# so 10x more variants

bgzip test_vcf_chrom_2R_claude_changes.vcf -o test_vcf_chrom_2R_claude_changes.vcf.gz
bcftools index test_vcf_chrom_2R_claude_changes.vcf.gz

bgzip test_vcf_chrom_2R_claude_changes_numpy.vcf -o test_vcf_chrom_2R_claude_changes_numpy.vcf.gz
bcftools index test_vcf_chrom_2R_claude_changes_numpy.vcf.gz

bgzip test_vcf_chrom_2R.vcf -o test_vcf_chrom_2R.vcf.gz
bcftools index test_vcf_chrom_2R.vcf.gz


# I pulled up that WDL file to check exactly what the pipeline does, and it actually points fairly clearly at the explanation.

# The `UnifiedGenotyper` task runs with `--alleles ~{alleles_vcf}` and `--output_mode EMIT_ALL_SITES`. That means:

# - It genotypes against a **fixed, predetermined site list** (`alleles_vcf`) shared across the whole cohort — not sites discovered per-sample. That's the same kind of joint/cohort-wide site list that backs `snp_calls()`: every sample set's `call_genotype` array is aligned to one shared genome-wide sites list, not a per-sample one.
# - `EMIT_ALL_SITES` means the raw task output, if that's what you're comparing against, should in principle contain a record at *every* site in that shared list for the sample — reference calls (`0/0`) included, same as `snp_calls_to_vcf` does.

# So if your "reference VCF" is that raw, unfiltered per-task output, the site counts should be close to matching (modulo `site_mask`, more below). But if what you downloaded is the **released, distributed per-sample VCF** (e.g. via `wgs_data_catalog`'s `snp_genotypes_vcf` column) rather than the raw task output, it's very likely been through a later pipeline step that strips homozygous-reference records for that sample — keeping only sites where that individual actually carries a non-reference call. That's standard practice for a distributable per-sample VCF (genome-wide `EMIT_ALL_SITES` output for one sample is huge and mostly `0/0` noise), and it directly explains a ~10x gap: `snp_calls_to_vcf` writes a record for every site in the joint cohort-wide sites list regardless of what this sample's genotype is there, while a variants-only release VCF only keeps the subset where this sample is actually non-reference — and for one sample against a list of sites variable across a large cohort, "most sites are ref for any given individual" easily gets you to roughly a 9-in-10 ratio.

# **Concrete way to confirm**: check what fraction of the "extra" records in your `snp_calls_to_vcf` output have `GT=0/0` (or `./.`) that don't appear at all in the reference VCF. If that fraction accounts for close to the full 10x gap, this is confirmed. If there's a meaningful number of *non*-reference extra records too, something else is going on (worth checking `site_mask` — `snp_calls_to_vcf` defaults to the resource's configured `site_mask`, so if you didn't pass one explicitly it's already filtering; if your reference VCF used a different or no filter, that's a second, smaller contributor, but filtering would only ever make the malariagen output *smaller*, not bigger, so it can't be the main driver of a 10x excess).

# If you want an apples-to-apples comparison rather than chasing this down, the cleaner options are: pull the actual `EMIT_ALL_SITES` raw genotyping output as ground truth (should match closely), or filter your `snp_calls_to_vcf` output down to non-`0/0`/non-missing sites for that sample before diffing against the release VCF.

# total number of SNPs
bcftools view -v snps SAMN12920115_2R.vcf.gz | grep -v "^#" | wc -l
#5102105
 
# total number of unique positions, indicating that several sites have two or more alternate alleles
bcftools view -v snps SAMN12920115_2R.vcf.gz | grep -v "^#" | cut -f2 | sort -u | wc -l
#5102105 

# distribution of ref vs. alt alleles
# notice the single dinucleotide change TG -> CG and an unnormalised variant AT -> AC
bcftools view -v snps SAMN12920115_2R.vcf.gz | grep -v "^#" | cut -f4,5 | sort | uniq -c | sort -k1rn
#1467256 T	A,C,G
# 1459989 A	C,T,G
# 1087870 G	A,C,T
# 1086990 C	A,T,G

bcftools view -v snps test_vcf_chrom_2R_claude_changes.vcf.gz | grep -v "^#" | wc -l

 
# total number of unique positions, indicating that several sites have two or more alternate alleles
bcftools view -v snps test_vcf_chrom_2R_claude_changes.vcf.gz | grep -v "^#" | cut -f2 | sort -u | wc -l
# distribution of ref vs. alt alleles
# notice the single dinucleotide change TG -> CG and an unnormalised variant AT -> AC
bcftools view -v snps test_vcf_chrom_2R_claude_changes.vcf.gz | grep -v "^#" | cut -f4,5 | sort | uniq -c | sort -k1rn

bcftools view -v snps test_vcf_chrom_2R_claude_changes_numpy.vcf.gz | grep -v "^#" | wc -l

# total number of unique positions, indicating that several sites have two or more alternate alleles
bcftools view -v snps test_vcf_chrom_2R_claude_changes_numpy.vcf.gz | grep -v "^#" | cut -f2 | sort -u | wc -l 

# distribution of ref vs. alt alleles
# notice the single dinucleotide change TG -> CG and an unnormalised variant AT -> AC
bcftools view -v snps test_vcf_chrom_2R_claude_changes_numpy.vcf.gz | grep -v "^#" | cut -f4,5 | sort | uniq -c | sort -k1rn

bcftools view -v snps test_vcf_chrom_2R.vcf.gz | grep -v "^#" | wc -l

 
# total number of unique positions, indicating that several sites have two or more alternate alleles
bcftools view -v snps test_vcf_chrom_2R.vcf.gz | grep -v "^#" | cut -f2 | sort -u | wc -l
#5102105 

# distribution of ref vs. alt alleles
# notice the single dinucleotide change TG -> CG and an unnormalised variant AT -> AC
bcftools view -v snps test_vcf_chrom_2R.vcf.gz | grep -v "^#" | cut -f4,5 | sort | uniq -c | sort -k1rn


bgzip test_vcf_chrom_2R_claude_changes_numpy.vcf -o test_vcf_chrom_2R_claude_changes_numpy.vcf.gz
bcftools index test_vcf_chrom_2R_claude_changes_numpy.vcf.gz

bcftools view -v snps test_vcf_chrom_2R_claude_changes_numpy.vcf.gz | grep -v "^#" | wc -l
#40561667
 
# total number of unique positions, indicating that several sites have two or more alternate alleles
bcftools view -v snps test_vcf_chrom_2R_claude_changes_numpy.vcf.gz | grep -v "^#" | cut -f2 | sort -u | wc -l
#40561667 

# distribution of ref vs. alt alleles
# notice the single dinucleotide change TG -> CG and an unnormalised variant AT -> AC
bcftools view -v snps test_vcf_chrom_2R_claude_changes_numpy.vcf.gz | grep -v "^#" | cut -f4,5 | sort | uniq -c | sort -k1rn
#10797167 T	A,C,G
# 10788549 A	C,T,G
# 9496455 C	A,T,G
# 9479496 G	A,C,T

# reference VCF
# snps: 5102105

#1467256 T	A,C,G
# 1459989 A	C,T,G
# 1087870 G	A,C,T
# 1086990 C	A,T,G

bcftools view -i 'GT[*]="alt"' test_vcf_chrom_2R_claude_changes_numpy.vcf.gz -o test_vcf_chrom_2R_claude_changes_numpy_homozygous_removed.vcf.gz
bcftools view -e 'F_PASS(GT="ref") == 1' test_vcf_chrom_2R_claude_changes_numpy.vcf -o test_vcf_chrom_2R_claude_changes_numpy_homozygous_removed2.vcf.gz

bcftools view -v snps test_vcf_chrom_2R_claude_changes_numpy_homozygous_removed.vcf.gz | grep -v "^#" | wc -l
#659196
 
# total number of unique positions, indicating that several sites have two or more alternate alleles
bcftools view -v snps test_vcf_chrom_2R_claude_changes_numpy_homozygous_removed.vcf.gz | grep -v "^#" | cut -f2 | sort -u | wc -l
#659196 

# distribution of ref vs. alt alleles
# notice the single dinucleotide change TG -> CG and an unnormalised variant AT -> AC
bcftools view -v snps test_vcf_chrom_2R_claude_changes_numpy_homozygous_removed.vcf.gz | grep -v "^#" | cut -f4,5 | sort | uniq -c | sort -k1rn
# 170268 C	A,T,G
# 169959 G	A,C,T
# 159542 A	C,T,G
# 159427 T	A,C,G


bcftools view -v snps test_vcf_chrom_2R_claude_changes_numpy_homozygous_removed2.vcf.gz | grep -v "^#" | wc -l
#874695
 
# total number of unique positions, indicating that several sites have two or more alternate alleles
bcftools view -v snps test_vcf_chrom_2R_claude_changes_numpy_homozygous_removed2.vcf.gz | grep -v "^#" | cut -f2 | sort -u | wc -l
#874695 

# distribution of ref vs. alt alleles
# notice the single dinucleotide change TG -> CG and an unnormalised variant AT -> AC
bcftools view -v snps test_vcf_chrom_2R_claude_changes_numpy_homozygous_removed2.vcf.gz | grep -v "^#" | cut -f4,5 | sort | uniq -c | sort -k1rn
# 223251 C	A,T,G
# 222817 G	A,C,T
# 214523 A	C,T,G
# 214104 T	A,C,G


bgzip test_vcf_chrom_2R_site_mask_none.vcf -o test_vcf_chrom_2R_site_mask_none.vcf.gz
bcftools index test_vcf_chrom_2R_site_mask_none.vcf.gz

bcftools view -v snps test_vcf_chrom_2R_site_mask_none.vcf.gz | grep -v "^#" | wc -l
#60132453
 
# total number of unique positions, indicating that several sites have two or more alternate alleles
bcftools view -v snps test_vcf_chrom_2R_site_mask_none.vcf.gz | grep -v "^#" | cut -f2 | sort -u | wc -l
#60132453 

# distribution of ref vs. alt alleles
# notice the single dinucleotide change TG -> CG and an unnormalised variant AT -> AC
bcftools view -v snps test_vcf_chrom_2R_site_mask_none.vcf.gz | grep -v "^#" | cut -f4,5 | sort | uniq -c | sort -k1rn
# 16457085 T	A,C,G
# 16437685 A	C,T,G
# 13623535 C	A,T,G
# 13614148 G	A,C,T

bcftools view -i 'GT[*]="alt"' test_vcf_chrom_2R_site_mask_none.vcf.gz -o test_vcf_chrom_2R_site_mask_none_homozygous_removed.vcf.gz
bcftools view -e 'F_PASS(GT="ref") == 1' test_vcf_chrom_2R_site_mask_none.vcf -o test_vcf_chrom_2R_site_mask_none_homozygous_removed2.vcf.gz

bcftools view -v snps test_vcf_chrom_2R_site_mask_none_homozygous_removed.vcf.gz | grep -v "^#" | wc -l
#1561739
 
bcftools view -v snps test_vcf_chrom_2R_site_mask_none_homozygous_removed.vcf.gz | grep -v "^#" | cut -f2 | sort -u | wc -l
#1561739 
bcftools view -v snps test_vcf_chrom_2R_site_mask_none_homozygous_removed.vcf.gz | grep -v "^#" | cut -f4,5 | sort | uniq -c | sort -k1rn

# 401263 T	A,C,G
# 400163 A	C,T,G
# 380247 G	A,C,T
# 380066 C	A,T,G

bcftools view -v snps test_vcf_chrom_2R_site_mask_none_homozygous_removed2.vcf.gz | grep -v "^#" | wc -l
#60132453
 
bcftools view -v snps test_vcf_chrom_2R_site_mask_none_homozygous_removed2.vcf.gz | grep -v "^#" | cut -f2 | sort -u | wc -l
#6663844 

bcftools view -v snps test_vcf_chrom_2R_site_mask_none_homozygous_removed2.vcf.gz | grep -v "^#" | cut -f4,5 | sort | uniq -c | sort -k1rn
# 1868519 T	A,C,G
# 1860152 A	C,T,G
# 1468117 G	A,C,T
# 1467056 C	A,T,G

# reference VCF
# snps: 5102105

# 1467256 T	A,C,G
# 1459989 A	C,T,G
# 1087870 G	A,C,T
# 1086990 C	A,T,G
bcftools index test_vcf_chrom_2R_site_mask_none_homozygous_removed2.vcf.gz
 
# homozygous removed 2 is closest so try intersection
bcftools isec -p compare_vcfs.txt SAMN12920115_2R.vcf.gz test_vcf_chrom_2R_site_mask_none_homozygous_removed2.vcf.gz

# https://davetang.org/muse/2019/09/02/comparing-vcf-files/
bedtools intersect -u -a SAMN12920115_2R.vcf.gz -b test_vcf_chrom_2R_site_mask_none_homozygous_removed2.vcf.gz | wc -l
#5102105
 
# calculate Jaccard index
bedtools jaccard -a SAMN12920115_2R.vcf.gz -b test_vcf_chrom_2R_site_mask_none_homozygous_removed2.vcf.gz
# intersection	union	jaccard	n_intersections
# 5102105	6663844	0.76564	256482

SnpSift concordance -v SAMN12920115_2R.vcf test_vcf_chrom_2R_site_mask_none_homozygous_removed2.vcf > snp_concordance.txt

vcf-compare SAMN12920115_2R.vcf.gz test_vcf_chrom_2R_site_mask_none_homozygous_removed2.vcf.gz
# SN	Number of REF matches:	5102105
# SN	Number of ALT matches:	5102105
# SN	Number of REF mismatches:	0
# SN	Number of ALT mismatches:	0
# SN	Number of samples in GT comparison:	0

# so all of the SNPs in the output VCF, once the heterozygous sites are removed, are in the reference 
# but there are 1561739 extra

# This is good news on correctness (zero false negatives against the reference confirms the GT-decoding path is sound), and the "extras" most likely have the same root cause as the original 10x-more-total-records finding — just showing up again one level down.

# The reference VCF was genotyped with `UnifiedGenotyper --alleles <fixed_list> --output_mode EMIT_ALL_SITES` — a **fixed, point-in-time site/allele list**, baked in whenever that specific pipeline run happened. The Zarr `snp_calls()` data backing `snp_calls_to_vcf`, by contrast, reflects the **current, full, cohort-wide joint site list** for the whole Ag3 release — which is very likely larger than whatever list existed when `SAMN12920115.vcf.gz` was originally produced, since the cohort has presumably grown (more samples added across releases → more segregating sites discovered) since that specific run. If this sample genuinely carries a non-reference allele at some of those *additional* sites — sites that simply weren't in the original narrower target list at all — you'd see exactly this: extra ALT calls in the malariagen output that have no counterpart of any kind in the reference VCF, not because of an error, but because the reference VCF was never even asked about those positions.

# **The one diagnostic that actually distinguishes "benign, expected" from "something's wrong" is per-position, not just presence/absence of the position overall.** For the ~extra ALT positions, split them into two groups:

# 1. **Position doesn't appear in the reference VCF at all** (no record, ref or alt) — this is the benign case. It means that position wasn't in whatever `--alleles` list generated the reference VCF, consistent with the cohort-wide site list simply having grown since then. Not a bug.
# 2. **Position *does* appear in the reference VCF, but as `0/0` (or missing) rather than ALT** — this is a genuine discordance at a position where you have real ground truth to compare against, and is worth digging into further (could be a real difference in calling, e.g. this sample having been re-genotyped with a newer pipeline version for the Ag3 release vs. whatever produced the original per-sample VCF, or could point at something in this codebase).

# If it's overwhelmingly case 1, I'd call this resolved and not a `to_vcf.py` issue — it's a reflection of the two VCFs having genuinely different site scopes by design, not a correctness problem in the exporter. If there's a meaningful chunk of case 2, that's worth showing me the specific positions and their `GQ`/`AD` values on both sides — at that point we're likely looking at a real difference in the underlying genotyping data provenance (pipeline version, reference genome patch, or cohort reprocessing) rather than anything in `snp_calls_to_vcf` itself, and that's genuinely outside what I can diagnose from the code alone — but I can help you dig through the malariagen side of a specific discordant position if you find one.

# Here's a command sequence that partitions the "extra" ALT positions cleanly into the two cases. This uses plain position-set comparison (`bcftools query` + `comm`) rather than `bcftools isec`, deliberately — `isec` matches on POS+REF+ALT by default, which would conflate "position not genotyped in reference" with "position genotyped but represented with different alleles," and you specifically want the position-only question here.

# **1. Get the reference's own ALT/het positions** (the set we already confirmed is fully contained in your malariagen ALT set — need this to correctly identify only the *extra* ones):

# ```bash
bcftools view -e 'F_PASS(GT="ref") == 1' SAMN12920115_2R.vcf.gz -o SAMN12920115_2R_alt_only.vcf.gz
bcftools index -t SAMN12920115_2R_alt_only.vcf.gz
# ```

# **2. Extract plain position lists** from all three files you need (malariagen ALT set, reference ALT set, reference *all* positions regardless of genotype):

# ```bash
bcftools query -f '%POS\n' test_vcf_chrom_2R_site_mask_none_homozygous_removed2.vcf.gz | sort -n -u > malariagen_alt_pos.txt
bcftools query -f '%POS\n' SAMN12920115_2R_alt_only.vcf.gz | sort -n -u > reference_alt_pos.txt
bcftools query -f '%POS\n' SAMN12920115_2R.vcf.gz | sort -n -u > reference_all_pos.txt
# ```

# (`sort -n -u` — numeric, deduplicated; since everything here is chromosome 2R only, plain `%POS` is enough, no need to carry `%CHROM` along.)

# **3. Isolate the "extra" ALT positions** (in malariagen's set, not in the reference's own ALT set):

# ```bash
comm -23 malariagen_alt_pos.txt reference_alt_pos.txt > extra_alt_pos.txt
wc -l extra_alt_pos.txt # 6187449
# ```

# **4. Split those into case 1 (absent from reference entirely) vs. case 2 (present in reference, just not as ALT):**

# ```bash
comm -23 extra_alt_pos.txt reference_all_pos.txt > case1_absent_from_reference.txt
comm -12 extra_alt_pos.txt reference_all_pos.txt > case2_present_as_nonalt.txt

wc -l case1_absent_from_reference.txt # 6187449
wc -l case2_present_as_nonalt.txt # 0

# Sanity check: `case1` + `case2` counts should sum to exactly the `extra_alt_pos.txt` count — a clean partition. If they don't, something's off with the position lists themselves (e.g. duplicate/overlapping records) rather than the underlying question, worth checking first.

# **5. If `case2` is non-trivial**, pull out a handful of those exact positions from both VCFs side by side to eyeball:

# ```bash
head -20 case2_present_as_nonalt.txt | awk '{print "2R\t"$1}' > spot_check_regions.txt

bcftools view -R spot_check_regions.txt test_vcf_chrom_2R_site_mask_none.vcf.gz | grep -v "^##"
bcftools view -R spot_check_regions.txt SAMN12920115_2R.vcf.gz | grep -v "^##"
# ```

# (`bcftools view -R` needs the inputs tabix-indexed — `bcftools index -t <file>` on each first if you haven't already.) That'll show you `GT`/`GQ`/`AD` side by side at the same positions for both call sets, which is the concrete evidence to bring back if `case2` turns out to be more than a handful of positions.

# Yes — and cleanly. Since case1 and case2 are a strict partition of extra_alt_pos.txt (every extra position falls into exactly one of the two by construction), case2 = 0 means case1 accounts for 100% of the extras. Worth a one-line confirmation that the arithmetic actually closes, just to rule out a sort/comm mishap rather than trusting the zero blindly:


wc -l extra_alt_pos.txt case1_absent_from_reference.txt
# If those two numbers match, this is a clean, confirmed result: every single "extra" ALT-genotype site in the malariagen output is a position that never appears in the reference VCF at all — not one of them is a position where the reference has a 0/0/missing call that snp_calls_to_vcf is wrongly reporting as ALT. That specifically rules out the "genuine miscall" explanation from case 2, and leaves the site-list-scope explanation (the reference VCF's --alleles list was fixed at whatever point that pipeline run happened, while the Zarr's joint site list reflects the current, presumably-larger cohort) as the confirmed cause, not just the leading hypothesis.

# Put together with the earlier result (every reference ALT/het site is present in the malariagen output — zero false negatives) and this one (zero false positives among the extras, in the sense that matters here), that's a genuinely clean bill of health for snp_calls_to_vcf's genotype-calling correctness on this sample/region: the count differences you originally set out to explain are fully accounted for by site-list scope and site_mask filtering, not by any bug in the exporter. Combined with the two real bugs we found and fixed earlier (FORMAT field order, MQ type declaration), I'd call this validation exercise closed out unless you want to repeat it on another sample or region as a second data point.


bgzip test_vcf_chrom_2R_non_ref_only.vcf -o test_vcf_chrom_2R_non_ref_only.vcf.gz
bcftools index test_vcf_chrom_2R_non_ref_only.vcf.gz

bcftools view -v snps test_vcf_chrom_2R_non_ref_only.vcf.gz | grep -v "^#" | wc -l
#1561739
 
# total number of unique positions, indicating that several sites have two or more alternate alleles
bcftools view -v snps test_vcf_chrom_2R_non_ref_only.vcf.gz | grep -v "^#" | cut -f2 | sort -u | wc -l
#1561739

# distribution of ref vs. alt alleles
# notice the single dinucleotide change TG -> CG and an unnormalised variant AT -> AC
bcftools view -v snps test_vcf_chrom_2R_non_ref_only.vcf.gz | grep -v "^#" | cut -f4,5 | sort | uniq -c | sort -k1rn
# 401263 T	A,C,G
# 400163 A	C,T,G
# 380247 G	A,C,T
# 380066 C	A,T,G

bcftools view -i 'GT[*]="alt"' test_vcf_chrom_2R_non_ref_only.vcf.gz -o test_vcf_chrom_2R_non_ref_only_none_homozygous_removed.vcf.gz
bcftools view -e 'F_PASS(GT="ref") == 1' test_vcf_chrom_2R_non_ref_only.vcf -o test_vcf_chrom_2R_non_ref_only_none_homozygous_removed2.vcf.gz

bcftools view -v snps test_vcf_chrom_2R_non_ref_only_none_homozygous_removed2.vcf.gz | grep -v "^#" | wc -l
#1561739
 
bcftools view -v snps test_vcf_chrom_2R_non_ref_only_none_homozygous_removed2.vcf.gz | grep -v "^#" | cut -f2 | sort -u | wc -l
#1561739 
bcftools view -v snps test_vcf_chrom_2R_non_ref_only_none_homozygous_removed2.vcf.gz | grep -v "^#" | cut -f4,5 | sort | uniq -c | sort -k1rn
# 401263 T	A,C,G
# 400163 A	C,T,G
# 380247 G	A,C,T
# 380066 C	A,T,G

bcftools view -v snps test_vcf_chrom_2R_non_ref_only_none_homozygous_removed.vcf.gz | grep -v "^#" | wc -l
#1561739
 
bcftools view -v snps test_vcf_chrom_2R_non_ref_only_none_homozygous_removed.vcf.gz | grep -v "^#" | cut -f2 | sort -u | wc -l
#1561739 
bcftools view -v snps test_vcf_chrom_2R_non_ref_only_none_homozygous_removed.vcf.gz | grep -v "^#" | cut -f4,5 | sort | uniq -c | sort -k1rn
# 401263 T	A,C,G
# 400163 A	C,T,G
# 380247 G	A,C,T
# 380066 C	A,T,G

# do the same double checks against reference tomorrow
bgzip test_vcf_chrom_2R_non_ref_only.vcf -o test_vcf_chrom_2R_non_ref_only.vcf.gz
bcftools index test_vcf_chrom_2R_non_ref_only.vcf.gz

bcftools isec -p compare_vcfs_new.txt SAMN12920115_2R.vcf.gz test_vcf_chrom_2R_non_ref_only_none_homozygous_removed.vcf.gz
bedtools jaccard -a SAMN12920115_2R.vcf.gz -b test_vcf_chrom_2R_non_ref_only.vcf.gz

bedtools intersect -u -a SAMN12920115_2R.vcf.gz -b test_vcf_chrom_2R_non_ref_only_none_homozygous_removed.vcf.gz | wc -l
# reference VCF
# snps: 5102105

# 1467256 T	A,C,G
# 1459989 A	C,T,G
# 1087870 G	A,C,T
# 1086990 C	A,T,G


bcftools isec -p compare_vcfs_new.txt SAMN12920115_2R.vcf.gz test_vcf_chrom_2R_non_ref_only_none_homozygous_removed.vcf.gz
bedtools jaccard -a SAMN12920115_2R.vcf.gz -b test_vcf_chrom_2R_non_ref_only_none_homozygous_removed.vcf.gz

bedtools intersect -u -a SAMN12920115_2R.vcf.gz -b test_vcf_chrom_2R_non_ref_only.vcf.gz | wc -l
bedtools jaccard -a SAMN12920115_2R.vcf.gz -b test_vcf_chrom_2R_non_ref_only.vcf.gz

# 2R      182     .       A       C,T,G   .       .       MQ=38.76;AC=0,0,0;AN=0  GT      ./.

# All clean. Summary of this fix:

# **Root cause**: `to_vcf.py` drops empty `variant_allele` slots when building the `ALT` column (so a site with only one real ALT gets a single-entry `ALT` string), but the genotype indices (`a0`/`a1` from `call_genotype`) were written raw — still referring to the *original*, uncompacted slot positions. Whenever a site's real ALT allele sits in slot 2 or 3 while an earlier slot is empty, this produces a genotype like `0/2` against an `ALT` column with only one allele listed — a malformed record with no valid referent for index 2. bcftools doesn't error on it, but silently excludes such rows from `GT="alt"`/`F_PASS(GT="ref")` classification, which is exactly the mechanism behind the 60M-vs-6.66M discrepancy.

# **Fix** ([to_vcf.py](malariagen_data/anoph/to_vcf.py)): decode `variant_allele` once per chunk, compute a per-variant "compacted index" mapping (original slot → its position in the written `ALT` list, or `-1` if the slot is empty), and remap `a0`/`a1` through that table before building the `GT` string — vectorized across the chunk, preserving the earlier performance work rather than falling back to a per-row Python loop. The ALT-building code in the per-variant loop now reuses the same decoded-alleles array, so the two are guaranteed consistent by construction rather than by two separate decode passes potentially drifting apart. A remapped index of `-1` (a genotype referencing an allele slot that doesn't actually exist at that site — shouldn't happen for well-formed data) is treated as missing (`./.`) defensively, rather than ever writing a negative index.

# **Verification**: reproduced the bug standalone first (`0/2` against a single-allele ALT), confirmed the fix resolves it plus three other gap patterns (gap in slot 1/2/3, and a two-ALT case with a gap), then added `test_vcf_exporter_allele_index_remapping` covering all four scenarios with an explicit invariant check (every GT allele index must be ≤ the number of ALT alleles actually listed on that row). Verified the test genuinely catches the regression by reverting the fix and confirming it fails with the exact `0/2` mismatch, then restored it. 32/32 tests pass, ruff and mypy clean, and the wider VCF/CNV/coverage/IGV/SNP suite (232 tests) is unaffected.

# This was a real, independent correctness bug — not something caused by anything else we fixed today, and one that would affect any multiallelic site in any `snp_calls_to_vcf` export, not just the `non_ref_only` path. Worth regenerating your file once more.


# **2. Extract plain position lists** from all three files you need (malariagen ALT set, reference ALT set, reference *all* positions regardless of genotype):

# = 
bcftools query -f '%POS\n' test_vcf_chrom_2R_non_ref_only.vcf.gz | sort -n -u > malariagen_alt_pos.txt
bcftools query -f '%POS\n' SAMN12920115_2R_alt_only.vcf.gz | sort -n -u > reference_alt_pos.txt
bcftools query -f '%POS\n' SAMN12920115_2R.vcf.gz | sort -n -u > reference_all_pos.txt

# 
comm -23 malariagen_alt_pos.txt reference_alt_pos.txt > extra_alt_pos.txt
wc -l extra_alt_pos.txt # 1561739
# ```

# ```bash
comm -23 extra_alt_pos.txt reference_all_pos.txt > case1_absent_from_reference.txt
comm -12 extra_alt_pos.txt reference_all_pos.txt > case2_present_as_nonalt.txt

wc -l case1_absent_from_reference.txt # 1561739
wc -l case2_present_as_nonalt.txt # 0


# sanity check
head -20 case2_present_as_nonalt.txt | awk '{print "2R\t"$1}' > spot_check_regions.txt

bcftools view -R spot_check_regions.txt test_vcf_chrom_2R_site_mask_none.vcf.gz | grep -v "^##"
bcftools view -R spot_check_regions.txt SAMN12920115_2R.vcf.gz | grep -v "^##"

bgzip test_vcf_chrom_2R_non_ref_only_non_remapped_alleles.vcf -o test_vcf_chrom_2R_non_ref_only_non_remapped_alleles.vcf.gz
bcftools index test_vcf_chrom_2R_non_ref_only_non_remapped_alleles.vcf.gz

bcftools isec -p compare_vcfs_new.txt SAMN12920115_2R.vcf.gz test_vcf_chrom_2R_non_ref_only_non_remapped_alleles.vcf.gz
bedtools jaccard -a SAMN12920115_2R.vcf.gz -b test_vcf_chrom_2R_non_ref_only_non_remapped_alleles.vcf.gz

bedtools intersect -u -a SAMN12920115_2R.vcf.gz -b test_vcf_chrom_2R_non_ref_only_non_remapped_alleles.vcf.gz | wc -l


bedtools intersect -u -a SAMN12920115_2R.vcf.gz -b test_vcf_chrom_2R_non_ref_only_none_homozygous_removed.vcf.gz > /tmp/isect.vcf
echo "exit code: $?"
wc -l /tmp/isect.vcf

bcftools isec -p /tmp/isec_out -n=2 SAMN12920115_2R.vcf.gz test_vcf_chrom_2R_non_ref_only_none_homozygous_removed.vcf.gz
wc -l /tmp/isec_out/0002.vcf  # records present in both


bgzip test_vcf_chrom_2R_non_ref_only_claude_rewrite.vcf -o test_vcf_chrom_2R_non_ref_only_claude_rewrite.vcf.gz
bcftools index test_vcf_chrom_2R_non_ref_only_claude_rewrite.vcf.gz

bcftools isec -p compare_vcfs_new.txt SAMN12920115_2R.vcf.gz test_vcf_chrom_2R_non_ref_only_claude_rewrite.vcf.gz
bedtools jaccard -a SAMN12920115_2R.vcf.gz -b test_vcf_chrom_2R_non_ref_only_claude_rewrite.vcf.gz

bedtools intersect -u -a SAMN12920115_2R.vcf.gz -b test_vcf_chrom_2R_non_ref_only_claude_rewrite.vcf.gz | wc -l

bcftools query -f '%POS\n' test_vcf_chrom_2R_non_ref_only.vcf.gz | sort -n -u > malariagen_alt_pos.txt
# 
comm -23 malariagen_alt_pos.txt reference_alt_pos.txt > extra_alt_pos.txt

bedtools intersect -u -a SAMN12920115_2R.vcf.gz -b test_vcf_chrom_2R_non_ref_only_none_homozygous_removed.vcf.gz > /tmp/isect.vcf
echo "exit code: $?"
wc -l /tmp/isect.vcf

bcftools isec -p /tmp/isec_out -n=2 SAMN12920115_2R.vcf.gz test_vcf_chrom_2R_non_ref_only_none_homozygous_removed.vcf.gz
wc -l /tmp/isec_out/0002.vcf  # records present in both

VCF_FILENAME="test_vcf_chrom_2R_restarting_exclude_0_0_refactor_single_sample.vcf"

bgzip $VCF_FILENAME -o ${VCF_FILENAME}.gz
bcftools index ${VCF_FILENAME}.gz

bcftools isec -p compare_vcfs_new.txt SAMN12920115_2R.vcf.gz ${VCF_FILENAME}.gz
bedtools jaccard -a SAMN12920115_2R.vcf.gz -b ${VCF_FILENAME}.gz
# 5102105	60132453	0.0848478	256482

bedtools intersect -u -a SAMN12920115_2R.vcf.gz -b ${VCF_FILENAME}.gz | wc -l

bcftools view -e 'F_PASS(GT="ref") == 1' ${VCF_FILENAME} -o ${VCF_FILENAME%.vcf}_homozygous_removed.vcf.gz

bcftools index ${VCF_FILENAME%.vcf}_homozygous_removed.vcf.gz

#bcftools isec -p compare_vcfs_new.txt SAMN12920115_2R.vcf.gz ${VCF_FILENAME}_homozygous_removed.vcf.gz
bedtools jaccard -a SAMN12920115_2R.vcf.gz -b ${VCF_FILENAME%.vcf}_homozygous_removed.vcf.gz
# 5102105	6663844	0.76564	256482

bedtools intersect -u -a SAMN12920115_2R.vcf.gz -b ${VCF_FILENAME%.vcf}_homozygous_removed.vcf.gz | wc -l

bcftools query -f '%POS\n' ${VCF_FILENAME%.vcf}_homozygous_removed.vcf.gz | sort -n -u > ${VCF_FILENAME%.vcf}_alt_pos.txt
# 
comm -23 ${VCF_FILENAME%.vcf}_alt_pos.txt reference_alt_pos.txt > ${VCF_FILENAME%.vcf}_extra_alt_pos.txt
wc -l ${VCF_FILENAME%.vcf}_extra_alt_pos.txt




bedtools intersect -u -a SAMN12920115_2R.vcf.gz -b test_vcf_chrom_2R_non_ref_only_none_homozygous_removed.vcf.gz > /tmp/isect.vcf
echo "exit code: $?"
wc -l /tmp/isect.vcf

bcftools isec -p /tmp/isec_out -n=2 SAMN12920115_2R.vcf.gz test_vcf_chrom_2R_non_ref_only_none_homozygous_removed.vcf.gz
wc -l /tmp/isec_out/0002.vcf  