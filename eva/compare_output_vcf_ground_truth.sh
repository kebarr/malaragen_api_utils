# need to compare vcf output to vcf downloaed based on accession
bgzip SAMN12920115.vcf -o SAMN12920115.vcf.gz # 60973065 SNPs in total

bcftools index SAMN12920115.vcf.gz
bcftools view -r "2R" SAMN12920115.vcf.gz -u -o SAMN12920115_2R.vcf
bcftools view -r "2L" SAMN12920115.vcf.gz -u -o SAMN12920115_2L.vcf

bcftools view -r "2R" SAMN12920115.vcf.gz -o SAMN12920115_2R_no_u_passed_to_bcftools_view.vcf


bgzip SAMN12920115_2R.vcf -o SAMN12920115_2R.vcf.gz # 5102105 SNPs in total

bcftools index SAMN12920115_2R.vcf.gz


bcftools isec -p compare_vcfs SAMN12920115.vcf.gz test_vcf.vcf.gz

# https://davetang.org/muse/2019/09/02/comparing-vcf-files/
bedtools intersect -u -a SAMN12920115.vcf.gz -b test_vcf.vcf.gz | wc -l

 
# calculate Jaccard index
bedtools jaccard -a SAMN12920115.vcf.gz -b test_vcf.vcf.gz


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


# I pulled up that WDL file to check exactly what the pipeline does, and it actually points fairly clearly at the explanation.

# The `UnifiedGenotyper` task runs with `--alleles ~{alleles_vcf}` and `--output_mode EMIT_ALL_SITES`. That means:

# - It genotypes against a **fixed, predetermined site list** (`alleles_vcf`) shared across the whole cohort — not sites discovered per-sample. That's the same kind of joint/cohort-wide site list that backs `snp_calls()`: every sample set's `call_genotype` array is aligned to one shared genome-wide sites list, not a per-sample one.
# - `EMIT_ALL_SITES` means the raw task output, if that's what you're comparing against, should in principle contain a record at *every* site in that shared list for the sample — reference calls (`0/0`) included, same as `snp_calls_to_vcf` does.

# So if your "reference VCF" is that raw, unfiltered per-task output, the site counts should be close to matching (modulo `site_mask`, more below). But if what you downloaded is the **released, distributed per-sample VCF** (e.g. via `wgs_data_catalog`'s `snp_genotypes_vcf` column) rather than the raw task output, it's very likely been through a later pipeline step that strips homozygous-reference records for that sample — keeping only sites where that individual actually carries a non-reference call. That's standard practice for a distributable per-sample VCF (genome-wide `EMIT_ALL_SITES` output for one sample is huge and mostly `0/0` noise), and it directly explains a ~10x gap: `snp_calls_to_vcf` writes a record for every site in the joint cohort-wide sites list regardless of what this sample's genotype is there, while a variants-only release VCF only keeps the subset where this sample is actually non-reference — and for one sample against a list of sites variable across a large cohort, "most sites are ref for any given individual" easily gets you to roughly a 9-in-10 ratio.

# **Concrete way to confirm**: check what fraction of the "extra" records in your `snp_calls_to_vcf` output have `GT=0/0` (or `./.`) that don't appear at all in the reference VCF. If that fraction accounts for close to the full 10x gap, this is confirmed. If there's a meaningful number of *non*-reference extra records too, something else is going on (worth checking `site_mask` — `snp_calls_to_vcf` defaults to the resource's configured `site_mask`, so if you didn't pass one explicitly it's already filtering; if your reference VCF used a different or no filter, that's a second, smaller contributor, but filtering would only ever make the malariagen output *smaller*, not bigger, so it can't be the main driver of a 10x excess).

# If you want an apples-to-apples comparison rather than chasing this down, the cleaner options are: pull the actual `EMIT_ALL_SITES` raw genotyping output as ground truth (should match closely), or filter your `snp_calls_to_vcf` output down to non-`0/0`/non-missing sites for that sample before diffing against the release VCF.



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

 
# reference VCF
# snps: 5102105

# 1467256 T	A,C,G
# 1459989 A	C,T,G
# 1087870 G	A,C,T
# 1086990 C	A,T,G


VCF_FILENAME="test_vcf_chrom_2R_restarting_exclude_0_0_refactor_single_sample.vcf"

bgzip $VCF_FILENAME -o ${VCF_FILENAME}.gz
bcftools index ${VCF_FILENAME}.gz


# total number of unique positions, indicating that several sites have two or more alternate alleles
bcftools view -v snps ${VCF_FILENAME}.gz | grep -v "^#" | cut -f2 | sort -u | wc -l

bcftools isec -p compare_vcfs_new.txt SAMN12920115_2R.vcf.gz ${VCF_FILENAME}.gz
bedtools jaccard -a SAMN12920115_2R.vcf.gz -b ${VCF_FILENAME}.gz
# 5102105	60132453	0.0848478	256482

bedtools intersect -u -a SAMN12920115_2R.vcf.gz -b ${VCF_FILENAME}.gz | wc -l

bcftools view -e 'F_PASS(GT="ref") == 1' ${VCF_FILENAME} -o ${VCF_FILENAME%.vcf}_homozygous_removed.vcf.gz

bcftools index ${VCF_FILENAME%.vcf}_homozygous_removed.vcf.gz

bedtools jaccard -a SAMN12920115_2R.vcf.gz -b ${VCF_FILENAME%.vcf}_homozygous_removed.vcf.gz
# 5102105	6663844	0.76564	256482

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


bedtools intersect -u -a SAMN12920115_2R.vcf.gz -b ${VCF_FILENAME%.vcf}_homozygous_removed.vcf.gz | wc -l

# **2. Extract plain position lists** from all three files you need (malariagen ALT set, reference ALT set, reference *all* positions regardless of genotype):


bcftools query -f '%POS\n' ${VCF_FILENAME%.vcf}_homozygous_removed.vcf.gz | sort -n -u > ${VCF_FILENAME%.vcf}_alt_pos.txt
bcftools query -f '%POS\n' SAMN12920115_2R_alt_only.vcf.gz | sort -n -u > reference_alt_pos.txt
bcftools query -f '%POS\n' SAMN12920115_2R.vcf.gz | sort -n -u > reference_all_pos.txt``
# 
# **3. Isolate the "extra" ALT positions** (in malariagen's set, not in the reference's own ALT set):

comm -23 ${VCF_FILENAME%.vcf}_alt_pos.txt reference_alt_pos.txt > extra_alt_pos.txt
wc -l extra_alt_pos.txt #

# **4. Split those into case 1 (absent from reference entirely) vs. case 2 (present in reference, just not as ALT):**


comm -23 extra_alt_pos.txt reference_all_pos.txt > case1_absent_from_reference.txt
comm -12 extra_alt_pos.txt reference_all_pos.txt > case2_present_as_nonalt.txt

wc -l case1_absent_from_reference.txt # 6187449
wc -l case2_present_as_nonalt.txt # 0

comm -23 ${VCF_FILENAME%.vcf}_alt_pos.txt reference_alt_pos.txt > ${VCF_FILENAME%.vcf}_extra_alt_pos.txt
wc -l ${VCF_FILENAME%.vcf}_extra_alt_pos.txt

# Sanity check: `case1` + `case2` counts should sum to exactly the `extra_alt_pos.txt` count — a clean partition. If they don't, something's off with the position lists themselves (e.g. duplicate/overlapping records) rather than the underlying question, worth checking first.

# **5. If `case2` is non-trivial**, pull out a handful of those exact positions from both VCFs side by side to eyeball:

#
head -20 case2_present_as_nonalt.txt | awk '{print "2R\t"$1}' > spot_check_regions.txt

bcftools view -R spot_check_regions.txt test_vcf_chrom_2R_site_mask_none.vcf.gz | grep -v "^##"
bcftools view -R spot_check_regions.txt SAMN12920115_2R.vcf.gz | grep -v "^##"
# ```

# (`bcftools view -R` needs the inputs tabix-indexed — `bcftools index -t <file>` on each first if you haven't already.) That'll show you `GT`/`GQ`/`AD` side by side at the same positions for both call sets, which is the concrete evidence to bring back if `case2` turns out to be more than a handful of positions.

# Yes — and cleanly. Since case1 and case2 are a strict partition of extra_alt_pos.txt (every extra position falls into exactly one of the two by construction), case2 = 0 means case1 accounts for 100% of the extras. Worth a one-line confirmation that the arithmetic actually closes, just to rule out a sort/comm mishap rather than trusting the zero blindly:


wc -l extra_alt_pos.txt case1_absent_from_reference.txt
# If those two numbers match, this is a clean, confirmed result: every single "extra" ALT-genotype site in the malariagen output is a position that never appears in the reference VCF at all — not one of them is a position where the reference has a 0/0/missing call that snp_calls_to_vcf is wrongly reporting as ALT. That specifically rules out the "genuine miscall" explanation from case 2, and leaves the site-list-scope explanation (the reference VCF's --alleles list was fixed at whatever point that pipeline run happened, while the Zarr's joint site list reflects the current, presumably-larger cohort) as the confirmed cause, not just the leading hypothesis.

# Put together with the earlier result (every reference ALT/het site is present in the malariagen output — zero false negatives) and this one (zero false positives among the extras, in the sense that matters here), that's a genuinely clean bill of health for snp_calls_to_vcf's genotype-calling correctness on this sample/region: the count differences you originally set out to explain are fully accounted for by site-list scope and site_mask filtering, not by any bug in the exporter. Combined with the two real bugs we found and fixed earlier (FORMAT field order, MQ type declaration), I'd call this validation exercise closed out unless you want to repeat it on another sample or region as a second data point.

bedtools intersect -u -a SAMN12920115_2R.vcf.gz -b test_vcf_chrom_2R_non_ref_only_none_homozygous_removed.vcf.gz > /tmp/isect.vcf
echo "exit code: $?"
wc -l /tmp/isect.vcf

bcftools isec -p /tmp/isec_out -n=2 SAMN12920115_2R.vcf.gz test_vcf_chrom_2R_non_ref_only_none_homozygous_removed.vcf.gz
wc -l /tmp/isec_out/0002.vcf  


## test bgzf compressed output

VCF_FILENAME="test_vcf_chrom_2R_no_pysam.vcf.gz"


# total number of unique positions, indicating that several sites have two or more alternate alleles
bcftools view -v snps ${VCF_FILENAME} | grep -v "^#" | cut -f2 | sort -u | wc -l # 6663844
bcftools view -v snps SAMN12920115_2R.vcf.gz | grep -v "^#" | cut -f2 | sort -u | wc -l # 5102105

bcftools index ${VCF_FILENAME}
bcftools isec -p compare_vcfs_no_pysam SAMN12920115_2R.vcf.gz ${VCF_FILENAME}


# 2. Sort the second VCF file

bcftools view -e 'F_PASS(GT="ref") == 1' ${VCF_FILENAME} -o ${VCF_FILENAME%.vcf.gz}_homozygous_removed.vcf.gz
bcftools view -e 'F_PASS(GT="ref") == 1' SAMN12920115_2R.vcf.gz -o SAMN12920115_2R_homozygous_removed.vcf.gz 

bcftools index ${VCF_FILENAME%.vcf.gz}_homozygous_removed.vcf.gz
bcftools index SAMN12920115_2R_homozygous_removed.vcf.gz

bcftools view -v snps SAMN12920115_2R_homozygous_removed.vcf.gz | grep -v "^#" | cut -f2 | sort -u | wc -l # 29302809
bcftools view -v snps ${VCF_FILENAME%.vcf.gz}_homozygous_removed.vcf.gz | grep -v "^#" | cut -f2 | sort -u | wc -l # 22573911


bcftools sort -O z -o ${VCF_FILENAME%.vcf.gz}_homozygous_removed.sorted.vcf.gz ${VCF_FILENAME%.vcf.gz}_homozygous_removed.vcf.gz
bcftools sort -O z -o SAMN12920115_2R_homozygous_removed.sorted.vcf.gz SAMN12920115_2R_homozygous_removed.vcf.gz

bcftools index ${VCF_FILENAME%.vcf.gz}_homozygous_removed.sorted.vcf.gz
bcftools index SAMN12920115_2R_homozygous_removed.sorted.vcf.gz

bedtools intersect -u -a SAMN12920115_2R_homozygous_removed.sorted.vcf.gz -b ${VCF_FILENAME%.vcf.gz}_homozygous_removed.sorted.vcf.gz | wc -l

bedtools jaccard -a SAMN12920115_2R_homozygous_removed.sorted.vcf.gz -b ${VCF_FILENAME%.vcf.gz}_homozygous_removed.sorted.vcf.gz # doesn't work even after sorting
# 5102105	6663844	0.76564	256482

bcftools isec -p compare_vcfs_full_ag3 SAMN12920115_2R_homozygous_removed.sorted.vcf.gz ${VCF_FILENAME%.vcf.gz}_homozygous_removed.sorted.vcf.gz

wc -l compare_vcfs_full_ag3/0002.vcf # Intesecting records: 27028131 compare_vcfs_full_ag3/0002.vcf
wc -l compare_vcfs_full_ag3/0000.vcf # Records in first file only: 11872282 compare_vcfs_full_ag3/0000.vcf
wc -l compare_vcfs_full_ag3/0001.vcf # Records in second file only: 18 compare_vcfs_full_ag3/0001.vcf

bedtools intersect -u -a SAMN12920115_2R_homozygous_removed.sorted.vcf.gz -b ${VCF_FILENAME%.vcf.gz}_homozygous_removed.sorted.vcf.gz | wc -l # 27028101
bedtools jaccard -a SAMN12920115_2R_homozygous_removed.sorted.vcf.gz -b ${VCF_FILENAME%.vcf.gz}_homozygous_removed.sorted.vcf.gz
# 3R	6	.	A	C,T,G	.	.	.	GT:GQ:AD:MQ	0/3:22:6,0,0,1:52
# still not working
bedtools jaccard \
  -a <(bcftools view -H SAMN12920115_2R_homozygous_removed.sorted.vcf.gz | sed 's/\r//g' | sort -k1,1 -k2,2n) \
  -b <(bcftools view -H ${VCF_FILENAME%.vcf.gz}_homozygous_removed.sorted.vcf.gz | sed 's/\r//g' | sort -k1,1 -k2,2n)


vcf-compare SAMN12920115_2R_homozygous_removed.sorted.vcf.gz ${VCF_FILENAME%.vcf.gz}_homozygous_removed.sorted.vcf.gz


# **2. Extract plain position lists** from all three files you need (malariagen ALT set, reference ALT set, reference *all* positions regardless of genotype):

bcftools view -e 'F_PASS(GT="ref") == 1' SAMN12920115_2R.vcf.gz -o SAMN12920115_2R_alt_only.vcf.gz

bcftools query -f '%POS\n' ${VCF_FILENAME%.vcf.gz}_homozygous_removed.sorted.vcf.gz | sort -n -u > ${VCF_FILENAME%.vcf.gz}_alt_pos_full.txt
bcftools query -f '%POS\n' SAMN12920115_2R_alt_only.vcf.gz | sort -n -u > reference_alt_pos_full.txt
bcftools query -f '%POS\n' SAMN12920115_2R.vcf.gz | sort -n -u > reference_all_pos_full.txt``
# 
# **3. Isolate the "extra" ALT positions** (in malariagen's set, not in the reference's own ALT set):

comm -23 ${VCF_FILENAME%.vcf.gz}_alt_pos_full.txt reference_alt_pos_full.txt > extra_alt_pos_full.txt
wc -l extra_alt_pos.txt # 1561739

# **4. Split those into case 1 (absent from reference entirely) vs. case 2 (present in reference, just not as ALT):**

comm -23 extra_alt_pos.txt reference_all_pos.txt > case1_absent_from_reference.txt
comm -12 extra_alt_pos.txt reference_all_pos.txt > case2_present_as_nonalt.txt

wc -l case1_absent_from_reference.txt # 1561739
wc -l case2_present_as_nonalt.txt # 0

# Sanity check: `case1` + `case2` counts should sum to exactly the `extra_alt_pos.txt` count — a clean partition. If they don't, something's off with the position lists themselves (e.g. duplicate/overlapping records) rather than the underlying question, worth checking first.

# **5. If `case2` is non-trivial**, pull out a handful of those exact positions from both VCFs side by side to eyeball:

#
head -20 case2_present_as_nonalt.txt | awk '{print "2R\t"$1}' > spot_check_regions.txt

bcftools view -R spot_check_regions.txt test_vcf_chrom_2R_site_mask_none.vcf.gz | grep -v "^##"
bcftools view -R spot_check_regions.txt SAMN12920115_2R.vcf.gz | grep -v "^##"
# ```

# (`bcftools view -R` needs the inputs tabix-indexed — `bcftools index -t <file>` on each first if you haven't already.) That'll show you `GT`/`GQ`/`AD` side by side at the same positions for both call sets, which is the concrete evidence to bring back if `case2` turns out to be more than a handful of positions.

# Yes — and cleanly. Since case1 and case2 are a strict partition of extra_alt_pos.txt (every extra position falls into exactly one of the two by construction), case2 = 0 means case1 accounts for 100% of the extras. Worth a one-line confirmation that the arithmetic actually closes, just to rule out a sort/comm mishap rather than trusting the zero blindly:


wc -l extra_alt_pos.txt case1_absent_from_reference.txt
# If those two numbers match, this is a clean, confirmed result: every single "extra" ALT-genotype site in the malariagen output is a position that never appears in the reference VCF at all — not one of them is a position where the reference has a 0/0/missing call that snp_calls_to_vcf is wrongly reporting as ALT. That specifically rules out the "genuine miscall" explanation from case 2, and leaves the site-list-scope explanation (the reference VCF's --alleles list was fixed at whatever point that pipeline run happened, while the Zarr's joint site list reflects the current, presumably-larger cohort) as the confirmed cause, not just the leading hypothesis.

# Put together with the earlier result (every reference ALT/het site is present in the malariagen output — zero false negatives) and this one (zero false positives among the extras, in the sense that matters here), that's a genuinely clean bill of health for snp_calls_to_vcf's genotype-calling correctness on this sample/region: the count differences you originally set out to explain are fully accounted for by site-list scope and site_mask filtering, not by any bug in the exporter. Combined with the two real bugs we found and fixed earlier (FORMAT field order, MQ type declaration), I'd call this validation exercise closed out unless you want to repeat it on another sample or region as a second data point.

# chroms in random example VCF:
##contig=<ID=2R>
##contig=<ID=2L>
##contig=<ID=3R>
##contig=<ID=3L>
##contig=<ID=X>

# chroms in SAMN12920115.vcf.gz:
##contig=<ID=2R,length=61545105>
##contig=<ID=3R,length=53200684>
##contig=<ID=2L,length=49364325>
##contig=<ID=UNKN,length=42389979>
##contig=<ID=3L,length=41963435>
##contig=<ID=X,length=24393108>
##contig=<ID=Y_unplaced,length=237045>
##contig=<ID=Mt,length=15363>
# so in ref VCF but not outpu, 42642387 sites
# sites in ref (before removal of homozygous sites) but not in output, 60973065 - 22573911= 38399154 sites
# sites in ref (after removal of homozygous sites) but not in output, 29302809 - 22573911= 6728898 sites

bcftools view -H SAMN12920115_homozygous_removed.sorted.vcf.gz | awk '{print $1}' | sort | uniq -c
# 6446528 2L
# 6663844 2R
# 4779460 3L
# 6257336 3R
#   30 Mt
# 11755964 UNKN
# 2880933 X
# 116258 Y_unplaced

# sites present in ref VCF but not output: 11872252

bcftools view -H ${VCF_FILENAME%.vcf.gz}_homozygous_removed.sorted.vcf.gz | awk '{print $1}' | sort | uniq -c

# 6446528 2L
# 6663844 2R
# 4779460 3L
# 6257336 3R
# 2880933 X

bcftools stats ${VCF_FILENAME%.vcf.gz}_homozygous_removed.sorted.vcf.gz > output_bcftools_stats2.txt
