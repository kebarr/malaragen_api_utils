#!/usr/bin/env bash

# Usage:
#   VCF_FILENAME="/path/to/input.vcf.gz" ./validate_vcf.sh
#
# Example:
#   VCF_FILENAME="/data/sample.vcf.gz" ./validate_vcf.sh

LOG_DIR="validation_logs"
RUN_TIMESTAMP=$(date '+%Y%m%d_%H%M%S')
RUN_NAME="${VCF_FILENAME##*/}"
RUN_NAME="${RUN_NAME%.vcf.gz}"
mkdir -p "$LOG_DIR"

CHECK_LOG="${LOG_DIR}/${RUN_NAME}_validation_${RUN_TIMESTAMP}.log"
PYTHON_RESOURCE_LOG="${LOG_DIR}/${RUN_NAME}_python_resources_${RUN_TIMESTAMP}.txt"
REFERENCE_PREFIX="SAMN12920115"

# Save all subsequent stdout and stderr from the checks to the run log,
# while also displaying it in the terminal.
exec > >(tee -a "$CHECK_LOG") 2>&1

: "${VCF_FILENAME:?Set VCF_FILENAME to the input VCF path}"
echo "Validation started: $(date)"
echo "Input VCF: $VCF_FILENAME"
echo "Expected generated reference VCF: ${REFERENCE_PREFIX}.vcf.gz"
echo "Python resource measurements: $PYTHON_RESOURCE_LOG"
#conda init
source activate base 
conda activate malariagen

echo
echo "========== Generate VCF from SNP calls =========="
if /usr/bin/time -l -o "$PYTHON_RESOURCE_LOG" \
    python3 make_vcf_from_snp_calls.py "$VCF_FILENAME"; then
    echo "Python script completed successfully."
else
    status=$?
    echo "Python script failed with exit status $status."
    echo "See $PYTHON_RESOURCE_LOG for elapsed time, CPU time, and peak memory."
    exit "$status"
fi
echo "Python resource measurements (macOS /usr/bin/time -l):"
cat "$PYTHON_RESOURCE_LOG"

echo
echo "========== Remaining VCF validation checks =========="

 
bcftools view -r '2R:2400000-2500000' \
  "${REFERENCE_PREFIX}.vcf.gz" \
  -Oz -o "${REFERENCE_PREFIX}_2R_2400000-2500000.vcf.gz"

bcftools index "${REFERENCE_PREFIX}_2R_2400000-2500000.vcf.gz"

# total number of unique positions, indicating that several sites have two or more alternate alleles
echo "Total unique positions in VCF: $(bcftools view -v snps ${VCF_FILENAME} | grep -v "^#" | cut -f2 | sort -u | wc -l)" >> ${LOG_DIR}/${RUN_NAME}_validation_${RUN_TIMESTAMP}.log
echo "Total unique positions in reference: $(bcftools view -v snps ${REFERENCE_PREFIX}_2R_2400000-2500000.vcf.gz | grep -v "^#" | cut -f2 | sort -u | wc -l)" >> ${LOG_DIR}/${RUN_NAME}_validation_${RUN_TIMESTAMP}.log

bcftools view -e 'F_PASS(GT="ref") == 1' ${VCF_FILENAME} -o ${VCF_FILENAME%.vcf.gz}_homozygous_removed.vcf.gz
bcftools view -e 'F_PASS(GT="ref") == 1' ${REFERENCE_PREFIX}_2R_2400000-2500000.vcf.gz -o ${REFERENCE_PREFIX}_2R_2400000-2500000_homozygous_removed.vcf.gz 

bcftools sort -O z -o ${VCF_FILENAME%.vcf.gz}_homozygous_removed.sorted.vcf.gz ${VCF_FILENAME%.vcf.gz}_homozygous_removed.vcf.gz
bcftools sort -O z -o ${REFERENCE_PREFIX}_2R_2400000-2500000_homozygous_removed.sorted.vcf.gz ${REFERENCE_PREFIX}_2R_2400000-2500000_homozygous_removed.vcf.gz

bcftools index ${VCF_FILENAME%.vcf.gz}_homozygous_removed.sorted.vcf.gz
bcftools index ${REFERENCE_PREFIX}_2R_2400000-2500000_homozygous_removed.sorted.vcf.gz

echo "Total unique positions in reference (homozygous removed): $(bcftools view -v snps ${REFERENCE_PREFIX}_2R_2400000-2500000_homozygous_removed.sorted.vcf.gz | grep -v "^#" | cut -f2 | sort -u | wc -l)" >> ${LOG_DIR}/${RUN_NAME}_validation_${RUN_TIMESTAMP}.log
echo "Total unique positions in VCF (homozygous removed): $(bcftools view -v snps ${VCF_FILENAME%.vcf.gz}_homozygous_removed.sorted.vcf.gz | grep -v "^#" | cut -f2 | sort -u | wc -l)" >> ${LOG_DIR}/${RUN_NAME}_validation_${RUN_TIMESTAMP}.log

mkdir ${LOG_DIR}/isec 
bcftools isec -p ${LOG_DIR}/isec ${REFERENCE_PREFIX}_2R_2400000-2500000_homozygous_removed.sorted.vcf.gz ${VCF_FILENAME%.vcf.gz}_homozygous_removed.sorted.vcf.gz

echo "Intersecting records: $(wc -l ${LOG_DIR}/isec/0002.vcf)" >> ${LOG_DIR}/${RUN_NAME}_validation_${RUN_TIMESTAMP}.log
echo "Records in first file only: $(wc -l ${LOG_DIR}/isec/0000.vcf)" >> ${LOG_DIR}/${RUN_NAME}_validation_${RUN_TIMESTAMP}.log
echo "Records in second file only: $(wc -l ${LOG_DIR}/isec/0001.vcf)" >> ${LOG_DIR}/${RUN_NAME}_validation_${RUN_TIMESTAMP}.log


echo "Bedtools intersect: $(bedtools intersect -u -a ${REFERENCE_PREFIX}_2R_2400000-2500000_homozygous_removed.sorted.vcf.gz -b ${VCF_FILENAME%.vcf.gz}_homozygous_removed.sorted.vcf.gz | wc -l)" >> ${LOG_DIR}/${RUN_NAME}_validation_${RUN_TIMESTAMP}.log

vcf-compare ${REFERENCE_PREFIX}_2R_2400000-2500000_homozygous_removed.sorted.vcf.gz ${VCF_FILENAME%.vcf.gz}_homozygous_removed.sorted.vcf.gz >> ${LOG_DIR}/${RUN_NAME}_validation_${RUN_TIMESTAMP}_vcf_compare.log

# **2. Extract plain position lists** from all three files you need (malariagen ALT set, reference ALT set, reference *all* positions regardless of genotype):

bcftools query -f '%POS\n' ${VCF_FILENAME%.vcf.gz}_homozygous_removed.sorted.vcf.gz | sort -n -u > ${VCF_FILENAME%.vcf.gz}_alt_pos_full.txt
bcftools query -f '%POS\n' ${REFERENCE_PREFIX}_2R_2400000-2500000_homozygous_removed.sorted.vcf.gz | sort -n -u > reference_alt_pos_full.txt
bcftools query -f '%POS\n' ${REFERENCE_PREFIX}_2R_2400000-2500000.vcf.gz | sort -n -u > reference_all_pos_full.txt
# **3. Isolate the "extra" ALT positions** (in malariagen's set, not in the reference's own ALT set):

comm -23 ${VCF_FILENAME%.vcf.gz}_alt_pos_full.txt reference_alt_pos_full.txt > extra_alt_pos_full.txt
echo "Extra ALT positions: $(wc -l extra_alt_pos_full.txt)" >> ${LOG_DIR}/${RUN_NAME}_validation_${RUN_TIMESTAMP}.log

# **4. Split those into case 1 (absent from reference entirely) vs. case 2 (present in reference, just not as ALT):**

comm -23 extra_alt_pos_full.txt reference_all_pos_full.txt > ${LOG_DIR}/case1_absent_from_reference.txt
comm -12 extra_alt_pos_full.txt reference_all_pos_full.txt > ${LOG_DIR}/case2_present_as_nonalt.txt

echo "Case 1 positions: $(wc -l ${LOG_DIR}/case1_absent_from_reference.txt)" >> ${LOG_DIR}/${RUN_NAME}_validation_${RUN_TIMESTAMP}.log
echo "Case 2 positions: $(wc -l ${LOG_DIR}/case2_present_as_nonalt.txt)" >> ${LOG_DIR}/${RUN_NAME}_validation_${RUN_TIMESTAMP}.log

# Sanity check: `case1` + `case2` counts should sum to exactly the `extra_alt_pos.txt` count — a clean partition. If they don't, something's off with the position lists themselves (e.g. duplicate/overlapping records) rather than the underlying question, worth checking first.

# **5. If `case2` is non-trivial**, pull out a handful of those exact positions from both VCFs side by side to eyeball:

echo "Spot check regions:" >> ${LOG_DIR}/${RUN_NAME}_validation_${RUN_TIMESTAMP}.log
echo "2R\t$(head -20 ${LOG_DIR}/case2_present_as_nonalt.txt | awk '{print $1}' | paste -sd " " -)" >> ${LOG_DIR}/${RUN_NAME}_validation_${RUN_TIMESTAMP}.log

bcftools view -R spot_check_regions.txt test_vcf_chrom_2R_site_mask_none.vcf.gz | grep -v "^##"
bcftools view -R spot_check_regions.txt ${REFERENCE_PREFIX}_2R.vcf.gz | grep -v "^##"

bcftools stats ${REFERENCE_PREFIX}_2R_2400000-2500000_homozygous_removed.sorted.vcf.gz > ${LOG_DIR}/${RUN_NAME}_ref_bcftools_stats.txt
bcftools stats ${VCF_FILENAME%.vcf.gz}_homozygous_removed.sorted.vcf.gz > ${LOG_DIR}/${RUN_NAME}_python_output_bcftools_stats.txt

echo "SNP counts:" >> ${LOG_DIR}/${RUN_NAME}_validation_${RUN_TIMESTAMP}.log
bcftools view -v snps ${VCF_FILENAME%.vcf.gz}_homozygous_removed.sorted.vcf.gz  | grep -v "^#" | cut -f4,5 | sort | uniq -c | sort -k1rn >> ${LOG_DIR}/${RUN_NAME}_validation_${RUN_TIMESTAMP}.log
bcftools view -v snps ${REFERENCE_PREFIX}_2R_2400000-2500000_homozygous_removed.sorted.vcf.gz   | grep -v "^#" | cut -f4,5 | sort | uniq -c | sort -k1rn >> ${LOG_DIR}/${RUN_NAME}_validation_${RUN_TIMESTAMP}.log