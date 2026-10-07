#!/usr/bin/env python3
# check output of to_vcf
import malariagen_data
import time

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


def load_ag3_write_vcf(output_path:str, simplecache_path:str,
                        results_cache_path:str)-> str:
    ag3 = malariagen_data.Ag3(
        "simplecache::gs://vo_agam_release_master_us_central1",
        simplecache=dict(cache_storage=simplecache_path),
        results_cache=results_cache_path,
    )

    res = check_wgs_data_catalog(ag3)
    sample_set_data = res[0]['bergey-2019']

    sample_id = sample_set_data["sample_id"].iloc[0]
    vcf_fileneme = sample_set_data["snp_genotypes_vcf"].iloc[0]
    # runs in a couple of minutes
    start = time.time()
    ag3.snp_calls_to_vcf(output_path=output_path, 
                         sample_sets = ['bergey-2019'],
                         sample_query = f"sample_id == '{sample_id}'",
                         region="2R:2,400,000-2,500,000",
                         fields = {"GT", "GQ", "AD", "MQ"},
                         overwrite = True)
    end = time.time()
    total_time = end-start
    print(f"Total time taken is: {total_time}")
    return vcf_fileneme


if __name__ == "__main__":
    import sys
    try:
        output_path = sys.argv[1]
    except IndexError:
        print("Please provide an output path for the VCF file.")
        sys.exit(1)
    if len(sys.argv) == 4:
        simplecache_path = sys.argv[2]
        results_cache_path = sys.argv[3]
    else:
        simplecache_path = "/Users/katie.barr/malariagen-data-python/gcs_cache_ag3"
        results_cache_path = "/Users/katie.barr/malariagen-data-python/results_cache"
    vcf_filename = load_ag3_write_vcf(output_path=output_path, simplecache_path=simplecache_path, results_cache_path=results_cache_path)
    print(f"{vcf_filename}")