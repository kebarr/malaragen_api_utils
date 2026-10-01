import malariagen_data


#### need to check each diagram is correct
Ag3 = malariagen_data.Ag3(
    "simplecache::gs://vo_agam_release_master_us_central1",
    simplecache=dict(cache_storage="../gcs_cache_ag3"),
    results_cache="../results_cache",
)
Af1 = malariagen_data.Af1(
        "simplecache::gs://vo_afun_release_master_us_central1",
        simplecache=dict(cache_storage="../gcs_cache_af1"),
        results_cache="../results_cache",
    )
As1=malariagen_data.As1(
        "simplecache::gs://vo_aste_release_master_us_central1",
        simplecache=dict(cache_storage="../gcs_cache_as1"),
        results_cache="../results_cache",
    )
Amin1= malariagen_data.Amin1(
        "simplecache::gs://vo_amin_release_master_us_central1",
        simplecache=dict(cache_storage="../gcs_cache_amin1"),
        results_cache="../results_cache",
    )
Adir1= malariagen_data.Adir1(
        "simplecache::gs://vo_adir_release_master_us_central1",
        simplecache=dict(cache_storage="../gcs_cache_adir1"),
        results_cache="../results_cache",
    )
Adar1= malariagen_data.Adar1(
        "simplecache::gs://vo_adar_release_master_us_central1",
        simplecache=dict(cache_storage="../gcs_cache_adar1"),
        results_cache="../results_cache",
    )

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

datasets = [Ag3, Af1, As1, Amin1, Adir1, Adar1]
dataset_names = ["Ag3", "Af1", "As1", "Amin1", "Adir1", "Adar1"]
dataset_accessions = {}
for dataset, dataset_name in zip(datasets, dataset_names):
    res = check_wgs_data_catalog(dataset)
    dataset_accessions[dataset_name] = res

def check_run_accessions(dataset)-> tuple:
    #sample_sets = dataset.sample_metadata()["sample_id"].to_list()
    sample_sets = dataset.sample_sets()["sample_set"].to_list() 
    accessions = {}
    samples_no_data = []
    for sample_set in sample_sets:
        try:
            result = dataset.wgs_run_accessions(sample_set)
            result.dropna(inplace=True)
            if result.shape[0] != 0:
                accessions[sample_set] = result
        except FileNotFoundError:
            samples_no_data.append(sample_set)
    return accessions, samples_no_data


run_accessions = {}
for dataset, dataset_name in zip(datasets, dataset_names):
    res = check_run_accessions(dataset)
    run_accessions[dataset_name] = res