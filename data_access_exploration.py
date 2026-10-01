import malariagen_data
import xarray

Ag3 = malariagen_data.Ag3(
    "simplecache::gs://vo_agam_release_master_us_central1",
    simplecache=dict(cache_storage="/Users/katie.barr/malariagen-data-python/gcs_cache_ag3"),
    results_cache="/Users/katie.barr/malariagen-data-python/results_cache",
)

filesystem = Ag3._fs
base_path = Ag3._base_path
config = Ag3.config

for key, value in Ag3.config.items():
    if "ZARR" in key.upper() or "GENOT" in key.upper():
        print(key, repr(value))

configured_zarr_path = 'reference/genome/agamp4/Anopheles-gambiae-PEST_SEQANNOTATION_AgamP4.12.zarr'

full_zarr_path = f"{Ag3._base_path}/{configured_zarr_path}"

store = Ag3._fs.get_mapper(full_zarr_path)

ds = xarray.open_zarr(
    store,
    consolidated=False,  # if .zmetadata is unavailable
)
## THIS WORKS!!!
import zarr
genotype_path = Ag3.config["SITE_ANNOTATIONS_ZARR_PATH"]

store = Ag3._fs.get_mapper(
    f"{Ag3._base_path}/{genotype_path}"
)

root = zarr.open_group(store=store, mode="r")
print(root.tree())