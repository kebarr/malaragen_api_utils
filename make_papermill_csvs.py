import csv
import json
from pathlib import Path

import malariagen_data

OUTPUT_DIR = Path(__file__).resolve().parent / "perf_testing_notebooks"
OUTPUT_PATH = OUTPUT_DIR / "papermill_parameter_file_filtered.csv"
FIELDNAMES = [
    "output_path",
    "region",
    "sample_sets",
    "sample_query",
    "sample_query_options",
    "site_mask",
    "max_region_size",
    "inline_array",
    "chunks",
    "overwrite",
    "fields",
]

af1 = malariagen_data.Af1()
sample_sets = [
    "1229-VO-GH-DADZIE-VMF00095",
    "1230-VO-GA-CF-AYALA-VMF00045",
    "1231-VO-MULTI-WONDJI-VMF00043",
    "1232-VO-KE-OCHOMO-VMF00044",
    "1240-VO-CD-KOEKEMOER-VMF00099",
    "1354-VO-KE-DONNELLY-VMF00281",
    "1330-VO-GN-LAMA-VMF00250",
    "1188-VO-NIANG-NIEL-SN-2304-VMF00259",
    "1281-VO-CM-CHRISTOPHE-VMF00227",
    "1236-VO-TZ-OKUMU-VMF00252",
    "1315-VO-NG-OMITOLA-OMOL-NG-2008",
    "1273-VO-ZM-MULEBA-VMF00176",
    "1354-VO-KE-DONNELLY-VMF00281",
]
end_indices = [2540000, 2640000, 2740000, 2840000, 2940000, 3040000, 3140000, 3240000, 3340000, 3440000]

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
with OUTPUT_PATH.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
    writer.writeheader()

    for i in range(10):
        for end_index in end_indices:
            for sample_set in sample_sets:
                print(f"Sample set: {sample_set}")
                sample_id = af1.sample_metadata(sample_sets=sample_set).iloc[i]["sample_id"]
                region = f"2RL:{end_index - 100000}-{end_index}"
                writer.writerow(
                    {
                        "output_path": f"local_igv_vcfs/{sample_id}_{region.replace(':', '_')}.vcf.gz",
                        "region": region,
                        "sample_sets": json.dumps([sample_set]),
                        "sample_query": f"sample_id == '{sample_id}'",
                        "inline_array": "true",
                        "overwrite": "true",
                        "fields": json.dumps(["GT", "GQ", "AD", "MQ"]),
                    }
                )


print(f"Wrote Papermill config CSV to {OUTPUT_PATH}")