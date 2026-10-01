

from typing import Any

import malariagen_data
from malariagen_data import adar1
from malariagen_data import adir1
ag3 = malariagen_data.Ag3(
    "simplecache::gs://vo_agam_release_master_us_central1",
     simplecache=dict(cache_storage="gcs_cache"), )
ag3.lookup_study_info


ag3_sample_sets = set(list(sample_metadata.sample_set))
study_infos = []
tos= []
for set in ag3_sample_sets:
    study_infos.append(ag3.lookup_study_info(set))
    tos.append(ag3.lookup_terms_of_use_info(set))



f = ag3.genome_features('X')

canonical = ag3.canonical_transcript("AGAP004707")
canonical = ag3.canonical_transcript("Pvr") # doesn'twork

q = ag3.cohort_geometries("admin1_quarter")

snp_effects = ag3.snp_effects("AGAP004707")
snp_calls = ag3.snp_calls('X')


# cnp_calls is a PITA, not sure if we need it, ditto snp_calls_to_vcf, plot_haplotype_sharing_arc

h = ag3.cohort_count_het("X")

# takes a million years don't run unless have to
h = ag3.cohort_heterozygosity("X", "admin1_quarter")


import malariagen_data
amin1 = malariagen_data.Amin1()


samples_metadata = amin1.sample_metadata()
amin1_sample_sets = set(list(samples_metadata.sample_set))

study_infos = []
tos= []
for s in amin1_sample_sets:
    study_infos.append(amin1.lookup_study_info(s))
    tos.append(amin1.lookup_terms_of_use_info(s))

gf = amin1.genome_features()
genes = gf.query("type == 'gene'")

valid_ids = genes["ID"]          # e.g. "AGAP004707"
valid_names = genes["Name"].dropna()  # e.g. "Pvr" (not every gene has a Name)

# then:
amin1.canonical_transcript(valid_ids.iloc[0])
q = amin1.cohort_geometries("admin1_quarter")


as1_name = "vo_aste_release_master_us_central1"


as1 = malariagen_data.As1(
    f"simplecache::gs://{as1_name}",
     simplecache=dict(cache_storage="gcs_cache"), )

sample_metadata = as1.sample_metadata()

as1_sample_sets = set(list(sample_metadata.sample_set))
study_infos = []
tos= []
for s in as1_sample_sets:
    study_infos.append(as1.lookup_study_info(s))
    tos.append(as1.lookup_terms_of_use_info(s))

Ag3 = malariagen_data.Ag3(
    "simplecache::gs://vo_agam_release_master_us_central1",
    simplecache=dict(cache_storage="gcs_cache_ag3"),
)
Af1 = malariagen_data.Af1(
        "simplecache::gs://vo_afun_release_master_us_central1",
        simplecache=dict(cache_storage="gcs_cache_af1"),
    )
As1=malariagen_data.As1(
        "simplecache::gs://vo_aste_release_master_us_central1",
        simplecache=dict(cache_storage="gcs_cache_as1"),
    )
Amin1= malariagen_data.Amin1(
        "simplecache::gs://vo_amin_release_master_us_central1",
        simplecache=dict(cache_storage="gcs_cache_amin1"),
    )
Adir1= malariagen_data.Adir1(
        "simplecache::gs://vo_adir_release_master_us_central1",
        simplecache=dict(cache_storage="gcs_cache_adir1"),
    )
Adar1= malariagen_data.Adar1(
        "simplecache::gs://vo_adar_release_master_us_central1",
        simplecache=dict(cache_storage="gcs_cache_adar1"),
    )

def can_get_canonical_transcript(species, name):
    gf = species.genome_features()
    genes = gf.query("type == 'gene'")

    valid_ids = genes["ID"]   
    print(f"Number of genes for {name}")
    print(len(valid_ids))      

    try:
        # then:
        species.canonical_transcript(valid_ids.iloc[0])
    except:
        print(f"No nanonical transcripts for {name}")


can_get_canonical_transcript(Ag3, "Ag3")
can_get_canonical_transcript(Af1, "Af1")
can_get_canonical_transcript(Amin1, "Amin1")
can_get_canonical_transcript(Adar1, "Adar1")
can_get_canonical_transcript(Adir1, "Adir1")
can_get_canonical_transcript(As1, "As1")

def can_get_study_info_tos(species, name):
    samples_metadata = species.sample_metadata()
    try:
        # then:
        sample_sets = set(list(samples_metadata.sample_set))

        study_infos = []
        tos= []
        for s in sample_sets:
            study_infos.append(species.lookup_study_info(s))
            tos.append(species.lookup_terms_of_use_info(s))
        print(f"Len lookup study info for {name}")
        print(len(study_infos))   
        print(f"Len lookup tos for {name}")
        print(len(tos))   
    except:
        print(f"No nanonical transcripts for {name}")


can_get_study_info_tos(Ag3, "Ag3")
can_get_study_info_tos(Af1, "Af1")
can_get_study_info_tos(Amin1, "Amin1")
can_get_study_info_tos(Adar1, "Adar1")
can_get_study_info_tos(Adir1, "Adir1")
can_get_study_info_tos(As1, "As1")


def wgs_run_accessions(species, name):
    samples_metadata = species.sample_metadata()
    try:
        # then:
        sample_sets = list(samples_metadata.sample_set)

        accs = species.wgs_run_accessions(sample_sets[0])
        print(f"wgs_run_accessions {name}")
        print(accs.shape)    
    except:
        print(f"No nanonical transcripts for {name}")


wgs_run_accessions(Ag3, "Ag3")
wgs_run_accessions(Af1, "Af1")
wgs_run_accessions(Amin1, "Amin1")
wgs_run_accessions(Adar1, "Adar1")
wgs_run_accessions(Adir1, "Adir1")
wgs_run_accessions(As1, "As1")


def get_random_sample_set(dataset):
    sample_sets =list(dataset.sample_sets()["sample_set"])
    return sample_sets[0]

def get_random_region(dataset):
    return dataset.contigs[0]


#### need to check each diagram is correct
Ag3 = malariagen_data.Ag3(
    "simplecache::gs://vo_agam_release_master_us_central1",
    simplecache=dict(cache_storage="gcs_cache_ag3"),
    results_cache="results_cache",
)
Af1 = malariagen_data.Af1(
        "simplecache::gs://vo_afun_release_master_us_central1",
        simplecache=dict(cache_storage="gcs_cache_af1"),
        results_cache="results_cache",
    )
As1=malariagen_data.As1(
        "simplecache::gs://vo_aste_release_master_us_central1",
        simplecache=dict(cache_storage="gcs_cache_as1"),
        results_cache="results_cache",
    )
Amin1= malariagen_data.Amin1(
        "simplecache::gs://vo_amin_release_master_us_central1",
        simplecache=dict(cache_storage="gcs_cache_amin1"),
        results_cache="results_cache",
    )
Adir1= malariagen_data.Adir1(
        "simplecache::gs://vo_adir_release_master_us_central1",
        simplecache=dict(cache_storage="gcs_cache_adir1"),
        results_cache="results_cache",
    )
Adar1= malariagen_data.Adar1(
        "simplecache::gs://vo_adar_release_master_us_central1",
        simplecache=dict(cache_storage="gcs_cache_adar1"),
        results_cache="results_cache",
    )

# Af1 should have all but AIM and phenotype data
Adar1.aim_ids # empty

#datasets = [Ag3, Af1, As1, Amin1,Adir1,Adar1]

datasets = [Af1, As1, Amin1,Adir1,Adar1]

def check_phenotype_data(dataset):
    phenotype_sample_sets = dataset.phenotype_sample_sets()
    print(f"Available phenotype sample sets: {phenotype_sample_sets}")

def test_xpehh_analysis(dataset):
    contig = dataset.contigs[0]
    try:
        x, xpehh = dataset.xpehh_gwss(contig)
        print(f"Results for dataaset {x}, xpehh:{xpehh}")
    except:
        print("No XPEHH")

# Needing appropriate inputs to validae:
#dipclust, hapclust, H1X - uses haplotypes so nothing without haplotypes will work, H12- same
#

def get_haplotype(dataset):
    region = get_random_region(dataset)
    haplotypes = dataset.haplotypes(region)
    if len(haplotypes) > 0:
        print(f"Found {len(haplotypes)} haplotypes for dataset")
    else:
        print("No haplotypes for dataset")

def get_cnv_calls(dataset, region, sample_set):
    cnv_calls = dataset.cnv_coverage_calls(region, sample_set)
    print(cnv_calls)


def get_snp_calls(dataset, region, sample_set):
    snp_calls = dataset.snp_calls(region=region, sample_sets=sample_set)
    print("SNP_CALLS SHAPE")
    print(snp_calls)

def get_bioallelic_snp_calls(dataset, region, sample_set):
    snp_calls = dataset.biallelic_snp_calls(region=region, sample_sets=sample_set)
    print("BIALELLIC SNP_CALLS SHAPE")
    print(snp_calls)


datasets = [Af1, As1, Amin1,Adir1,Adar1]

# If I can't get haplotypes at all then can't get anything that depends on haplotypes 

for dataset in datasets:
    print(f" checking for dataset:{dataset.config['GENOME_REF_NAME']}")
    # only implemented for Ag3, as expected
    #check_phenotype_data(dataset)
    #print(f"AIM data: {dataset.aim_ids}")
    #test_xpehh_analysis(dataset)
    sample_set = get_random_sample_set(dataset)
    region = get_random_region(dataset)
    #get_cnv_calls(dataset, region, sample_set)
    get_snp_calls(dataset, region, sample_set)
    get_bioallelic_snp_calls(dataset, region, sample_set)
    #get_haplotype(dataset)
    print("#############################################")
    print("Next dataset")    
    print("#############################################")
    print("\n\n\n\n")

Amin1_region = get_random_region(Amin1)
Amin1_sample = get_random_sample_set(Amin1)

amin1_haps = Amin1.haplotypes(region='KB663766', sample_sets=Amin1_sample, analysis='minimus')
# tried for a few regions, no luck




class Region:
    """A region of a reference genome, i.e., a contig or contig interval."""

    def __init__(self, contig, start=None, end=None):
        self._contig = contig
        self._start = start
        self._end = end

    @property
    def contig(self):
        return self._contig

    @property
    def start(self):
        return self._start

    @property
    def end(self):
        return self._end

    def __hash__(self):
        return hash((self.contig, self.start, self.end))

    def __eq__(self, other):
        return (
            isinstance(other, Region)
            and (self.contig == other.contig)
            and (self.start == other.start)
            and (self.end == other.end)
        )

    def __repr__(self):
        return f"Region({self._contig!r}, {self._start!r}, {self._end!r})"

    def __str__(self):
        out = self._contig
        if self._start is not None or self._end is not None:
            out += ":"
            if self._start is not None:
                out += f"{self._start:,}"
            out += "-"
            if self.end is not None:
                out += f"{self._end:,}"
        return out

    def to_dict(self):
        return dict(
            contig=self.contig,
            start=self.start,
            end=self.end,
        )
    

# can't do from ipython without some twiddling
def test_igv(dataset, region, sample_set):
    dataset.igv(region=region)
    sample_id = dataset.sample_metadata(sample_sets=sample_set).iloc[0]["sample_id"]
    ag3.view_alignments(sample=sample_id, region=region)


def test_pca(dataset, region, sample_set):
    df_pca, evr = dataset.pca(
        region=region,
        sample_sets=sample_set,
        n_snps=10_000,
    )
    print("PCA DF SUMMMARY")
    print(df_pca.describe())


def test_heterozygosity(dataset, sample_set):
    # df_cohort_samples = dataset.sample_metadata(
    #     sample_sets=sample_set,
    #     sample_query="cohort_admin2_month",
    # )
    # region_normalised = Region(region, None, None)
    # result = dataset.cohort_count_het(region_normalised, df_cohort_samples)
    # print("WINDOWS:")
    # list(result.keys())[:10]
    # print("counts")
    # list(result.values())[:10]
    all_sample_sets = dataset.sample_sets()["sample_set"].to_list()
    sample_set = str(np.random.choice(all_sample_sets))

    cohort_params = dict(
        region=str(np.random.choice(dataset.contigs)),
        cohorts="taxon",
        sample_sets=sample_set,
        window_size=20_000,
    )

    # Run function under test.
    df = dataset.cohort_heterozygosity(**cohort_params)

    # Check results.
    print(df)
    


# from test code
import numpy as np
def test_fst(dataset):
    all_sample_sets = dataset.sample_sets()["sample_set"].to_list()
    all_countries = dataset.sample_metadata()["country"].dropna().unique().tolist()
    try:
        countries = np.random.choice(all_countries, size=2, replace=False).tolist()
        cohort1_query = f"country == {countries[0]!r}"
        cohort2_query = f"country == {countries[1]!r}"
        fst_params = dict(
            region=np.random.choice(dataset.contigs, size=2, replace=False).tolist(),
            sample_sets=all_sample_sets,
            cohort1_query=cohort1_query,
            cohort2_query=cohort2_query,
            site_mask=str(np.random.choice(dataset.site_mask_ids)),
            min_cohort_size=1,
            n_jack=int(np.random.randint(10, 201)),
        )

        # Run function under test.
        fst, se = dataset.average_fst(**fst_params)
        print(f"fst: {fst}")
        print(f"se: {se}")
    except ValueError:
        print("Cannot compute without two cohorts")

def test_g123_snp_calls_only(dataset):
    try:
        all_sample_sets = dataset.sample_sets()["sample_set"].to_list()
        g123_params = dict(
            contig=str(np.random.choice(dataset.contigs)),
            sample_sets=[str(np.random.choice(all_sample_sets))],
            window_size=int(np.random.randint(100, 501)),
            min_cohort_size=10,
        )

        x, g123 = dataset.g123_gwss(**g123_params)
        print(f"X: {x}")
        print(f"g123: {g123}")
    except Exception as e:
        print(e)


for dataset in datasets:
    print(f" checking for dataset:{dataset.config['GENOME_REF_NAME']}")
    # only implemented for Ag3, as expected
    #check_phenotype_data(dataset)
    #print(f"AIM data: {dataset.aim_ids}")
    #test_xpehh_analysis(dataset)
    sample_set = get_random_sample_set(dataset)
    region = get_random_region(dataset)
    test_fst(dataset)
    #test_pca(dataset, region, sample_set)
    #test_heterozygosity(dataset, sample_set)
    test_g123_snp_calls_only(dataset)
    #test_igv(dataset, region, sample_set)
    #get_haplotype(dataset)
    print("#############################################")
    print("Next dataset")    
    print("#############################################")
    print("\n\n\n\n")


#Check the implementation and dataset configuration with:

# rg -n "def karyotype|inversion|karyotype|KARYOTYPE"  malariagen_data tests

#Then inspect the configured inversion IDs:
print(Adar1.karyotype(inversion="...", sample_sets=...))

#Use only inversion names explicitly supported by that dataset’s configuration. The notebook’s Ag3 examples demonstrate Ag3 support; they do not establish support for every dataset.

adar1_sequence = Adar1.open_genome()
adir1_sequence = Adir1.open_genome()
amin1_sequence = Amin1.open_genome()
as1_sequence = As1.open_genome()
af1_sequence = Af1.open_genome()

features = []
for dataset in datasets:
    print(f" checking for dataset:{dataset.config['GENOME_REF_NAME']}")
    region = get_random_region(dataset)
    print(region)
    feature = dataset.genome_features(region)
    features.append(feature)

region = get_random_region(Af1)
hap = Af1.haplotypes(region)


def list_available_cohorts(dataset):
    metadata = dataset.sample_metadata()
    cohort_columns = [
    column for column in metadata.columns
    if "cohort" in column.lower()
    ]
    return cohort_columns

def check_fst_support(
    dataset: Any,
    *,
    region: str,
    cohort1_query: str,
    cohort2_query: str,
    sample_sets: str | list[str] | None = None,
    site_mask: str | None = None,
) -> dict[str, Any]:
    """Check whether Fst can be calculated for a dataset and cohort pair."""
    if not callable(getattr(dataset, "average_fst", None)):
        return {
            "dataset": type(dataset).__name__,
            "supported": False,
            "reason": "average_fst() is not available",
        }

    try:
        fst, se = dataset.average_fst(
            region=region,
            cohort1_query=cohort1_query,
            cohort2_query=cohort2_query,
            sample_sets=sample_sets,
            site_mask=site_mask,
            n_jack=1,
        )
    except Exception as exc:
        return {
            "dataset": type(dataset).__name__,
            "supported": False,
            "reason": f"{type(exc).__name__}: {exc}",
        }

    return {
        "dataset": type(dataset).__name__,
        "supported": True,
        "fst": float(fst),
        "standard_error": float(se),
    }


result = check_fst_support(
    Adar1,
    region="2:1,000,000-2,000,000",
    cohort1_query="country == 'Colombia'",
    cohort2_query="country == 'Guyana'",
)

print(result)

# get countries:
set(Adar1.sample_metadata()["country"].to_list())
#and contigs
Adar1.contigs

set(Adir1.sample_metadata()["country"].to_list())


result = check_fst_support(
    Adir1,
    region="KB673646:1,000-2,000",
    cohort1_query="country == 'Bangladesh'",
    cohort2_query="country == 'Cambodia'",
)

result = check_fst_support(
    As1,
    region="2RL:1,000,000-2,000,000",
    cohort1_query="country == 'Iran'",
    cohort2_query="country == 'India'",
)

result = check_fst_support(
    As1,
    region="2RL:1,000,000-2,000,000",
    cohort1_query="country == 'Iran'",
    cohort2_query="country == 'India'",
)

result = check_fst_support(
    Af1,
    region="2RL:1,000,000-2,000,000",
    cohort1_query="country == 'Benin'",
    cohort2_query="country == 'Gambia'",
)

### for igv, can't get it running in VSCode but can from terminal with "jupyter lab"