import pandas as pd
import numpy as np
mut_path = "prj00215-variant-calls-annotated.tsv"

#note this script requires the mut file to be annotated with the metadata file

mut = pd.read_csv(mut_path, sep="\t", low_memory=False)
#loading... annotated mut file

print(mut.head())

mut["indv_read"] = mut["end"] - mut["start"]
#used to calculate total bases read per sample (calculates length of interval)

mut["is_mutation"] = mut["variation_type"].ne("no_variant")
#used to determine if mutation is present
regions = ["intergenic","intronic","exon_non_coding","coding"]
mut["region_categorise"] = np.select(
    [
    mut["region"].eq ("intergenic"),
    mut["region"].eq("intronic"),
    mut["region"].eq("exon_non_coding"),
    mut["region"].eq("coding"),
    ],
    regions, default = "other"
) #assigns regions to categories and if not recognised assigns other

mut["var_class"] = np.select(
    [mut["variation_type"].eq("no_variant"),
     mut["variation_type"].eq("snv"),
     mut["variation_type"].eq("mnv"),
     mut["variation_type"].eq("deletion"),
     mut["variation_type"].eq("insertion"),
     mut["variation_type"].eq("complex"),
     mut["variation_type"].eq("symbolic"),
     ],
    ["no_variant", "snv", "mnv", "deletion", "insertion", "complex", "symbolic"], default = "other"
) #assigns variation types to categories and if not recognised assigns other

base_summary = (mut.groupby(["sample"], as_index=False)
           .agg(
            total_read = ("indv_read", "sum"),
            mutation_event = ("is_mutation", "sum"),
            total_rows = ("sample", "size"),
)
) #creates summary by sample, calculates total bases read and mutations (whole sample summary)

region_summary = (mut.groupby(["sample", "region_categorise"], as_index=False)
                  .agg(region_assessed_bases = ("indv_read", "sum"),
                      region_mutation_event = ("is_mutation", "sum"),
                       )) #summarises assessed bases by region and sanple (i.e. interval lengths mutation events)
region_bases = (
    region_summary.pivot_table(index="sample", columns="region_categorise", values="region_assessed_bases", aggfunc="sum", fill_value=0)
    .reindex(columns= regions, fill_value=0).add_suffix("_assessed_bases").reset_index()
)  #summarises assessed bases by region and sample (used in later model offsets)

region_counts = (region_summary.pivot_table(index="sample", columns="region_categorise",
                                                     values="region_mutation_event",
                                                     aggfunc="sum",
                                                     fill_value=0)
                 .reindex(columns=regions, fill_value=0).add_suffix("_mutation_event").reset_index()
                 ) #counts mutations by region

base_summary = base_summary.merge(region_counts, how="left", on="sample")
base_summary = base_summary.merge(region_bases, how="left", on="sample")

base_summary["total_mutations_mf"] = (base_summary["mutation_event"]
                                      / base_summary["total_read"].replace(0, np.nan))
base_summary["intergenic_mutations_mf"] = (base_summary["intergenic_mutation_event"]
                                           / base_summary["intergenic_assessed_bases"].replace(0, np.nan))
base_summary["intronic_mutations_mf"] = (base_summary["intronic_mutation_event"]
                                         / base_summary["intronic_assessed_bases"].replace(0, np.nan))
base_summary["exon_non_coding_mf"] = (base_summary["exon_non_coding_mutation_event"]
                                      / base_summary["exon_non_coding_assessed_bases"].replace(0, np.nan))
base_summary["coding_mutations_mf"] = (base_summary["coding_mutation_event"] /
                                       base_summary["coding_assessed_bases"].replace(0, np.nan))
#calculates mutation frequencies (mutation event per assessed base) per sample and by region

meta_data = (
    mut.groupby("sample", as_index=False)
    .agg(
        tissue=("tissue", "first"),
        treatment = ("treatment", "first"),
        dose = ("dose_mgkg", "first"),
    )
) #extracts first sample tissue, treatment and dose value for sample

var_counts = (
    mut[mut["is_mutation"]].pivot_table(index="sample", columns= "var_class", values= "variation_type", aggfunc="size", fill_value=0).reset_index()
) #counts mutation classes i.e. SNV, MNV, deletion etc
subtype_counts = (mut[mut["variation_type"].eq("snv")].pivot_table(index="sample", columns= "subtype", values= "variation_type", aggfunc="size", fill_value=0).reset_index())
#counts SNV mutations by subtype

base_summary = base_summary.merge(meta_data, how="left", on="sample")
base_summary = base_summary.merge(var_counts, how="left", on="sample")
base_summary = base_summary.merge(subtype_counts, how="left", on="sample")
#merges all dataframes for downstream analysis

base_summary.to_csv("data_set_mf.tsv", sep="\t", index=False)
#saves sample summaries to tsv