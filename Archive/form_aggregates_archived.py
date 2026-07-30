import pandas as pd
import numpy as np
#archived version of this script, note regional MFs are calculated relative to total assessed bases unlike the updated version which is relative to regional assessed bases
mut_path = "prj00215-variant-calls-annotated.tsv"

mut = pd.read_csv(mut_path, sep="\t", low_memory=False) #loads prev annotated mut file
print(mut.head())
mut["indv_read"] = mut["end"] - mut["start"]
#used to calculate total bases read per sample (calculates length of interval)
mut["is_mutation"] = mut["variation_type"].ne("no_variant")
#used to determine if mutation is present

mut["region_categorise"] = np.select(
    [
    mut["region"].eq ("intergenic"),
    mut["region"].eq("intronic"),
    mut["region"].eq("exon_non_coding"),
    mut["region"].eq("coding"),
    ],
    ["intergenic", "intronic", "exon_non_coding", 'coding'], default = "other"
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
)
#assigns variation types to categories and if not recognised assigns other

base_summary = (mut.groupby(["sample"], as_index=False)
           .agg(
            total_read = ("indv_read", "sum"),
            mutation_event = ("is_mutation", "sum"),
            total_rows = ("sample", "size"),
)
)#creates summary by sample, calculates total bases read and mutations (whole sample summary)
region_counts = (mut[mut["is_mutation"]].pivot_table(index="sample", columns="region_categorise", values="variation_type", aggfunc="size", fill_value=0).reset_index()
                 )
base_summary = base_summary.merge(region_counts, how="left", on="sample")
#merges region count df for downstream analysis

base_summary["total_mutations_mf"] = (base_summary["mutation_event"] / base_summary["total_read"].replace(0, np.nan))
base_summary["intergenic_mutations_mf"] = (base_summary["intergenic"] / base_summary["total_read"].replace(0, np.nan))
base_summary["intronic_mutations_mf"] = (base_summary["intronic"] / base_summary["total_read"].replace(0, np.nan))
base_summary["exon_non_coding_mf"] = (base_summary["exon_non_coding"] / base_summary["total_read"].replace(0, np.nan))
base_summary["coding_mutations_mf"] = (base_summary["coding"] / base_summary["total_read"].replace(0, np.nan))
#calculates mutation frequencies (mutation event/overall assessed bases) per sample and by region
#note region-specific mutation frequencies are done relative to toal bases read here NOT regional assessed bases

base_summary["coding_mutation_total_fraction"] = (base_summary["coding"] / base_summary["mutation_event"])
#calculates coding mutation fraction (wasn't used in analysis due to low number of coding mutations)

#add columns ie dose etc and tissue
meta_data = (
    mut.groupby("sample", as_index=False)
    .agg(
        tissue=("tissue", "first"),
        treatment = ("treatment", "first"),
        dose = ("dose_mgkg", "first"),
    )
)

#additional by sample mutation type and substitution type counts
var_counts = (
    mut[mut["is_mutation"]].pivot_table(index="sample", columns= "var_class", values= "variation_type", aggfunc="size", fill_value=0).reset_index()
) #counts mutation classes i.e. SNV, MNV, deletion etc
subtype_counts = (mut[mut["variation_type"].eq("snv")].pivot_table(index="sample", columns= "subtype", values= "variation_type", aggfunc="size", fill_value=0).reset_index())
#counts SNV mutations by subtype
base_summary = base_summary.merge(meta_data, how="left", on="sample")
base_summary = base_summary.merge(var_counts, how="left", on="sample")
base_summary = base_summary.merge(subtype_counts, how="left", on="sample")
base_summary.to_csv("data_set_mf.tsv", sep="\t", index=False)
#merges dataframes into one table and saves to tsv