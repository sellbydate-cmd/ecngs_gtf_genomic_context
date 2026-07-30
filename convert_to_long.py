import pandas as pd
import numpy as np


df = pd.read_csv("dataset_cleaned.tsv", sep="\t") #uses pyrimidine converted dataset (conducted in excel prior)
df.columns = df.columns.str.strip()
#reads tab delimited file and removes white space from columns

ids = ["sample", "tissue", "treatment", "dose"]

count_cols = ["intergenic_mutation_event", "intronic_mutation_event", "exon_non_coding_mutation_event", "coding_mutation_event"]
base_cols = ["intergenic_assessed_bases", "intronic_assessed_bases","exon_non_coding_assessed_bases","coding_assessed_bases"]
#for genomic region dataset

counts = df.melt( #.melt essentially just converts to long format from wide
    id_vars=ids, #but does not change these columns in the resulting format
    value_vars=count_cols,
    var_name="region", #new column with region names for long format
    value_name="mutation_count" #new column with mutation counts for long format
)

counts["region"] = counts["region"].str.replace(
    "_mutation_event", "", regex=False
) #removes _mutation_event from region names for long format

bases = df.melt(
    id_vars=ids,
    value_vars=base_cols,
    var_name="region",
    value_name="assessed_bases"
) #same as above, but for assessed bases, converting to long format

bases["region"] = bases["region"].str.replace(
    "_assessed_bases", "", regex=False
) #removes _assessed_bases from region names for long format

regional_long = counts.merge(
    bases,
    on=ids + ["region"]
) #merges the two melted dataframes on the id columns and region column

regional_long["log_assessed_bases"] = np.log(regional_long["assessed_bases"]) #takes natural log of assesed bases

regional_long = regional_long.sort_values(["sample", "region"]).reset_index(drop=True)
#sorts data by sample and region

base_substitutions = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]

base_sub_long = df.melt(
    id_vars=ids + ["total_assessed_bases"],
    value_vars=base_substitutions,
    var_name="base_substitution",
    value_name="mutation_count"
) #melts yet again, but this time for base substitutions (i.e. converting to long format)

base_sub_long["log_total_assessed_bases"] = np.log(base_sub_long["total_assessed_bases"]) #natural log of total assessed bases
base_sub_long = base_sub_long.sort_values(["sample", "base_substitution"]).reset_index(drop=True) #sorts data by sample and base substitution

subtypes = ["complex", "indel", "mnv", "snv", "symbolic"]

subtype_long = df.melt(id_vars=ids + ["total_assessed_bases"], value_vars=subtypes, var_name="mutation_subtype", value_name="mutation_count")
subtype_long["log_total_assessed_bases"] = np.log(subtype_long["total_assessed_bases"]) #natural log of total assessed bases for model offset
subtype_long = subtype_long.sort_values(["sample", "mutation_subtype"]).reset_index(drop=True) #sorts data by sample and mutation subtype

regional_long.to_csv("regional_mutation_counts_long.tsv", sep="\t", index=False)
base_sub_long.to_csv("base_substitution_counts_long.tsv", sep="\t", index=False)
subtype_long.to_csv("mutation_subtype_counts_long.tsv", sep="\t", index=False)
#saves the resulting long format dataframes to tsv files for use in modeling

print(regional_long.head(10))
print(base_sub_long.head(10))
print(subtype_long.head(10)) #prints first 10 rows of each dataframe to ensure conversion was as expected