import re
import pyranges as pr
import numpy as np
import pandas as pd
#note this is an archived version, not currently in use

mut_path = "prj00215-variant-calls.mut"
gene_path = "gencode.vM38.basic.annotation.gtf"
metadata_path = "manifest.txt"
output_mut_annotated = "prj00215-variant-calls-annotated.tsv"
#file paths

# reading the files :)
mut = pd.read_csv(mut_path, sep="\t", dtype= {"sample" : str}, low_memory=False)
metadata = pd.read_csv(metadata_path, sep="\t")
gtf = pd.read_csv(
    gene_path, sep="\t", #tab-delimited
    comment="#", header=None, #ignore comments
    names = ["chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"],
    dtype={"chr" : str},
    low_memory=False
)
# Cleaning columns
metadata.columns = [c.strip() for c in metadata.columns]
metadata = metadata.rename(columns={"ImportName": "sample_type"})
#renames column

def extract_dose(desc):
    d = (desc or "").lower() #ensures lower case
    if "vehicle" in d:
        return 0 #assigns 0 for vehicle control
    if "bap" in d:
        return -1 #tagged the positive control as -1, easiest way to differentiate
    match = re.search(r"(\d+)\s*mg/kg", desc) #searches for one or more digits followed by mg/kg and handles whitespace
    return int(match.group(1)) if match else None #extracts found dose (adjustments required for other formatting)

def extract_treatment(desc):
    d = (desc or "").lower() #handles capitalisation
    if "vehicle" in d: #this is the naming convention for the vehicle control, would need to handle other names if used with other peoples code and gui
        return "vehicle_control"
    if "positive control" in d or "bap" in d: #similar here, could parse a standard variable into the function, make it callable,
        return "positive_control"  #
    if "test compound" in d:
        return "test_compound"
    return None

def extract_tissue(desc):
    d = (desc or "").lower() #ensures lower case
    if "liver" in d: #searches for liver text in parsed column and returns "liver"
        return "liver"
    if "femur" in d:
        return "femur"
    return "unknown"
#searches and returns liver/femur tissue types else returns unknown

metadata["dose_mgkg"] = metadata["Description"].apply(extract_dose) #parses the description column for dose information
metadata["treatment"] = metadata["Description"].apply(extract_treatment)
metadata["tissue"] = metadata["sample_type"].apply(extract_tissue)
#runs the above functions on specified columns (Description/sample_type)

mut = mut.merge(metadata, on="sample", how="left") #merges new columns into mut dataframe
mut = mut.drop(columns=["sample_type", "Description"]) #drops old columns


mut["start"] = pd.to_numeric(mut["start"], errors = "coerce")
mut["end"] = pd.to_numeric(mut["end"], errors = "coerce")
mut["contig"] = mut["contig"].astype(str)
#ensuring right datatypes
mut = mut.dropna(subset=["contig", "start", "end"]).copy()
#removes missing values

mut["start"] = mut["start"].astype(int)
mut["end"] = mut["end"].astype(int)
mut = mut.reset_index(drop=True)  # <-- ensure mut index is 0, 1, 2, ...

#GTF files are + 1 relative to .PyRanges
gtf["start0"] = gtf["start"]- 1 #thus, subtracting 1 from start column to match pyranges 0-based
gtf["end0"] = gtf["end"]

def gtf_intervals(gtf:pd.DataFrame, features, chr_col = "chr", start_col = "start0", end_col = "end0"):
    #returns PyRanges-ready DataFrame :)
    if isinstance(features, str):
        features = [features]
    out = gtf.loc[gtf["feature"].isin(features), [chr_col, start_col, end_col]].copy()
    out = out.rename(columns={chr_col: "Chromosome", start_col: "Start", end_col: "End"})
    #if feature occurs, selects rows with feature type and copies coordinates

    #type handling
    out["Chromosome"] = out["Chromosome"].astype(str) #ensures chromosomes are strings
    out["Start"] = pd.to_numeric(out["Start"], errors = "coerce")
    out["End"] = pd.to_numeric(out["End"], errors = "coerce") #converts coordinate starts/ends to numeric
    out = out.dropna(subset=["Chromosome","Start", "End"]).copy() #drops invalids :)
    out["Start"] = out["Start"].astype(int)
    out["End"] = out["End"].astype(int)
    return out

genes = gtf_intervals(gtf, "gene")
exons = gtf_intervals(gtf, "exon")
coding = gtf_intervals(gtf, ["CDS", "start_codon", "stop_codon", "Selenocysteine"])
#parses features into the above function, checks for instances

mut_interval = mut[["contig", "start", "end"]].copy()
mut_interval["id"] = np.arange(len(mut)) #adds an ID column for later use
mut_interval = mut_interval.rename(columns = {"contig": "Chromosome","start": "Start",  "end": "End"})
#renames columns for PyRanges

#print(mut_interval.head()) #test print

print("mut chromosomes:", mut_interval["Chromosome"].unique()[:10])
print("gene chromosomes:", genes["Chromosome"].unique()[:10]) #test print for sanity check

_mut = pr.PyRanges(mut_interval)
_gene = pr.PyRanges(genes)
_exon = pr.PyRanges(exons)
_coding = pr.PyRanges(coding)
#creates PyRanges objects for each interval for later overlap

crossref_gene = _mut.join(_gene, how="left").df
crossref_exon = _mut.join(_exon, how="left").df
crossref_coding = _mut.join(_coding, how="left").df #joins mut and pyranges objects

gene_rows = set(crossref_gene.loc[crossref_gene["Start_b"] != -1, "id"]) #!= -1 otherwise it classifies no overlap as an overlap and incorrectly classifies everything
exon_rows = set(crossref_exon.loc[crossref_exon["Start_b"] != -1, "id"])
coding_rows  = set(crossref_coding.loc[crossref_coding["Start_b"] != -1, "id"]) #

print(f"gene hits: {len(gene_rows)}, exon hits: {len(exon_rows)}, coding hits: {len(coding_rows)}, total muts: {len(mut)}")

intronic = gene_rows - exon_rows #so basically, it minuses off the IDs that are exons = resulting in intronic :)
mut["region"] = "intergenic" #default class
mut.loc[mut_interval["id"].isin(intronic), "region"] = "intronic"
mut.loc[mut_interval["id"].isin(exon_rows), "region" ]= "exon_non_coding"
mut.loc[mut_interval["id"].isin(coding_rows), "region"] = "coding"
#overwrites region in hierarchy of coding > exon_non_coding > intronic > intergenic

print(mut["region"].value_counts(dropna=False))

mut.to_csv(output_mut_annotated, sep="\t", index=False)
#saves annotated mut file
