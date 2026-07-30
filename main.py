import re
import pyranges as pr
import numpy as np
import pandas as pd

def read_mut(mut_path):
    return pd.read_csv(mut_path, sep="\t", dtype={"sample": str}, low_memory=False)
def read_gtf(gtf_path):

    return pr.read_gtf(
        gtf_path, as_df=True
    ) #pyranges reads gtf file, returns pyranges compatible coordinates (0-based)
    #archived version did not load using pyranges and adjusted coordinates manually

def read_metadata(meta_path):
    return pd.read_csv(meta_path, sep="\t", dtype={"sample": str}, low_memory=False)
#initial functions used to read the files, taking arguments of the file paths (so can be called from GUI or other)
#these files are returned as pandas dataframes

def extract_dose(desc):
     d = (desc or "").lower() #handles capitilisation
     if "vehicle" in d: #tags vehicle control as 0 for dose
         return 0
     if "bap" in d:
         return -1 #tagged the control as -1, easiest way to differentiate from dose groups
     match = re.search(r"(\d+)\s*mg/kg", d)
     return int(match.group(1)) if match else None #tags corresponding dose (adjustments required for other formatting)

def extract_treatment(desc):
     d = (desc or "").lower()
     if "vehicle" in d: #tags vehicle control as vehicle_control
         return "vehicle_control"
     if "positive control" in d or "bap" in d: #searches for positive control or bap, returns postive_control
         return "positive_control"  #
     if "test compound" in d:
         return "test_compound"
     return None
#currently this function is used to search for the above terms in the description column and assigns labels
#future versions could use user-defined search variables

def extract_tissue(desc):
     d = (desc or "").lower()
     if "liver" in d: #similar to above function, could use user-defined search variables to handle other tissue types
         return "liver"
     if "femur" in d:
         return "femur"
     return "unknown"

def process_metadata(metadata, sample_col = "sample", description_col = "Description", import_name_col = "ImportName"):

    metadata = metadata.copy() #makes a copy of dataframe to avoid modifying original
    metadata.columns = [c.strip() for c in metadata.columns] #removes lead/trail spaces

    if sample_col not in metadata.columns:
        raise KeyError(f"Metadata file must contain a '{sample_col}' column.")
    if description_col not in metadata.columns:
        raise KeyError(f"Metadata file must contain a '{description_col}' column.")
    if import_name_col in metadata.columns and "sample_type" not in metadata.columns:
        metadata = metadata.rename(columns={import_name_col: "sample_type"})
    if "sample_type" not in metadata.columns:
        metadata["sample_type"] = metadata[description_col]
    metadata[sample_col] = metadata[sample_col].astype(str)
    metadata["dose_mgkg"] = metadata[description_col].apply(extract_dose)
    metadata["treatment"] = metadata[description_col].apply(extract_treatment)
    metadata["tissue"] = metadata["sample_type"].apply(extract_tissue)
    keep_cols = [sample_col, "dose_mgkg", "treatment", "tissue"]
    return metadata[keep_cols].rename(columns={sample_col: "sample"})
#generally, runs the above functions on metadata file (i.e. extract functions), returns df with columns sample, dose_mgkg, treatment, tissue
#throws error if sample or desc columns are missing, otherwise returns a df with the above columns

def attach_metadata(mut, metadata):
    mut = mut.copy()
    mut["sample"] = mut["sample"].astype(str)
    parsed_metadata = process_metadata(metadata)
    merged = mut.merge(parsed_metadata, on="sample", how="left", validate="many_to_one")
    missing_samples = sorted(merged.loc[merged["tissue"].isna(), "sample"].dropna().unique())
    if missing_samples:
        print(
            "Warning: metadata was not found for "
            f"{len(missing_samples)} sample(s): {missing_samples[:10]}"
        )
    return merged
#merges extracted metadata with mut, returns mut with additional columns (dose, treatment, tissue) (mut = variant calls)

def mut_processing (mut):
    mut["start"] = pd.to_numeric(mut["start"], errors="coerce")
    mut["end"] = pd.to_numeric(mut["end"], errors="coerce")
    mut["contig"] = mut["contig"].astype(str)
    mut = mut.dropna(subset=["contig", "start", "end"]).copy()
    mut["start"] = mut["start"].astype(int)
    mut["end"] = mut["end"].astype(int)
    mut = mut.reset_index(drop=True)
    return mut
#cleans up mut dataframe, ensures start/end columns are numeric, drops missing or non-num coor rows and resets index

def gtf_intervals(gtf:pd.DataFrame, features):
    if isinstance(features, str):
        features = [features]
    out = gtf.loc[gtf["Feature"].isin(features), ["Chromosome", "Start", "End"]].copy()
    #type handling
    out["Chromosome"] = out["Chromosome"].astype(str)
    out["Start"] = pd.to_numeric(out["Start"], errors = "coerce")
    out["End"] = pd.to_numeric(out["End"], errors = "coerce")
    out = out.dropna(subset=["Chromosome","Start", "End"]).copy()
    out["Start"] = out["Start"].astype(int)
    out["End"] = out["End"].astype(int)
    return out
#returns dataframe ready for PyRanges :)

def annotation_mut(mut: pd.DataFrame, gtf: pd.DataFrame, output_path: str ) -> pd.DataFrame:
    mut = mut_processing(mut)

    genes = gtf_intervals(gtf, "gene")
    exons = gtf_intervals(gtf, "exon")
    coding = gtf_intervals(gtf, ["CDS", "start_codon", "stop_codon", "Selenocysteine"])
    #extracts region-associated intervals from gtf file

    mut_interval = mut[["contig", "start", "end"]].copy()
    mut_interval["id"] = np.arange(len(mut))
    mut_interval = mut_interval.rename(columns = {"contig": "Chromosome","start": "Start",  "end": "End"})

    _mut = pr.PyRanges(mut_interval)
    _gene = pr.PyRanges(genes)


    _exon = pr.PyRanges(exons)
    _coding = pr.PyRanges(coding)
    #creates PyRanges objects for each interval type

    gene_rows = set(_mut.overlap(_gene).df["id"])
    exon_rows = set(_mut.overlap(_exon).df["id"])
    coding_rows = set(_mut.overlap(_coding).df["id"]) #changed from .join to .overlap, retuns .mut rows overlapping specified feature
    #IDs of .mut rows overlapping with annotation, non overlapping are classified as intergenic (default)

    print(
        f"gene hits: {len(gene_rows)}, exon hits: {len(exon_rows)}, coding hits: {len(coding_rows)}, total rows: {len(mut)}")
    #prints number of rows overlapping each annotation type (genes, exons, coding)
    #NOTE this isn't number of mutation events, it includes no_variant and variant

    intronic = gene_rows - exon_rows  #defines rows as intronic if in genes but not inside exons
    mut["region"] = "intergenic"  # default class
    mut.loc[mut_interval["id"].isin(intronic), "region"] = "intronic"
    mut.loc[mut_interval["id"].isin(exon_rows), "region"] = "exon_non_coding"
    mut.loc[mut_interval["id"].isin(coding_rows), "region"] = "coding"
    #overwrites region in hierarchy of coding > exon_non_coding > intronic > intergenic
    print(mut["region"].value_counts(dropna=False))
    mut.to_csv(output_path, sep="\t", index=False)
    return mut
#annotates all rows in the .mut, this could be adjusted to only apply regional context to mutation events
#Future iterations could include GUI-based options to annotate only mutation events
#Since implementing this broke pivot table functionality, it was not included
def run_code(
    mut_path: str,
    gtf_path: str,
    meta_path: str,
    output_path: str = "prj00215-variant-calls-annotated.tsv"):
    mut = read_mut(mut_path)
    metadata = read_metadata(meta_path)
    gtf = read_gtf(gtf_path)
    print("GTF feature counts")
    print(gtf["Feature"].value_counts(dropna=False))
    mut = attach_metadata(mut, metadata)
    annotated = annotation_mut(mut, gtf, output_path=output_path)
    return annotated
#runs the above functions, returns annotated mut dataframe

if __name__ == "__main__":
    run_code(
        mut_path="prj00215-variant-calls.mut",
        gtf_path="gencode.vM38.basic.annotation.gtf",
        meta_path="manifest.txt",
        output_path="prj00215-variant-calls-annotated.tsv")
#uses these paths if file is run directly
#if using GUI will use gui paths