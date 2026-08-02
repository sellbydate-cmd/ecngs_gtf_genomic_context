#################################
Requirements: 
#################################
Context assignment: 
1. Python 3.12.3
2. Windows subsystem for Linux (WSL: Ubuntu)
3. packages pyranges, numpy, pandas, tkinter

Downstream analysis:

4. R 4.3.2 
5. packages readxl, ggplot2, glmmTMB, emmeans

#################################
Data availability:
#################################

Note: This Git repository does not contain any of the input files required to reproduce the analysis.

Files required:

1. ecNGS manifest
2. ecNGS variant call file (.mut)
3. GTF file

Please contact Simon et al. (2025) for the ecNGS file and manifest.

GTF Mouse annotation file used in this project was sourced from GENCODE (GRCm39 release m38, basic GTF annotation)

#################################
Running the analysis:
#################################

Context assignment:
1. Initially run either 'gui_.py' or 'main.py'
This extracts metadata from the manifest and uses the GTF file to assign genomic context to the variant call file.
2. Run 'form_aggregates.py'
This creates a sample level summary of the annotated ecNGS file, used in downstream analysis.

Downstream analysis:
3. Complete excel pre-processing (listed below in downstream analysis script notes)
4. Run 'convert_to_long.py'
5. Run 'negbinom.R'

#################################
Archived files: 
#################################

The 'Archive' folder contains earlier working versions of the code for reference only, which were subsequently updated to the current version. Please use the current versions of 'main.py' and 'form_aggregates' in the main project directory for analysis. 

##############################
Downstream analysis script notes:
##############################

All scripts in this section were conducted after initial analysis and aggregation. 

convert_to_long.py: converts the aggregated ecNGS variant call file to long format for use in generalized linear mixed models

NOTE: This script requires pyrimidine context for base substitutions (which was performed in excel prior to running the script)

i.e. C>A (+ G>T); C>G (+ G>C); C>T (+ G>A); T>A (+ A>T); T>C (+ A>G); T>G (+ A>C). 

C>A now represents C>A and G>T and so on.

Additionally, insertions and deletions were collapsed into an indel column 

Also note positive control groups were also removed in excel prior to running the script and total_read was renamed to total_assessed_bases 

Although pandas would be more suitable for this purpose, these transformations were already conducted for graphing and thus not repeated. 

Prior to using the R script, the long format files were opened in excel and saved as .xlsx files. 

negbinom.R: R script used to fit negative binomial and zero-inflated negative binomial models and post-hoc tests

ecNGS/manifest files were attained from: 

Simon, S., Jörg Schlingemann, Johnson, G., Brenneis, C., & Dieckhoff, J. (2025). Deriving safe limits for N-nitroso-bisoprolol by error-corrected next-generation sequencing (ecNGS) and benchmark dose (BMD) analysis, integrated with QM modeling and CYP-docking analysis. Archives of Toxicology, 99(10), 3935–3962. https://doi.org/10.1007/s00204-025-04103-2
