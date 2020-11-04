library(tidyverse, quietly = TRUE)
library(argparser, quietly = TRUE)
library(here, quietly = TRUE)
source(here::here('src' , 'R', 'utils', 'IO.R'))
source(here::here('src', 'R', 'utils', 'TCRClonaltiyFxns.R'))
options(future.globals.maxSize= 100000000000)  # Override memory limit

# Parse commandline arguments
p = arg_parser("TCR sequence clonality visualization program")
p = add_argument(p, "--input", help = "Input TCR sequence file")
p = add_argument(p, "--output", help = "Path to save output files")
p = add_argument(p, "--select_pt_loc", help = "Select patients for location plots")
p = add_argument(p, "--select_pt_sev", help = "Select patients for severity plots")
argv = parse_args(p)

cat("Running Analysis...\n")
seq_data = read.csv(argv$input, stringsAsFactors = FALSE)
savedir = argv$output
create_file_directory(savedir)

TCRcols = colnames(seq_data)

seq_data = filter(seq_data, cell_type %in% c("CD4", "CD8") & location %in% c("Blood", "LM", "LAD", "LCX", "RCA", "mixed")) %>%
  mutate(cell_type=factor(cell_type, levels=c("CD4", "CD8")),
         phenotype=factor(phenotype, levels=c("Blood", "Lipid-rich", "Fibroatheroma", "Complex", "Fibrocalcific")),
         location=factor(location, levels=c("Blood", "LM", "LAD", "LCX", "RCA", "mixed")))

getPalette = colorRampPalette(c(brewer.pal(8, "Set2"),
                                brewer.pal(12, "Paired"),
                                brewer.pal(7, "Dark2")))

# Set clonality colors
# multiple clonality groups
# myClonColor = data.frame(clon_cut=c(1, 2, 5, 10, 20, 30),
#                          clonality=c("Unique", ">= 2", ">= 5", ">= 10", ">= 20", ">= 30"),
#                          clon_color=c("#EFEFEF", "yellow", "red", "darkblue", "brown", "purple"),
#                          stringsAsFactors=FALSE)
# Clonal vs non clonal
myClonColor = data.frame(clon_cut=c(1, 2),
                         clonality=c("Unique", ">= 2"),
                         clon_color=c("#EFEFEF", "yellow"),
                         stringsAsFactors=FALSE)

# Set plotting parameters
paras = data.frame(cellType=c("CD8", "CD4"), 
                   brch2=c("ab", "ab"), 
                   brch1=c("b", "b"),
                   branch2=c("alpha-beta", "alpha-beta"), 
                   branch1=c("beta", "beta"))

# Select patients
if (!is.na(argv$select_pt_loc)) {
  select_patients_location = Map(trimws, strsplit(argv$select_pt_loc, ',')[[1]])
  other_patients_location = setdiff(unique(seq_data$masked_ID), select_patients_location)
} else {
  select_patients_location = NA
  other_patients_location = unique(seq_data$masked_ID)
}
if (!is.na(argv$select_pt_sev)) {
  select_patients_severity = Map(trimws, strsplit(argv$select_pt_sev, ',')[[1]])
  other_patients_severity = setdiff(unique(seq_data$masked_ID), select_patients_severity)
} else {
  select_patients_severity = NA
  other_patients_severity = unique(seq_data$masked_ID)
}

# Generate Plots
for(i in 1:nrow(paras)) {
  # CDR3a/b pair sequences
  ab_subset = filter(seq_data, cell_type==paras$cellType[i] & !is.na(CDR3a.g) & !is.na(CDR3b.d))
  if (nrow(ab_subset)==0) next
  
  # CDR3ab across artery location
  gridPieBar(filter(ab_subset[, TCRcols], location %in% c("LM", "LAD", "LCX", "RCA")),
             masked_ID, location, CDR3_pair,
             paste0(paras$cellType[i], paras$brch2[i], "_pie_PtLoc_c"),
             paras$branch2[i], B.D_variable, B.D_joining, A.G_variable, A.G_joining,
             ovlpColor=list(2,3), filCol1s=other_patients_location)
  
  # CDR3ab across artery location - select patients
  if (!is.na(select_patients_location)) {
    gridPieBar(filter(ab_subset[, TCRcols], location %in% c("LM", "LAD", "LCX", "RCA")),
               masked_ID, location, CDR3_pair,
               paste0(paras$cellType[i], paras$brch2[i], "_pie_PtLoc_selected_c"),
               paras$branch2[i], B.D_variable, B.D_joining, A.G_variable, A.G_joining,
               ovlpColor=list(2,3), filCol1s=select_patients_location)
  }
  
  # CDR3ab across artery location & Blood
  gridPieBar(filter(ab_subset[, TCRcols], location != "mixed"),
             masked_ID, location, CDR3_pair,
             paste0(paras$cellType[i], paras$brch2[i], "_pie_PtLoc_Blood_c"),
             paras$branch2[i], B.D_variable, B.D_joining, A.G_variable, A.G_joining,
             ovlpColor=list(2,3), filCol1s=other_patients_location)
  
  # CDR3ab across artery location & Blood - select patients
  if (!is.na(select_patients_location)) {
    gridPieBar(filter(ab_subset[, TCRcols], location != "mixed"),
               masked_ID, location, CDR3_pair,
               paste0(paras$cellType[i], paras$brch2[i], "_pie_PtLoc_Blood_selected_c"),
               paras$branch2[i], B.D_variable, B.D_joining, A.G_variable, A.G_joining,
               ovlpColor=list(2,3), filCol1s=select_patients_location)
  }
  
  # CDR3ab across phenotype
  gridPieBar(filter(ab_subset[, TCRcols], phenotype %in% c("Lipid-rich", "Fibroatheroma", "Complex", "Fibrocalcific")),
             masked_ID, phenotype, CDR3_pair,
             paste0(paras$cellType[i], paras$brch2[i], "_pie_Ptphenotype_c"),
             paras$branch2[i], B.D_variable, B.D_joining, A.G_variable, A.G_joining,
             ovlpColor=list(2,3), filCol1s=other_patients_severity)
  
  # CDR3ab across phenotype - select patients
  if (!is.na(select_patients_severity)) {
    gridPieBar(filter(ab_subset[, TCRcols], phenotype %in% c("Lipid-rich", "Fibroatheroma", "Complex", "Fibrocalcific")),
               masked_ID, phenotype, CDR3_pair,
               paste0(paras$cellType[i], paras$brch2[i], "_pie_Ptphenotype_selected_c"),
               paras$branch2[i], B.D_variable, B.D_joining, A.G_variable, A.G_joining,
               ovlpColor=list(2,3), filCol1s=select_patients_severity)
  }
  
  # CDR3b sequences
  b_subset = filter(seq_data, cell_type==paras$cellType[i] & !is.na(CDR3b.d))
  if (nrow(b_subset)==0) next
  
  # CDR3b across artery location
  gridPieBar(filter(b_subset[, TCRcols], location %in% c("LM", "LAD", "LCX", "RCA")),
             masked_ID, location, CDR3b.d,
             paste0(paras$cellType[i], paras$brch1[i], "_pie_PtLoc_c"),
             paras$branch1[i], B.D_variable, B.D_joining,
             ovlpColor=list(2,3), filCol1s=other_patients_location)
  
  # CDR3b across artery location - select patients
  if (!is.na(select_patients_location)) {
    gridPieBar(filter(b_subset[, TCRcols], location %in% c("LM", "LAD", "LCX", "RCA")),
               masked_ID, location, CDR3b.d,
               paste0(paras$cellType[i], paras$brch1[i], "_pie_PtLoc_selected_c"),
               paras$branch1[i], B.D_variable, B.D_joining,
               ovlpColor=list(2,3), filCol1s=select_patients_location)
  }
  
  # CDR3b across artery locations & Blood
  gridPieBar(filter(b_subset[, TCRcols], location != "mixed"),
             masked_ID, location, CDR3b.d,
             paste0(paras$cellType[i], paras$brch1[i], "_pie_PtLoc_Blood_c"),
             paras$branch1[i], B.D_variable, B.D_joining,
             ovlpColor=list(2,3), filCol1s=other_patients_location)
  
  # CDR3b across artery locations & Blood - select patients
  if (!is.na(select_patients_location)) {
    gridPieBar(filter(b_subset[, TCRcols], location != "mixed"),
               masked_ID, location, CDR3b.d,
               paste0(paras$cellType[i], paras$brch1[i], "_pie_PtLoc_Blood_selected_c"),
               paras$branch1[i], B.D_variable, B.D_joining,
               ovlpColor=list(2,3), filCol1s=select_patients_location)
  }
  
  # CDR3b across phenotype
  gridPieBar(filter(b_subset[, TCRcols], phenotype %in% c("Lipid-rich", "Fibroatheroma", "Complex", "Fibrocalcific")),
             masked_ID, phenotype, CDR3b.d,
             paste0(paras$cellType[i], paras$brch1[i], "_pie_Ptphenotype_c"),
             paras$branch1[i], B.D_variable, B.D_joining,
             ovlpColor=list(2,3), filCol1s=other_patients_severity)
  
  # CDR3b across phenotype - select patients
  if (!is.na(select_patients_severity)) {
    gridPieBar(filter(b_subset[, TCRcols], phenotype %in% c("Lipid-rich", "Fibroatheroma", "Complex", "Fibrocalcific")),
               masked_ID, phenotype, CDR3b.d,
               paste0(paras$cellType[i], paras$brch1[i], "_pie_Ptphenotype_selected_c"),
               paras$branch1[i], B.D_variable, B.D_joining,
               ovlpColor=list(2,3), filCol1s=select_patients_severity)
  }
}

warnings() # See warnings
cat("Analysis Finished\n")
cat("All output is saved in: ", savedir, "\n")