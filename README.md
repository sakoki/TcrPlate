# TCR Plate Data Analysis

Folder structure

    TcrPlate
    ├── data
    │  ├── fullTCRseq.csv  # Plate data sequence
    │  └── output          # Save output
    ├── R
    │  ├── TCR_clonality_figures.R
    │  └── utils
    │     ├── IO.R
    │     └── TCRClonaltiyFxns.R
    └── TCR_analysis.sh

TCR_clonality_figures.R is a analysis script that generates plots for sequence clonality by various metadata categories. 

The following arguments are accepted

    -x, --opts           RDS file containing argument values
    -i, --input          Input TCR sequence file
    -o, --output         Path to save output files
    -s, --select_pt_loc  Select patients for location plots
    --select_pt_sev      Select patients for severity plots

The script was developed with the following:

R version 4.0.1 (2020-06-06)  
Platform: x86_64-apple-darwin17.0 (64-bit)  
Running under: macOS Catalina 10.15.7  

Packages: 
1. car 3.0-10
2. scales 1.1.1
3. ggrepel 0.8.2
4. ggpubr 0.4.0
5. RColorBrewer 1.1-2
6. here 0.1
7. argparser 0.6
8. tidyverse 1.3.0   

To run TCR_clonality_figures.R in commandline:

```bash
Rscript ./R/TCR_clonality_figures.R \
--input "${inputDir}/fullTCRseq.csv" \
--output $outputDir \
--select_pt_loc "pt_id_1, pt_id_2, pt_id_3" \
--select_pt_sev "pt_id_1, pt_id_2, pt_id_3"
```
