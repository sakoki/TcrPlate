# Original Author: Houyin
# Modified & adapted by: Koki Sasagawa

library(RColorBrewer)
library(ggpubr)
library(ggrepel)
library(tidyverse)
library(scales)
library(car)

# ---------- Clonality calculation and graphing functions ----------
CalculateClonality = function(subSeqData, col1, col2, colseq, clonColor, ...) {
  #' Calculates the clonality of CDR3 sequences by counting the number of times it occurs after grouping the data by specified columns
  #' @param subSeqData TCR clonality dataframe 
  #' @param col1 first column to group data 
  #' @param col2 second column to group data 
  #' @param colseq third column to group data (CDR3 a.a. sequence) 
  #' @param clonColor dataframe assigning color based on clonality 
  TCRclon = group_by(subSeqData, !!enquo(col1), !!enquo(col2), !!enquo(colseq), !!!enquos(...)) %>%
    summarise(clon=n()) %>% as.data.frame() %>%
    mutate(clonality="Unique")
  for (i in 2:nrow(clonColor)) {
    TCRclon = mutate(TCRclon, clonality=replace(clonality, clon >= clonColor$clon_cut[i], clonColor$clonality[i]))
  }
  TCRclon = mutate(TCRclon, clonality=factor(clonality, levels=clonColor$clonality))
  return(TCRclon)
}


CalculateOverlap = function(TCRclon, col1, colseq, col2) {
  #' Count the number of locations a CDR3 sequence appears after grouping the data
  #' @param TCRclon TCR clonality dataframe 
  #' @param col1 first column to group data 
  #' @param colseq second column to group data (CDR3 a.a. sequence) 
  #' @param col2 third column to group data 
  ovlpSeq = group_by(TCRclon, !!enquo(col1), !!enquo(colseq), !!enquo(col2)) %>%
    summarise(duplSeq=n()) %>% # duplicated seqs with different variables (such as Va, Ja, ...)
    summarise(nOvlp=n(), ovlpCols=paste(!!enquo(col2), collapse=";")) %>%
    as.data.frame()
  return(ovlpSeq)
}


FilterOvlpSeq = function(ovlpSeq, noContain, knOvlpEqual) {
  #' Filters data based on specified paramters 
  #' @param noContain if specified, removes all sequences that appear in this location (e.g., 'LM' --> remove all sequences found in LM)
  #' @param knOvlpEqual if specified, keeps all sequences overlapping in exactly knOvlpEqual locations (e.g., '2' --> keep all sequences found in two locations)
  if (!is.null(noContain)) {
    ovlpSeq = filter(ovlpSeq, !grepl(noContain, ovlpCols))
  }
  if (!is.null(knOvlpEqual)) {
    ovlpSeq = filter(ovlpSeq, nOvlp==knOvlpEqual)
  }
  return(ovlpSeq)
}


MergeClonOvlp = function(TCRclon, ovlpSeq, allx) {
  return(merge(x=TCRclon, y=ovlpSeq, sort=F, all.x=allx))
}


FilterCol1 = function(TCRclon, col1, filCol1s) {
  #' Filters the data by removing entries containing any of the target labels in a specified column 
  #' @param TCRclon TCR clonality dataframe 
  #' @param col1 target column to filter by
  #' @param filCol1s vector containing target labels to look for in col1 
  if (!is.null(filCol1s)) {
    return(filter(TCRclon, !!enquo(col1) %in% filCol1s))
  } else {
    return(TCRclon)
  }
}


CalculateTCRclonSZ = function(TCRclon, ...) {
  #' Sums together the clonality counts in any specified groupings 
  TCRclon_SZ = group_by(TCRclon, !!!enquos(...)) %>%
    summarise(totClon=sum(clon)) %>%
    as.data.frame()
  return(TCRclon_SZ)
}


TCRclon.xyPercfillC = function(TCRclonOv, TCRclon_SZ, col1, col2, colseq) {
  #' Calculates the values to fill the gridPieBar
  #' @param TCRclonOv TCR clonality dataframe 
  #' @param TCRclon_SZ The total number of clonal sequnces by grouping (as calculated by CalculateTCRclonSZ)
  #' @param col1 first column to group by
  #' @param col2 second column to group by 
  #' @param colseq (CDR3 a.a. sequence) 
  TCRclon = merge(TCRclonOv, TCRclon_SZ, sort=FALSE) %>% ungroup() %>%
    group_by(!!enquo(col1), !!enquo(col2)) %>%
    arrange(clon, !!enquo(colseq), .by_group=TRUE) %>%
    mutate(x=cumsum(clon)-clon/2,
           y=nOvlp) %>%
    mutate(x=x/totClon,
           perc=clon/totClon,
           angle=-360*x,
           fillC=!!enquo(colseq),
           motif=NA,
           motif=replace(motif, grepl("QEG", fillC), "QEG"),
           motif=replace(motif, grepl("SQE", fillC), "SQE"),
           motif=replace(motif, grepl("SQEG", fillC), "SQEG")) 
  return(TCRclon)
}


colorSeq = function(TCRclon, ovlpColor, ncol2) {
  #' Assign color to fill the pie chart 
  #' @param TCRclon TCR clonality dataframe 
  #' @param ovlpColor list of pair numbers where first number is the index speciying the item in a switch() function. 
  #' Second number is optional parameter for an expression in the list 
  TCRclon = mutate(TCRclon, fillC=replace(fillC, 
                                          switch(ovlpColor[[1]],
                                                 nOvlp!=ncol2,                                     # list(1):           case1, overlap all ncol2 rows
                                                 nOvlp<ovlpColor[[2]],                             # list(2, nn):       case2, overlap at least nn (2 for default) rows
                                                 nOvlp==1 | !(grepl(ovlpColor[[2]], ovlpCols)),    # list(3, "Blood"):  case3, overlap at least two rows and must have "Blood" row
                                                 ovlpCols!=ovlpColor[[2]],                         # list(4, "Plaque"): case4, seqs enriched in "Plaque" row
                                                 nOvlp!=1),                                        # list(5):           case5, non overlapped seqs (enriched in each row)
                                          "white")) %>% ungroup()
  return(TCRclon)
}


# Used in gridPieBar to specify row and column variables of the figure
SetFacetVars = function(ncol1, ncol2, col1, col2) {
  if(ncol1 >= ncol2) {
    varRow = vars(!!enquo(col2))
    varCol = vars(!!enquo(col1))
  } else {
    varRow = vars(!!enquo(col1))
    varCol = vars(!!enquo(col2))
  }
  return(list(varRow, varCol))
}


ClonalityPerc = function(TCRclon, col1, col2, y0) {
  clonPercGrp = group_by(TCRclon, !!enquo(col1), !!enquo(col2), clonality) %>%
    summarise(x=weighted.mean(x, clon),
              y=y0,
              perc=sum(perc),
              angle=-360*x)
  return(clonPercGrp)
}


ExtractClon = function(clonPercGrp, col1, col2) {
  clonPercG = filter(clonPercGrp, clonality != "Unique") %>%
    ungroup() %>%
    group_by(!!enquo(col2), !!enquo(col1)) %>%
    summarise(perc=sum(perc)) %>%
    mutate(seqType="EnrichedExpanded")
  return(clonPercG)
}

# ---------- Main Function ----------
gridPieBar = function(subSeqData, col1, col2, colseq, nmStr, branch, ..., 
                      clonColor=myClonColor, ovlpColor=list(2,2), 
                      noContain=NULL, knOvlpEqual=NULL, allx=FALSE,
                      filCol1s=NULL) {
  # ---------- Clonality Calculation ----------
  TCRclon = CalculateClonality(subSeqData, !!enquo(col1), !!enquo(col2), !!enquo(colseq), clonColor, !!!enquos(...))
  ovlpSeq = CalculateOverlap(TCRclon, !!enquo(col1), !!enquo(colseq), !!enquo(col2))
  ovlpSeq = FilterOvlpSeq(ovlpSeq, noContain, knOvlpEqual)
  TCRclonOv = MergeClonOvlp(TCRclon, ovlpSeq, allx)
  TCRclonOv = FilterCol1(TCRclonOv, !!enquo(col1), filCol1s)
  TCRclon_SZ = CalculateTCRclonSZ(TCRclonOv, !!enquo(col1), !!enquo(col2))
  TCRclonOv = TCRclon.xyPercfillC(TCRclonOv, TCRclon_SZ, !!enquo(col1), !!enquo(col2), !!enquo(colseq))
  ncol1 = length(unique(TCRclonOv[[deparse(substitute(col1))]]))
  ncol2 = length(unique(TCRclonOv[[deparse(substitute(col2))]]))
  TCRclonOv = colorSeq(TCRclonOv, ovlpColor, ncol2)
  write.csv(TCRclonOv, file=paste0(savedir, nmStr, ovlpColor[[1]], "-indiv.csv"))
  facetVars = SetFacetVars(ncol1, ncol2, !!enquo(col1), !!enquo(col2))
  
  # -------- New clonality plot ----------
  # highlighted in grey - sequences that do not meet the overlap threshold to be colored, but overlap in at least one other location 
  tmp = mutate(TCRclonOv, fillC=replace(fillC, clon >= 1 & nOvlp > 1 & fillC == "white", "grey20"))
  
  if ("grey20" %in% tmp$fillC){
    ovlpSeqColor = c(getPalette(length(unique(tmp$fillC))-2), "grey20", "white")
  } else {
    ovlpSeqColor = c(getPalette(length(unique(tmp$fillC))-1), "white")
  }
  
  f1 = ggplot() + theme_bw(base_size=20) +
    geom_bar(data=tmp, mapping=aes(x=x, y=length(unique(tmp[[2]])), fill=fillC), color="grey80", width=tmp$perc*0.99, stat="identity", position="dodge", size=0.1) +
    scale_fill_manual(values=ovlpSeqColor) +
    geom_text(data=TCRclon_SZ, mapping=aes(x=0, y=-4, label=paste0("N=", totClon)), size=6) +
    facet_grid(rows=facetVars[[1]], cols=facetVars[[2]]) +
    coord_polar() +
    labs(x=branch, y="", fill="color codes for sequences") +
    theme(axis.text=element_blank(),
          axis.ticks=element_blank(),
          plot.title=element_text(color="blue", size=14, face="bold", hjust=0.5),
          panel.grid.minor=element_blank(),
          # panel.grid.major=element_line(color="black", size=0.1, linetype="dashed"),
          panel.grid  = element_blank(),
          panel.border = element_blank(),
          legend.position="bottom",
          legend.direction = "vertical",
          legend.key.size=unit(.3, "cm"),
          legend.text=element_text(size=5.5)) +
    scale_x_continuous(breaks=c()) +
    guides(fill=guide_legend(ncol=5, byrow=TRUE)) +
    scale_y_continuous(breaks=seq(0, ncol2, 2), limits=c(-4, ncol2), expand=c(0,0))
  ggsave(filename=paste0(savedir, nmStr, ovlpColor[[1]], "-Fig1.pdf"), f1, width=min(45, 3*ncol1), height=22, dpi=300, units="in")

  # ---------- Clonality Plots within locations ----------
  clonPercGrp = ClonalityPerc(TCRclonOv, !!enquo(col1), !!enquo(col2), ncol2+0.5)
  clonGrpColor = clonColor$clon_color[match(sort(unique(TCRclonOv$clonality)), clonColor$clonality)]
  f2 = ggplot() + theme_bw(base_size=20) +
    geom_bar(data=TCRclonOv, mapping=aes(x=x, y=1.2, fill=clonality), width=TCRclonOv$perc, stat="identity", position="dodge") +
    scale_fill_manual(values=clonGrpColor) +
    geom_text(data=clonPercGrp, mapping=aes(x=x, y=1.3, label=ifelse(100 > round(perc*100,1) & round(perc*100,1) > 5.0, paste0(round(perc*100,1), "%"), ''), vjust=0, angle=angle)) +
    geom_label(data=TCRclon_SZ, mapping=aes(x=0.5, y=0.7, label=paste0("N=", totClon)), vjust=0.5, hjust=0.5) +
    facet_grid(rows=facetVars[[1]], cols=facetVars[[2]]) +
    coord_polar() +
    labs(x=branch, y="") +
    theme(axis.text=element_blank(),
          axis.ticks=element_blank(),
          panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          panel.border = element_blank()) +
    scale_y_continuous(breaks=seq(0, 1.3, 1), limits=c(0, 1.3))
  ggsave(filename=paste0(savedir, nmStr, ovlpColor[[1]], "-Fig2.pdf"), f2, width=min(45, 3*ncol1), height=17, dpi=300, units="in")
}



