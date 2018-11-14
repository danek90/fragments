library(bsseq)
library(VennDiagram)

### findDMROverlaps
# Find overlaps between DMRs (must include # of overlapping CpGs)
# Basically uses the findOverlaps function but also takes into
# account the number of CpGs in the overlap
# 
#' @param DMR1    granges of DMRs #1
#' @param DMR2    granges of DMRs #2
#' @param bsseq1  bsseq object with to CpGs in DMRs #1 
#' @param bsseq2  bsseq object with to CpGs in DMRs #2 (optional)
#' @param minCpG  Minimum number of overlapping CpGs in both DMRs 
#'                to be called an overlap (default = 1)
#'                Must have this # of CpGs in BOTH bsseq1 & bsseq2
#' 
#' @return granges object with intersecting granges
#'         metadata columns include: from & to
#'         > these are obtained from findOverlaps()
#' 
#' @import bsseq
# 
findDMROverlaps <- function(DMR1, DMR2, bsseq1, bsseq2 = NULL, minCpG = 1) {
  # Find overlap of DMRs with at least 2 base overlap
  ovl <- findOverlaps(DMR1, DMR2, minoverlap = 2)
  DMRintersect <- pintersect(DMR1[from(ovl)], DMR2[to(ovl)])
  # Quit of no overlap
  if (length(DMRintersect) == 0) {
    values(DMRintersect) <- data.frame(from=NULL, to=NULL)
    return(DMRintersect)
  } else {
    # Check if it contains more than minCpG CpGs
    DMRCpGcount1 <- CountCpGsGranges(DMRintersect, bsseq1)
    if(!is.null(bsseq2)) {
      DMRCpGcount2 <- CountCpGsGranges(DMRintersect, bsseq2)
    } else {
      DMRCpGcount2 <- DMRCpGcount1
    }
    DMRCpGmin <- pmin(DMRCpGcount1, DMRCpGcount2)
    DMRkeep <- DMRCpGmin >= minCpG
    # Generate intersecting DMRs and annotations
    DMRintersect <- DMRintersect[as.vector(DMRkeep)]
    DMRmatch <- data.frame(from=from(ovl)[DMRkeep],
                           to=to(ovl)[DMRkeep])
    values(DMRintersect) <- DMRmatch
    return(DMRintersect)
  }
}

### CountCpGsGranges
# Count number of CpGs in granges object
#' @param granges  granges object
#' @param bsseq    bsseq object with CpGs
#' @return df of CpG count
CountCpGsGranges <- function(granges, bsseq) {
  CpGcount <- data.frame(CpGcount=countOverlaps(granges, granges(bsseq)))
  return(CpGcount)
}

### findDMRintersect
# Get a list of DMR1/DMR2 with only intersecting DMRs
#' @param DMR1          DMR1
#' @param DMR2          DMR2
#' @param DMRintersect  DMRintersect output
#' @param unique        (default=FALSE) only return unique DMRs
#'                      But DMRs will NOT match by rows
#' @return List of DMR1 DMR2 intersection
findDMRintersect <- function (DMR1, DMR2, DMRintersect, unique=FALSE) {
  if (unique == FALSE) {
    out <- list(DMR1=DMR1[DMRintersect$from],
                DMR2=DMR2[DMRintersect$to])
  } else {
    out <- list(DMR1=DMR1[unique(DMRintersect$from)],
                DMR2=DMR2[unique(DMRintersect$to)])
  }
  return(out)
}

### find3DMROverlaps
# Finding overlap for 3 DMRs
#' @param    DMR1
#' @param    DMR2
#' @param    DMR2
#' @param    bsseq
find3DMROverlaps <- function(DMR1, DMR2, DMR3, bsseq) {
  # Find all overlaps between 2 DMRs
  DMR.ovl12 <- findDMROverlaps(DMR1, DMR2, bsseq, minCpG = 3)
  DMR.gr12 <- findDMRintersect(DMR1, DMR2, DMR.ovl12)
  DMR.ovl13 <- findDMROverlaps(DMR1, DMR3, bsseq, minCpG = 3)
  DMR.gr13 <- findDMRintersect(DMR1, DMR3, DMR.ovl13)
  names(DMR.gr13) <- c("DMR1", "DMR3")
  DMR.ovl23 <- findDMROverlaps(DMR2, DMR3, bsseq, minCpG = 3)
  DMR.gr23 <- findDMRintersect(DMR2, DMR3, DMR.ovl23)
  names(DMR.gr23) <- c("DMR2", "DMR3")
  # Find overlap for all 3 DMRs
  DMR.ovl123 <- findDMROverlaps(DMR.gr12$DMR1, DMR3, bsseq, minCpG = 3)
  DMR.gr123 <- findDMRintersect(DMR.gr12$DMR1, DMR3, DMR.ovl123)
  DMR.grtemp <- findDMRintersect(DMR.gr12$DMR2, DMR3, DMR.ovl123)
  DMR.gr123 <- list(DMR1=DMR.gr123$DMR1,
                    DMR2=DMR.grtemp$DMR1,
                    DMR3=DMR.gr123$DMR2)
  #
  DMR1only <- list(DMR1 = subtractGRintersect(subtractGRintersect(DMR1, DMR.gr12$DMR1), DMR.gr13$DMR1),
                   DMR1U2 = subtractGRintersect(DMR.gr12$DMR1, DMR.gr123$DMR1),
                   DMR1U3 = subtractGRintersect(DMR.gr13$DMR1, DMR.gr123$DMR1),
                   DMR1U2U3 = removeduplicatedDMRs(DMR.gr123$DMR1))
  DMR1only <- removeduplicatedDMRlist(DMR1only)
  DMR2only <- list(DMR2 = subtractGRintersect(subtractGRintersect(DMR2, DMR.gr12$DMR2), DMR.gr23$DMR2),
                   DMR1U2 = subtractGRintersect(DMR.gr12$DMR2, DMR.gr123$DMR2),
                   DMR2U3 = subtractGRintersect(DMR.gr23$DMR2, DMR.gr123$DMR2),
                   DMR1U2U3 = DMR.gr123$DMR2)
  DMR2only <- removeduplicatedDMRlist(DMR2only)
  DMR3only <- list(DMR3 = subtractGRintersect(subtractGRintersect(DMR3, DMR.gr13$DMR3), DMR.gr23$DMR3),
                   DMR1U3 = subtractGRintersect(DMR.gr13$DMR3, DMR.gr123$DMR3),
                   DMR2U3 = subtractGRintersect(DMR.gr23$DMR3, DMR.gr123$DMR3),
                   DMR1U2U3 = DMR.gr123$DMR3)
  DMR3only <- removeduplicatedDMRlist(DMR3only)
  # Count overlaps
  # DMR.ovl <- data.frame(DMR1 = length(DMR1),
  #                       DMR2 = length(DMR2),
  #                       DMR3 = length(DMR3),
  #                       DMR1U2 = countDMRoverlap(DMR.ovl12),
  #                       DMR1U3 = countDMRoverlap(DMR.ovl13), 
  #                       DMR2U3 = countDMRoverlap(DMR.ovl23), 
  #                       DMR1U2U3 = countDMRoverlap(DMR.ovl123))
  # Count overlap for drawing venn
  # This is a horrible way to do it - think of a better way?
  DMR.venn <- data.frame(DMR1 = length(DMR1only$DMR1), 
                         DMR2 = length(DMR2only$DMR2), 
                         DMR3 = length(DMR3only$DMR3),
                         DMR1U2 = floor(mean(c(length(DMR1only$DMR1U2), length(DMR2only$DMR1U2)))),
                         DMR1U3 = floor(mean(c(length(DMR1only$DMR1U3), length(DMR3only$DMR1U3)))),
                         DMR2U3 = floor(mean(c(length(DMR2only$DMR2U3), length(DMR3only$DMR2U3)))),
                         DMR1U2U3 = floor(mean(c(length(DMR1only$DMR1U2U3), length(DMR2only$DMR1U2U3), length(DMR3only$DMR1U2U3))))
                         )
  # Compile output
  DMR.list <- list(DMR1=DMR1only, DMR2=DMR2only, DMR3=DMR3only)
  output <- list(DMRcount=DMR.venn, DMRlist=DMR.list)
  return(output)
}

### find3DMRintersect
#
#' @param DMR1    DMR1
#' @param DMR2    DMR2
#' @param DMR3    DMR3
#' @param bsseq1  bsseq for DMR1
#' @param bsseq2  bsseq for DMR2
#' @param bsseq3  bsseq for DMR3
#' @return intersection for all 3 DMRs
find3DMRintersect <- function (DMR1, DMR2, DMR3, bsseq1, bsseq2, bsseq3) {
  # Find intersect between DMR1 and DMR2
  int12 <- findonlyDMRintersect(DMR1, DMR2, bsseq1, bsseq2)
  # Find intersect between DMR1(intersect with DMR2) and DMR3
  int12.3 <- findonlyDMRintersect(int12$DMR1, DMR3, bsseq1, bsseq3)
  # Find intersect between DMR2(intersect with DMR1) and DMR3
  int21.3 <- findonlyDMRintersect(int12$DMR2, DMR3, bsseq2, bsseq3)
  # Find true intersect for DMR3 (with the other 2)
  int3 <- unique(c(int12.3$DMR2, int21.3$DMR2))
  # Find true intersects for DMR1 and DMR2
  int1 <- findonlyDMRintersect(DMR1, int3, bsseq1, bsseq3)$DMR1
  int2 <- findonlyDMRintersect(DMR2, int3, bsseq2, bsseq3)$DMR1
  # Return
  return(list(DMR1 = int1, DMR2 = int2, DMR3 = int3))
}

### find4DMRintersect
#
#' @param DMR1    DMR1
#' @param DMR2    DMR2
#' @param DMR3    DMR3
#' @param DMR4    DMR4
#' @param bsseq1  bsseq for DMR1
#' @param bsseq2  bsseq for DMR2
#' @param bsseq3  bsseq for DMR3
#' @param bsseq4  bsseq for DMR4
#' @return intersection for all 4 DMRs
find4DMRintersect <- function (DMR1, DMR2, DMR3, DMR4, bsseq1, bsseq2, bsseq3, bsseq4) {
  # Find intersects for all pairs
  int12 <- findonlyDMRintersect(DMR1, DMR2, bsseq1, bsseq2)
  int13 <- findonlyDMRintersect(DMR1, DMR3, bsseq1, bsseq3)
  int14 <- findonlyDMRintersect(DMR1, DMR4, bsseq1, bsseq4)
  int23 <- findonlyDMRintersect(DMR2, DMR3, bsseq2, bsseq3)
  int24 <- findonlyDMRintersect(DMR2, DMR4, bsseq2, bsseq4)
  int34 <- findonlyDMRintersect(DMR3, DMR4, bsseq3, bsseq4)
  # Find only DMRs that intersect with all other sets
  int1 <- int12$DMR1[int12$DMR1 %in% int13$DMR1 & int12$DMR1 %in% int14$DMR1]
  int2 <- int12$DMR2[int12$DMR2 %in% int23$DMR1 & int12$DMR2 %in% int24$DMR1]
  int3 <- int13$DMR2[int13$DMR2 %in% int23$DMR2 & int13$DMR2 %in% int34$DMR1]
  int4 <- int14$DMR2[int14$DMR2 %in% int24$DMR2 & int14$DMR2 %in% int34$DMR2]
  # Return
  return(list(DMR1 = int1, DMR2 = int2, DMR3 = int3, DMR4 = int4))
}

findonlyDMRintersect <- function (DMR1, DMR2, bsseq1, bsseq2, unique=TRUE) {
  overlap <- findDMROverlaps(DMR1, DMR2, bsseq1, bsseq2)
  intersect <- findDMRintersect(DMR1, DMR2, overlap, unique = unique)
  return(intersect)
}

### find4DMROverlaps
# ... description ...
#' @param DMRs    list of 4 DMR objects
#' @param bsseqs  list of 4 corresponding bsseq objects
#' @return list of: 
#'         DMRlist = List of DMRs for each section
#'                   (DMR1 includes DMR1 and intersects)
#'         DMRlistunique = List of DMRs unique to each section
#'         DMRcount = DMR counts from DMRlist
#'         DMRcountunique = DMR counts from DMRlistunique
#'         
find4DMROverlaps <- function(DMRs, bsseqs) {
  # Append values with names
  for (i in 1:length(DMRs)) {
    values(DMRs[[i]]) <- data.frame(list.name=names(DMRs[i]), values(DMRs[[i]]))
  }
  # Get intersect of everything
  DMR.gr12 <- findonlyDMRintersect(DMRs[[1]], DMRs[[2]], bsseqs[[1]], bsseqs[[2]])
  DMR.gr13 <- findonlyDMRintersect(DMRs[[1]], DMRs[[3]], bsseqs[[1]], bsseqs[[3]])
  DMR.gr14 <- findonlyDMRintersect(DMRs[[1]], DMRs[[4]], bsseqs[[1]], bsseqs[[4]])
  
  DMR.gr23 <- findonlyDMRintersect(DMRs[[2]], DMRs[[3]], bsseqs[[2]], bsseqs[[3]])
  DMR.gr24 <- findonlyDMRintersect(DMRs[[2]], DMRs[[4]], bsseqs[[2]], bsseqs[[4]])
  
  DMR.gr34 <- findonlyDMRintersect(DMRs[[3]], DMRs[[4]], bsseqs[[3]], bsseqs[[4]])
  
  DMR.gr123 <- find3DMRintersect(DMRs[[1]], DMRs[[2]], DMRs[[3]],
                                 bsseqs[[1]], bsseqs[[2]], bsseqs[[3]])
  DMR.gr124 <- find3DMRintersect(DMRs[[1]], DMRs[[2]], DMRs[[4]],
                                 bsseqs[[1]], bsseqs[[2]], bsseqs[[4]])
  DMR.gr134 <- find3DMRintersect(DMRs[[1]], DMRs[[3]], DMRs[[4]],
                                 bsseqs[[1]], bsseqs[[3]], bsseqs[[4]])
  DMR.gr234 <- find3DMRintersect(DMRs[[2]], DMRs[[3]], DMRs[[4]],
                                 bsseqs[[2]], bsseqs[[3]], bsseqs[[4]])
  DMR.gr1234 <- find4DMRintersect(DMRs[[1]], DMRs[[2]], DMRs[[3]], DMRs[[4]],
                                  bsseqs[[1]], bsseqs[[2]], bsseqs[[3]], bsseqs[[4]])
  
  # Combine everything
  DMR.list <- list(DMR1=DMRs[[1]],
                   DMR2=DMRs[[2]],
                   DMR3=DMRs[[3]],
                   DMR4=DMRs[[4]],
                   DMR12=Reduce("c", DMR.gr12),
                   DMR13=Reduce("c", DMR.gr13),
                   DMR14=Reduce("c", DMR.gr14),
                   DMR23=Reduce("c", DMR.gr23),
                   DMR24=Reduce("c", DMR.gr24),
                   DMR34=Reduce("c", DMR.gr34),
                   DMR123=Reduce("c", DMR.gr123),
                   DMR124=Reduce("c", DMR.gr124),
                   DMR134=Reduce("c", DMR.gr134),
                   DMR234=Reduce("c", DMR.gr234),
                   DMR1234=Reduce("c", DMR.gr1234))
  # Get count, adjust for duplications (crude)
  DMR.countlist <- sapply(DMR.list, length)
  for (i in 1:length(DMR.countlist)) {
    DMR.countlist[i] <- floor(DMR.countlist[i] / (nchar(names(DMR.countlist[i]))-3))
  }
  
  # Find only DMRs in each section of venn
  DMR.venn <- DMR.list
  removeCommonDMRs <- function(DMR1, DMR2) {
    if(!class(DMR2) == "list") DMR2 <- list(DMR2)
    for (i in DMR2) {
      DMR1 <- DMR1[!DMR1 %in% i]
    }
    return (DMR1)
  }
  DMR.venn$DMR1 <- removeCommonDMRs(DMR.venn$DMR1, DMR.venn[c("DMR12", "DMR13", "DMR14")])
  DMR.venn$DMR2 <- removeCommonDMRs(DMR.venn$DMR2, DMR.venn[c("DMR12", "DMR23", "DMR24")])
  DMR.venn$DMR3 <- removeCommonDMRs(DMR.venn$DMR3, DMR.venn[c("DMR13", "DMR23", "DMR34")])
  DMR.venn$DMR4 <- removeCommonDMRs(DMR.venn$DMR4, DMR.venn[c("DMR14", "DMR24", "DMR34")])
  DMR.venn$DMR12 <- removeCommonDMRs(DMR.venn$DMR12, DMR.venn[c("DMR123", "DMR124")])
  DMR.venn$DMR13 <- removeCommonDMRs(DMR.venn$DMR13, DMR.venn[c("DMR123", "DMR134")])
  DMR.venn$DMR14 <- removeCommonDMRs(DMR.venn$DMR14, DMR.venn[c("DMR124", "DMR134")])
  DMR.venn$DMR23 <- removeCommonDMRs(DMR.venn$DMR23, DMR.venn[c("DMR123", "DMR234")])
  DMR.venn$DMR24 <- removeCommonDMRs(DMR.venn$DMR24, DMR.venn[c("DMR124", "DMR234")])
  DMR.venn$DMR34 <- removeCommonDMRs(DMR.venn$DMR34, DMR.venn[c("DMR134", "DMR234")])
  DMR.venn$DMR123 <- removeCommonDMRs(DMR.venn$DMR123, DMR.venn["DMR1234"])
  DMR.venn$DMR124 <- removeCommonDMRs(DMR.venn$DMR124, DMR.venn["DMR1234"])
  DMR.venn$DMR134 <- removeCommonDMRs(DMR.venn$DMR134, DMR.venn["DMR1234"])
  DMR.venn$DMR234 <- removeCommonDMRs(DMR.venn$DMR234, DMR.venn["DMR1234"])
  # Get count, adjust for duplications (crude)
  DMR.countvenn <- sapply(DMR.venn, length)
  for (i in 1:length(DMR.countvenn)) {
    DMR.countvenn[i] <- floor(DMR.countvenn[i] / (nchar(names(DMR.countvenn[i]))-3))
  }
  return(list(DMRlist=DMR.list, 
              DMRlistunique=DMR.venn,
              DMRcount=DMR.countlist,
              DMRcountunique=DMR.countvenn))
}


### subtractGRintersect
# subtract the intersect of granges 
# (Returns gr1 without overlapping ranges)
#' @param gr1   granges 1    
#' @param gr2   granges 2
#' @return gr1 without overlap
subtractGRintersect <- function (gr1, gr2) {
  ovl <- findOverlaps(gr1, gr2)
  out <- gr1[!seq_len(length(gr1)) %in% unique(from(ovl))]
  return(out)
}

### removeduplicatedDMRs
removeduplicatedDMRs <- function(DMRs) {
  return(DMRs[!duplicated(DMRs)])
}
### removeduplicatedDMRlist
removeduplicatedDMRlist <- function(DMRlist) {
  return(lapply(DMRlist, removeduplicatedDMRs))
}

### filterDMRbasic
# DMR filter
#' Basic filter includes:
#'  * CpGcount >= 3
#'  * !seqnames(DMRs) %in% c("chrX", "chrY", "chrM")
#' roadmapFilter = TRUE will also filter
#'  * DMRs$roadmapFlagged == "N"
#'  
#'  #'  *** mean.CpGcov.all >= 3 (filter taken out due to issues)
#'  
#' @param DMRs              granges object with DMRs
#' @param beta.diff         beta difference filter (use 0-1!)
#' @param roadmapFilter     (default = TRUE)
#' @return granges object - filtered
filterDMRbasic <- function(DMRs, beta.diff = 0, roadmapFilter = TRUE) {
  # Basic filters
  DMRs <- DMRs[DMRs$CpGcount >= 3,]
  # DMRs <- DMRs[DMRs$mean.CpGcov.all >= 3,]
  DMRs <- DMRs[!seqnames(DMRs) %in% c("chrX", "chrY", "chrM"),]
  
  # Filter by p-value (metilene only)
  DMRvalues <- values(DMRs)
  if (any(grepl("MWU", names(DMRvalues)))) {
    DMRs <- DMRs[DMRvalues[,grepl("MWU", names(DMRvalues))] < 0.05]
  }
  
  # Optional roadmapFilter
  if (roadmapFilter == TRUE) DMRs <- DMRs[DMRs$roadmapFlagged == "N",]
  
  ## Filter by difference in beta
  # Remove NA (also if multiple columns match, use first one)
  DMRvalues <- values(DMRs)
  mean.diff <- DMRvalues[,grepl("mean.difference",names(DMRvalues))]
  if (class(mean.diff) == "DataFrame") mean.diff <- mean.diff[,1]
  DMRs <- DMRs[!is.na(mean.diff),]
  # Filter by beta
  DMRvalues <- values(DMRs)
  mean.diff <- DMRvalues[,grepl("mean.difference",names(DMRvalues))]
  if(max(mean.diff) > 2) beta.diff <- beta.diff * 100
  DMRs <- DMRs[abs(mean.diff) > beta.diff,]
  
  return(DMRs)
}

### filterDMRlistbasic
# Filters a DMR list
#'  
#' @param DMRs            granges object with DMRs
#' @param roadmapFilter   (default = TRUE)
#' @return granges object - filtered
filterDMRlistbasic <- function(DMRlist, beta.diff = 0, roadmapFilter = TRUE) {
  DMRlist <- lapply(DMRlist, filterDMRbasic, 
                    beta.diff = beta.diff, roadmapFilter = roadmapFilter)
  return(DMRlist)
}

### drawDMRvenn2
# Draw venn diagram of overlap between 2 DMRs
# WARNING: does NOT give true intersect.
#          because DMR intersect is NOT 1-to-1!!!
# 
#' @param DMR1           granges of DMRs #1
#' @param DMR2           granges of DMRs #2
#' @param DMRintersect   findDMROverlaps output
#' @param label1         label for DMR1
#' @param label2         label for DMR2
#' 
#' @return grobTree object of venn diagram
#' 
#' @import VennDiagram
# 
drawDMRvenn2 <- function(DMR1, DMR2, DMRintersect, label1 = "DMR1", label2 = "DMR2") {
  nDMR1 <- length(DMR1)
  nDMR2 <- length(DMR2)
  # This is NOT a true intersect, since DMR union is NOT 1 to 1
  nDMR12 <- floor(mean(c(length(unique(DMRintersect$from)),
                         length(unique(DMRintersect$to)))))
  venn.rot <- 0 + 180 * (nDMR2 > nDMR1)
  venn <- draw.pairwise.venn(nDMR1, nDMR2, nDMR12,
                             euler.d = FALSE, scaled = FALSE,
                             category = c(label1, label2),
                             rotation.degree = venn.rot)
  return(venn)
}

drawDMRvenn2_values <- function(nDMR1, nDMR2, nDMR12, label1 = "DMR1", label2 = "DMR2") {
  venn.rot <- 0 + 180 * (nDMR2 > nDMR1)
  venn <- draw.pairwise.venn(nDMR1, nDMR2, nDMR12,
                             euler.d = FALSE, scaled = FALSE,
                             category = c(label1, label2),
                             rotation.degree = venn.rot)
  return(venn)
}

### makeDMRvenntable
# Get a table (instead of venn) of overlap between 2 DMRs
# WARNING: does NOT give true intersect.
#          because DMR intersect is NOT 1-to-1!!!
#' @param DMR1           granges of DMRs #1
#' @param DMR2           granges of DMRs #2
#' @param DMRintersect   findDMROverlaps output
makeDMRvenntable <- function(DMR1, DMR2, DMRintersect) {
  nDMR1 <- length(DMR1)
  nDMR2 <- length(DMR2)
  # This is NOT a true intersect, since DMR union is NOT 1 to 1
  nDMR12 <- floor(mean(c(length(unique(DMRintersect$from)),
                         length(unique(DMRintersect$to)))))
  out <- data.frame(DMR1=nDMR1, DMR2=nDMR2, intersect=nDMR12)
  return(out)
}

### countDMRoverlap
# Get overlap count between two DMRs
# WARNING: does NOT give true intersect.
#          because DMR intersect is NOT 1-to-1!!!
#' @param DMRintersect   findDMROverlaps output
countDMRoverlap <- function(DMRintersect) {
  # This is NOT a true intersect, since DMR union is NOT 1 to 1
  nDMR12 <- floor(mean(c(length(unique(DMRintersect$from)),
                         length(unique(DMRintersect$to)))))
  return(nDMR12)
}

### standardizeDMRcolumns
# Standardize columns for metilene and DSS output
# Works if no additional columns are added in front of standard output
#' @param   DMR  granges object
#' @return  DMR - formatted
standardizeDMRcolumns <- function (DMR) {
  # Get DMR values
  vals <- values(DMR)
  # If DSS:
  if(any(grepl("score", names(vals)))) {
    newVals <- data.frame(
      vals[,1:2],
      q.value = NA,
      p.MWU. = NA,
      p.2D.KS. = NA,
      vals[,3:ncol(vals)])
  } else if (any(grepl("p.MWU", names(vals)))) {
    newVals <- data.frame(
      score = NA,
      name = NA,
      vals[,c(1,4,5,8,2,9:ncol(vals))]
    )
    names(newVals)[7] <- gsub("\\.\\.", ".", names(newVals)[7])
    names(newVals) <- gsub("\\.1", "", names(newVals))
  } else {
    message("Unrecognized DMR list type!")
    newVals <- vals
  }
  values(DMR) <- newVals
  return(DMR)
}

### crudefindoverlap
crudefindoverlap <- function(granges, chr, start, end) {
  x <- seqnames(granges) %in% chr & 
    !(end(granges) < start | start(granges) > end)
  return(which(x))
}


