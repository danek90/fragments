library("openxlsx")
library("readxl")
library("GenomicFeatures")
library("bsseq")
source("~/Documents/pipelines/R/CompareDMRs.R")

#Import files (bsseq object and list of DMRs)
bsseq <- readRDS("Komen.bsseq.rds")
DMRs <- data.frame(read_excel("ERPlus_HER2Minus_DSS_hemi.xlsx"))
DMRs_grange <- makeGRangesFromDataFrame(DMRs, keep.extra.columns = TRUE)
bismarkBSseq <- read.bismark(files = "ERPlus_HER2Minus_S1_L001_R1_001_00_bismark_bt2_pe.cat.rmdup.bismark.cov.gz", 
                             sampleNames = "ERPlus_HER2Minus", 
                             rmZeroCov = FALSE, 
                             strandCollapse = FALSE, 
                             fileType = "cov", 
                             verbose = TRUE)

#Find overlaps between bsseq and DMR list
overlap <- findOverlaps(query = bismarkBSseq, subject = DMRs_grange, maxgap = -1L, minoverlap = 0L, type = "within")

#Merge data from bsseq and DMR lists using their index
#Rename columns for easier identification
overlap_df <- data.frame(overlap)
names(overlap_df)[1] <- "hemiMethyl_index"
names(overlap_df)[2] <- "DMR_index"


filteredDMRs <- DMRs[to(overlap),2:4]
names(filteredDMRs)[1] <- "DMR_start"
names(filteredDMRs)[2] <- "DMR_end"
names(filteredDMRs)[3] <- "DMR_width"
cpgSites <- granges(bismarkBSseq[from(overlap), 1])

overlap_df <- cbind(overlap_df, filteredDMRs, cpgSites)

overlap_df <- overlap_df[c(1, 2, 3, 4, 5, 6, 7, 10)]
names(overlap_df)[7] <- "CpG_position"

#Gather hemimethylation data from bsseq object
#Merge methylation data to overlap data
hemiMethyl <- getMeth(bismarkBSseq, region = granges(bismarkBSseq)[overlap_df[,1]], type = "raw")
hemiMethyl <- data.frame(unlist(hemiMethyl))
names(hemiMethyl)[1] <- "hemi_Methylation"

output <- cbind(overlap_df, hemiMethyl)
output <- output[c(1, 2, 6, 3, 4, 5, 7, 8, 9)]

difference <- data.frame(diff(output$CpG_position))
names(difference)[1] <- "forward"
difference <- rbind(difference, 0)

reverse <- data.frame(reverse = rev(c(-diff(rev(output$CpG_position)),0)))

output_2 <- cbind(output, difference, reverse)
cpgOnlyOutput <- output_2[c(3, 6, 7, 8)]
output_grange <- makeGRangesFromDataFrame(cpgOnlyOutput, 
                                          keep.extra.columns = TRUE,
                                          seqnames.field = "seqnames",
                                          start.field = "CpG_position",
                                          end.field = "CpG_position",
                                          strand.field = "strand")

find_start <- findOverlaps(query = output_grange, subject = reference, maxgap = -1L, minoverlap = 0L, type = "start")
find_end <- findOverlaps(query = output_grange, subject = reference, maxgap = -1L, minoverlap = 0L, type = "end")

#Find a way to get pluses on the outfile. Look here, this part works, just need to update the output file
##############################################
start_reference <- output_grange[from(find_start)]
start_plus <- data.frame(strand_update = (matrix("+", ncol = 1, length(start_reference))))
mcols(start_reference, use.names = TRUE) <- c(start_reference, start_plus)

end_reference <- output_grange[from(find_end)]
###############################################


if (4102517 %in% start(reference)){
  output_2$strand[2,8] <- "+"
}

for (i in output_2$CpG_position){
  
  if (i %in% start(reference)){
    output_2$strand <- "+"
  } else (i %in% end(reference)){
    output_2$strand <- "-"
  }
  
}

output_2$strand <-  ifelse(output_2$forward == 1, "+", "*")
output_2$strand <-  ifelse(output_2$reverse == 1, "-", output_2$strand)

output_3 <- within(output_2, output_2$strand[output_2$`lengthOneForward` == 1] <- "+")
output_4 <- within(output_3, output_3$strand[output_3$`lengthOneReverse` == 1] <- "-")

write.xlsx(output_2, "test2.xlsx")



