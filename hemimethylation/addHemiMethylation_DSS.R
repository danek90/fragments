library("openxlsx")
library("readxl")
library("GenomicFeatures")
library("bsseq")
library("ggplot2")
library("reshape2")
library("scales")
library("data.table")

###Build human reference bsseq object - needed for comparison of hemimeth.
buildEmptyBsseqHg19 <- function () {
  require(bsseq)
  require(BSgenome.Hsapiens.UCSC.hg19)
  chrs <- names(Hsapiens)[1:24]
  cgs <- lapply(chrs, function(x) start(matchPattern("CG", Hsapiens[[x]])))
  cpgr <- do.call(c, lapply(1:24, function(x) GRanges(names(Hsapiens)[x], IRanges(cgs[[x]], width = 2))))
}
reference <- buildEmptyBsseqHg19()


###Import files (bsseq object and filtered_DMRs)
###bsseqObject = the bsseq
###dmrCalls = excel file with DMRs for hemi comparison
###coverageObject = coverage object from bismark run
###sample = sample name from coverage object (in pooled runs)
import_files <- function(bsseqObject, dmrCalls, coverageObject, sample) {
  bsseq <- readRDS(bsseqObject)
  bsseq <- updateObject(bsseq, verbose = TRUE)
  filtered_DMRs <- data.frame(read_excel(dmrCalls))
  filtered_DMRs_grange <- makeGRangesFromDataFrame(filtered_DMRs, keep.extra.columns = TRUE)
  #add one to end for DSS
  end(filtered_DMRs_grange) <-  end(filtered_DMRs_grange) + 1 
  bismarkBSeq <- read.bismark(files = coverageObject, 
                             sampleNames = sample, 
                             rmZeroCov = FALSE, 
                             strandCollapse = FALSE, 
                             fileType = "cov", 
                             verbose = TRUE)

}


#Find overlaps between bsseq and DMR list
overlap <- findOverlaps(query = bismarkBSeq, subject = filtered_DMRs_grange, maxgap = -1L, minoverlap = 0L, type = "within")

#Merge data from bsseq and DMR lists using their index
#Rename columns for easier identification
overlap_df <- data.frame(overlap)
names(overlap_df)[1] <- "hemiMethyl_index"
names(overlap_df)[2] <- "DMR_index"

filtered_DMRs <- filtered_DMRs[to(overlap),4:6]
names(filtered_DMRs)[1] <- "start"
names(filtered_DMRs)[2] <- "end"
names(filtered_DMRs)[3] <- "width"
cpgSites <- granges(bismarkBSeq[from(overlap), 1])

#Create a meta list of each DML and the index from each imported list
overlap_df <- cbind(overlap_df, filtered_DMRs, cpgSites)
names(overlap_df)[7] <- "CpG_position"

#Rearrange metadata
overlap_df <- overlap_df[c(1, 2, 3, 4, 5, 6, 7, 10)]


#Gather hemimethylation data from bsseq object
#Merge methylation data to overlap data
hemiMethyl <- getMeth(bismarkBSeq, region = granges(bismarkBSeq)[overlap_df[,1]], type = "raw")
hemiMethyl <- data.frame(unlist(hemiMethyl))
names(hemiMethyl)[1] <- "hemi_Methylation"

output <- cbind(overlap_df, hemiMethyl)

output_grange <- makeGRangesFromDataFrame(output, 
                                          keep.extra.columns = FALSE,
                                          seqnames.field = "seqnames",
                                          start.field = "CpG_position",
                                          end.field = "CpG_position",
                                          strand.field = "strand")

### Find DML's in reference genome
find_start <- findOverlaps(query = output_grange, subject = reference, maxgap = -1L, minoverlap = 0L, type = "start")
find_end <- findOverlaps(query = output_grange, subject = reference, maxgap = -1L, minoverlap = 0L, type = "end")

strand(output_grange)[from(find_start)] <- "+"
strand(output_grange)[from(find_end)] <- "-"

final_output <- data.frame(output_grange)
final_output <- cbind(output, final_output)
final_output <- final_output[c(6, 3, 4, 5, 7, 9, 14)]

write.xlsx(final_output, "./output/allCPG_ERminusHER2plus_DSS.xlsx")


##Sperate positive and negative strands to average
positive_strand <- final_output[final_output[,"strand"] == "+",]
negative_strand <- final_output[final_output[,"strand"] == "-",]


filtered_tag <- data.frame(read_excel("ERMinus_HER2Plus_DSS_75percent.xlsx"))

average_pos <- aggregate(hemi_Methylation ~ start, positive_strand, mean)
names(average_pos)[2] <- "forward_strand_methylation" 

average_neg <- aggregate(hemi_Methylation ~ end, negative_strand, mean)
names(average_neg)[2] <- "reverse_strand_methylation"

hemi_merged_1 <-  merge(filtered_tag, average_pos, by.x = "start", by.y = "start")

hemi_merged_2 <-  merge(hemi_merged_1, average_neg, by.x = "end", by.y = "end")

hemi_data_out<- hemi_merged_2[c(3, 4, 5, 2, 1, 6, 7, 10, 11, 25, 26, 12, 13, 14, 15, 
                                16, 17, 18, 19, 20, 21, 22, 23, 24)]
write.xlsx(hemi_data_out, "./output/hemi_full_test.xlsx")

#Scraps
###############################################

output$strand <- ifelse(output$CpG_position == start_ref$start, "+", "*")
output_2$strand <-  ifelse(output_2$forward == 1, "+", "*")
output_2$strand <-  ifelse(output_2$reverse == 1, "-", output_2$strand)

output_3 <- within(output_2, output_2$strand[output_2$`lengthOneForward` == 1] <- "+")
output_4 <- within(output_3, output_3$strand[output_3$`lengthOneReverse` == 1] <- "-")





