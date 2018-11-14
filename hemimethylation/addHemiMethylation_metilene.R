library("readxl")
library("openxlsx")
library("GenomicFeatures")
library("bsseq")

###Build human reference bsseq object - needed for comparison of hemimeth.
buildEmptyBsseqHg19 <- function () {
  require(bsseq)
  require(BSgenome.Hsapiens.UCSC.hg19)
  chrs <- names(Hsapiens)[1:24]
  cgs <- lapply(chrs, function(x) start(matchPattern("CG", Hsapiens[[x]])))
  cpgr <- do.call(c, lapply(1:24, function(x) GRanges(names(Hsapiens)[x], IRanges(cgs[[x]], width = 2))))
}

###Import files (bsseq object and DMRs, and bismark.cov)
import_data <- function(path_to_bsseq, path_to_dmr_file, path_to_cov, sample_name) {
  bsseq <<- readRDS(path_to_bsseq)
  bsseq <<- updateObject(bsseq, verbose = TRUE)
  filtered_DMRs_data <- data.frame(read_excel(path_to_dmr_file))
  #filtered_DMRs_data <<- sample_n(filtered_DMRs_import, 10000)
  #filtered_DMRs <<- readRDS(path_to_dmr_file)
  filtered_DMRs <<- makeGRangesFromDataFrame(filtered_DMRs_data, keep.extra.columns = TRUE)
  bismarkBSeq <<- read.bismark(files = path_to_cov, 
                               sampleNames = sample_name, 
                               rmZeroCov = FALSE, 
                               strandCollapse = FALSE, 
                               fileType = "cov", 
                               verbose = TRUE)
}

#Find overlaps between bsseq and DMR list
find_overlap <- function(){
  overlap <- findOverlaps(query = bismarkBSeq, subject = filtered_DMRs, 
                          maxgap = -1L, minoverlap = 0L, type = "within")
  
  #Merge data from bsseq and DMR lists using their index
  #Rename columns for easier identification
  overlap_df <- data.frame(overlap)
  names(overlap_df)[1] <- "hemiMethyl_index"
  names(overlap_df)[2] <- "DMR_index"

  filtered_DMRs_2 <- filtered_DMRs[to(overlap),2]
  cpgSites <- granges(bismarkBSeq[from(overlap), 1])
  
  #Create a meta list of each DML and the index from each imported list
  overlap_df <- cbind(overlap_df, filtered_DMRs_2, cpgSites)
  names(overlap_df)[10] <- "CpG_position"
  
  #Rearrange metadata
  overlap_df <<- overlap_df[c(1, 2, 3, 4, 5, 6, 7, 8, 10)]
}

#Gather hemimethylation and coverage data from bsseq object
attach_hemimethylation <- function(){
  #Merge methylation data to overlap data
  hemiMethyl <- getMeth(bismarkBSeq, region = granges(bismarkBSeq)[overlap_df[,1]], type = "raw")
  hemiMethyl <- data.frame(unlist(hemiMethyl))
  names(hemiMethyl)[1] <- "hemi_Methylation"
  
  #merge coverage data to overlap data
  hemiCov <- getCoverage(bismarkBSeq, region = granges(bismarkBSeq)[overlap_df[,1]], type = "Cov")
  hemiCov <- data.frame(unlist(hemiCov))
  names(hemiCov)[1] <- "hemi_Coverage"
  
  output <- cbind(overlap_df, hemiMethyl, hemiCov)
  
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
  final_output <<- final_output[c(3, 4, 5, 6, 8, 9, 10, 11, 16)]
  cpg_sorted <<- final_output[order(final_output$CpG_position),]
  
  ##Sperate positive and negative strands to average
  positive_strand <- final_output[final_output[,"strand"] == "+",]
  negative_strand <- final_output[final_output[,"strand"] == "-",]
  names(negative_strand)[5] <- "neg_MP"
  names(negative_strand)[7] <- "neg_methyl"
  names(negative_strand)[8] <- "neg_cov"
  negative_strand$CpG_position <- negative_strand$CpG_position - 1
  merged_strands <-  merge(x = positive_strand, y = negative_strand, by = "CpG_position")
  merged_strands <- merged_strands[c(2, 3, 4, 5, 1, 6, 7, 8, 15, 16)]
  names(merged_strands) <- c("seqnames", "start", "end", "width", 
                             "CpG_position", "mean.ERM.HERP", "pos_methyl", "pos_cov", "neg_methyl", "neg_cov")
  
  meth_diff <- cbind(merged_strands, data.frame(abs(merged_strands$neg_methyl - merged_strands$pos_methyl), 
                                                data.frame(abs(merged_strands$neg_cov - merged_strands$pos_cov))))
  names(meth_diff)[11] <- "meth_difference"
  names(meth_diff)[12] <- "cov_difference"
}

reference <- buildEmptyBsseqHg19()

import_data(path_to_bsseq = "reference_files/Komen.bsseq.rds", 
            path_to_dmr_file = "cov_vs_methyl/EmHp/EmHp_high_cov.xlsx",
            path_to_cov = "reference_files/ERMinus_HER2Plus_S3_L001_R1_001_00_bismark_bt2_pe.cat.rmdup.bismark.cov.gz", 
            sample_name = "ERplusHER2minus")

find_overlap()


attach_hemimethylation()

write.xlsx(meth_diff, "EmHp_metilene_hemi_high_cov.xlsx")


