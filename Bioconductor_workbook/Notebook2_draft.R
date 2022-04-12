#-------------------------------------------------------------------------------------------
## GenomicRanges: Basic Usage
#-------------------------------------------------------------------------------------------


##-----Installing or loading GenomicRanges package-----
if (!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("GenomicRanges")

#Load the GenomicRanges Package
library(GenomicRanges)

##-----Creating a GRange----- using GRanges(). 
gr = GRanges(seqnames = c("chr1"), strand = c("+", "-", "+"), ranges = IRanges(start = c(1,3,5), width = 3)) #chromosome name = seqnames, strand = DNA Strand identity
gr #returns 3 sequences for two different strands of the same chromosome
#chromosome strands can be of 3 types: "+", "-", "*" <* indicates unknown strand/entity present on both strand irrespective of direction> 

##-----What is the flanking sequence on the Grange?-----  
flank(gr, 5) #to get a GRanges object containing the ranges that include the 5 bases upstream of the ranges
#region that flanks would be relative to direction of transcription (explains the + and negative values returned to tell if the region is flanking from the left or the right side of my GRange)

##-----What would be the Promoter on this range------
promoters(gr)
#this gives a 2200 base pair interval where 200 bases are downstream and 2000 bases are upstream of the transcription start site

##-----What is known of my chromosome of interest-----
seqinfo(gr) #Not much. Let's fill more info and check later. Note that 'my chromosome' = chromosome of interest being investigated at all times

##-----Adding information about my chromosome-----
##--- ---
##--- ---
##--- ---