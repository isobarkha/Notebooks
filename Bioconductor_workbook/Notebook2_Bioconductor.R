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
##---Define Length---
seqlengths(gr) = c("chr1" = 10)
seqinfo(gr)
seqlevels(gr)
#since a length end of chromosome is also defined now, we can talk of attribtes that are relative to chromosome as a whole

#--- Adding more chromosomes---
seqnames(gr) = c("chr1", "chr2", "chr1") #gives an error because it has recorded that this hypothetical organism only has a single chromosome <Error in .normalize_seqnames_replacement_value(value, x) : levels of supplied 'seqnames' must be identical to 'seqlevels(x)'>
#Resolve this by either redefining a new object from start OR BETTER, by defining new seqlevels as given below
seqlevels(gr) = c("chr1", "chr2") #Now it nknows that seqnames can take 2 different values/levels
seqnames(gr) = c("chr1", "chr2", "chr1")
gr

#--- Is my genome circular ---
isCircular(gr) = c(FALSE, TRUE) 
seqinfo(gr)

#--- Assigning a genome to this GRange ---
genome(gr) = "hg19" #nice info label to have on our genome since sequence/DNA might even be coming from different organisms or like a control/test sitation later. Protects from any errors later
seqinfo(gr)

#---Add a new genome---
gr2 = gr #This made a copy of gr and stored it as gr2 which can be modified as below to make changes on desired positions on gr
genome(gr2) = "hg18" #so gr2 is a copy of gr but I changed what info this genome is from
seqinfo(gr2)

##-----FindOverlaps between GRanges-----
findOverlaps(gr,gr2) #Now this gives an error saying "hey! we tried to compare both genes against each other but you're trying to compare 2 genes of different chromosomes so we can't find "overlaps". This is not a sequence similarity run like BLAST.

##-----View in Asc/Descending order-----
sort(gr)

##-----Gaps on my chromosome-----
#Gaps is a function that will give all the area on chromosome that is not covered by ranges in the GRanges
gaps(gr)


