#-------------------------------------------------------------------------------------------
## GenomicRanges: 2.0. 
#-------------------------------------------------------------------------------------------
#Load the GenomicRanges Package
BiocManager::install("GenomicRanges")
library(GenomicRanges)

##-----DataFrame-----
##---Creating DataFrames to store IRanges
ir = IRanges(start = c(1:3), width = 2)
df = DataFrame(ir = ir, score = rnorm(3))
df 

#---more structured output than data.frame(so, don't use for IRanges--
df2 = data.frame(ir = ir)
df2

##-----Subsetting-----
#---Using SelectBrackets---
df[1,1]
#---Using $---
df$ir

##-----Adding MetaColumns to DataFrame-----
#aka element Meta Data/ mcols
#GRanges can have these additional columns, unlike IRanges, called MetaColumns;
#contains additional data or other supplementary things
####-TLDR; DataFrame with additional values in it-

#---Construct a GRange---
gr <- GRanges(seqnames = "chr1", strand = c("+", "-", "+"), ranges = IRanges(start = c (1,3,5), width = 3))
gr
#---Adding additional MetaData column in DataFrame (along with a GRange)
#add MetaData column called SCORE (as in GeneScore) in gr and storing both mcol and gr in a single DataFrame 
values(gr) = DataFrame(score = rnorm(3))
gr

##-----Subsetting the metacolumn-----
#---calling only values with brackets---
values(gr) #way1
mcols(gr) #way 2
#---calling columns with $---
gr$score
