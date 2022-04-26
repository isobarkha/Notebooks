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

#---calling a column to modify---
gr$score_modified = score/3 #Error in score/3 : non-numeric argument to binary operator
gr$score_modified = gr$score/3 #gr$score_modified 
gr$score_modified #could perform vector operations now.

gr #print gr #Now the GRange shows both score and an additional metacolumn of "score_modified"


#-Construct gr2-
gr2 <- GRanges(seqnames = c("chr1", "chr2", "chr1"), strand = "*", ranges = IRanges(start= c(1,3,5), width = 3))
gr2

##-----Finding Overlaps-----
#The main work of GRanges class is this function

#---Find overlaps gr VS gr2---
findOverlaps(gr,gr2)
#-How to read this overlap result?-
# > 1. There are 4 overlaps. 
# > 2. gr[1] overlaps with gr2 [1]; gr[2] overlaps with gr2[1] and gr2[3]; gr[3] overlaps with gr2[3]
# > 3. The star strands of gr2 overlapped with both + and - of gr (as they should). NOTE: If strand direction is known (aka no '*' strands) and there is strand incompatibility, but you still want to see if + will overlap with -, use (ignore.strand = TRUE) as shown below. This will completely ignore strand direction and allow all strands to behave like '*' i.e. where positive strands overlap with both positive and negative strands i.e. irrespective of the direction. By doing this, there will be no case where where we don't have overlap because of strand incompatibility.

#---Select/Subset by Overlaps---
#To only select ranges that overlap with some range

subsetByOverlaps(gr, gr2) #3 overlaps returned
#versus
subsetByOverlaps(gr2, gr) #2 overlaps returned because gr2[2] doesn't overlap with anything as the query


##-----makeGRangesfromDataFrame
#for when you find odd looking 'data.frames' that seem to be containing data values like a GRange, so you want to convert them into proper-looking GRanges

#---step1: make data.frame
df = data.frame(chr = "chr1", start = 1:3, end = 4:6, score = rnorm(3))
df

#---step2: data.frame -> GRanges
makeGRangesFromDataFrame(df) #note that the SCORE column gets dropped as default

#---step3: retain the extra columns lost during conversion; DataFrame made with GRange info given to this
makeGRangesFromDataFrame(df, keep.extra.columns = TRUE) #returns DataFrame with an intact metacolumn of SCORE

