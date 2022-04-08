-------------------------------------------------------------------------------------------
  ## IRanges: Basic Usage
-------------------------------------------------------------------------------------------
  
  
#-----Installing or loading IRanges package-----
if (!require("BiocManager"))
  + install.packages("BiocManager")
BiocManager::install("IRanges")
#Load the IRanges Package
library(IRanges)


##-----Constructing IRanges----- 
#Construct 2 IRanges vectors 'ir' using the function:
ir1 <- IRanges(start = c(1,3,5), end = c(3,5,7)) 
ir1 #print #Note that you provided the start and end argument but a width estimation will be returned on its own
#more on metadata columns later

ir2 <- IRanges(start = c(1,3,5), width = 3) #any 2 arguments out of start,end,width work fine
ir2 #print


##-----Naming/Labeling the ranges-----
names(ir1) <- paste("A", 1:3, sep = "")
ir1
names(ir2) <- paste("B", 1:3, sep = "")
ir2

##-----Accessing information/components from IRanges----- 
#---Using colnames---
start(ir1)

#---Access rows/an individual range---
ir1[1]
ir1["A3"]

#---Editing information in the column accessed---
width(ir2) <- 2 #Resizes the width, changing the width of ir2 from 3 to 2 #the resize function later(more flexible)
ir2

dim(ir2) #NULL because ir1 is a vector so it has length not dimension
length(ir2) #Note that this indicates the no. of vertical columns not the width of each range


##-----Concatenate IRanges-----
c(ir1,ir2)


##----------What is a Normal IRange??----------

##-----Plotting an IRange-----
plotRanges <- function(x, xlim = x, main = deparse(substitute(x)), col = "black", sep = 0.5, ...) 
{height <- 1
if (is(xlim, "Ranges"))
  xlim <- c(min(start(xlim)), max(end(xlim)))
bins <- disjointBins(IRanges(start(x),end(x) + 1))
plot.new()
plot.window(xlim, c(0, max(bins)*(height + sep)))
ybottom <- bins * (sep + height) - height
rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col = col, ...)
title(main) 
axis(1)
}
par(mfrow = c(1,1)) #(1,1) if margins don't allow c(2,1)

ir <- IRanges(start = c(1,3,7,9), end = c(4,4,8,10))
plotRanges(ir) #gives a plot showing 2 IRanges (on the left) that are overlapping and the other two non overlapping IRanges. There was no requirement for these intervals to be non-overlapping.
#The plot can be thought of like Exons in the genome

##-----Plotting a "NORMAL" IRange-----

#---reduce()---
#created by the reduce() function.
#It's a minimal representation of the original IRange as a : set. (each integer that belongs to one or more IRanges, belong to a set)
plotRanges(reduce(ir))

#---disjoin()---
#Inverse a reduce() and disjoin the non-overlapping ranges of ir, use disjoin()
plotRanges(disjoin(ir))

##-----Manipulating IRanges-----
#This will take all original ranges and produce a single changed range

##---Resize---
ir
#---Resizing to a width of 1 where the start is fixed---
resize(ir, width = 1, fix = "start")

#---Resizing to a width of 1 where the center is fixed---
resize(ir, width = 3, fix = "center")

##---Union or Intersection of IRanges---
#because IRanges are just a set of integers, arithmatic manipulation should hold too. Just convert them to normal IRanges first.

#construct ir1, ir2
ir1 <- IRanges(start = c(1,3,5), width = 1)
ir2 <- IRanges(start = c(4,5,6), width = 1)
ir1
ir2
#---union/add ---
union(ir1, ir2) #union works like concatenating followed by reduction aka reduce()

#---reducing---
reduce(c(ir1, ir2))

#---intersect---
intersect(ir1, ir2)

##-----Finding Overlaps-----
#great power of IRanges Library.
#allows us to relate 2 sets of IRanges to each other

#WHAT ARE THE OVERLAPS/WHERE ARE THEY LOCATED
#constructing 2 IRanges 
ir1 <- IRanges(start = c(1,4,8), end = c(3,7,10))
ir2 <- IRanges(start = c(3,4), width = 3)
ir1
ir2

#---finding overlaps--- between ir1 and ir2
ov <- findOverlaps(ir1, ir2) #contents of the function will be treated like "subject vs query" i.e. "ir1 vs ir2" to match against each ither
ov #gives a 2D matrix with query hits and subject hits(????) 
#ov will basically be an adjacency matrix of indices of different overlaps
#ov above shows 3 overlap cases exist when ir1 compared to ir2

#---reading the OV output---
#---Read the OV like---
#1 of ir1 (query) overlaps with 1 of ir2 (subject)
#2 of ir1 (query) also overlaps with 1 of (subject)
#2 of ir1 (query) also overlaps with 2 of (subject)

#The OV shows that range 1 in query overlaps range 1 in subject. Can verify this by printing ir1[1] and ir2[1]
ir1[1]
ir2[1] #the ranges indeed overlap thrice


#---Access the Query hits and Subject hits-----
#by the following accessor functions
queryHits(ov) #range no. 2 of query overlaps multiple ranges of subject(2 and 3), so 2 comes twice in what is returned here

#calling out only the unique hits below
unique(queryHits(ov))

#---args(findOverlaps)---
#the many arguments of this complicated function i.e. findOverlaps
args(findOverlaps) #more uses later


#OVERLAPPING HOW OFTEN---HOW MANY OVERLAPS ARE THERE (Interested in number)
#this is more useful when we are looking at large ranges
countOverlaps(ir1, ir2) #results reads as: range1 of ir1 overlaps range 1 of ir2, element 2 of ir1 overlaps 2 elements of ir2, element 3 of ir1 overlaps nothing.

##-----FIND NEIGHBOURING RANGES-----
ir1
ir2
#what ranges in ir2 are closer to the ones in ir1?
#useful when say, you have a peak/region of interest and you want to find the nearest gene
nearest(ir1, ir2) 

##-----session info-----
sessionInfo()
