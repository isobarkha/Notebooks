#Do i have bioconductor already installed?
BiocManager::available()

#I didn't have it installed already. So I used the code below
#Install Bioconductor using BiocManager instead of BiocLite since version>3.7
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

BiocManager::install("Biostrings")

BiocManager::install("IRanges")


#Load the IRanges Package
library(IRanges)

##Constructing IRanges vector 'ir' using the function:
ir1 <- IRanges(start = c(1,3,5), end = c(3,5,7)) 
ir1 #Note that you provided the start and end argument but a width estimation will be returned on its own
#more on metadata columns later


ir2 <- IRanges(start = c(1,3,5), width = 3) #any 2 arguments out of start,end,width work fine
ir2

##Naming/Labelling the ranges
names(ir1) <- paste("A", 1:3, sep = "")
ir1
names(ir2) <- paste("B", 1:3, sep = "")
ir2

##Accessing information/components from IRanges 
#Using colnames
start(ir1)

#Access rows/an individual range
ir1[1]
ir1["A3"]

#Editing information in the column accessed
width(ir2) <- 2 #Resizes the width, changing the width of ir2 from 3 to 2 #the resize function later(more flexible)
ir2

dim(ir2) #NULL because ir1 is a vector so it has length not dimension
length(ir2) #Note that this indicates the no. of vertical columns not the width of each range


##Concatenate IRanges
c(ir1,ir2)


#!!!!!!!!!!!!!!!!!!!!!#What is a Normal IRange?? mentions around 2:30


##Plotting an IRange
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
par(mfrow = c(2,1))
ir <- IRanges(start = c(1,3,7,9), end = c(4,4,8,10))
plotRanges(ir)
