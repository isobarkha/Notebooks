#Do i have bioconductor already installed?
BiocManager::available()

#I didn't have it installed already. So I used the code below
#Install Bioconductor using BiocManager instead of BiocLite since version>3.7
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

#Update all/some/none? [a/s/n]
a

##Learning basic R codes

##----------I. VECTORS----------
#Creating vectors and checking their class type

##-----numeric-----
x <- 1:10

#print x
x

#naming vector elements of x vector
names(x) = letters[1:10]
#check class
class(x)

#subsetting the vector with 1D subsetting
x[1:3]
x[c("a","b")]

##-----uNames (unique names)-----
#if i have used non unique names for my vector elements, only the first match is returned when i seek for it with the select brackets
x = 1:3
names(x) = c("a", "a", "b")
#now when i ask my code to return what is stored as "a" (below), it gives back 1 but not 2. This approach can create problems later and with large dataset size
x["a"]

##-----uNames2-----
#to check if we have any duplicated vectors
anyDuplicated(names(x))
#returns "2" indicating 2 copies are present in total
names(x)

#to get a "0" return, change it to all unique names, reset vector names(as below)
names(x) = letters[1:3]
#check if i have duplicates now
anyDuplicated(names(x))
#Zero duplicates present now

##-----integerVSnumerics-----
x = 1
class(x) #class() shows that 1 is numeric

#if i want to assign x as an integer, put an L at the end to tell R that
x = 1L
class(x) #x is an integer now

##-----machine-----
#There's a limit to how big integers you can represent in R. Easy to run into these limits when working with genomic data with billions of base pairs
#The limit is given by a machine constant((dot)machine). you can get these limits by writing out the following
.Machine$integer.max
#rounding off to see what the limit  is roughly equal to
round(.Machine$integer.max/10^6)
#return value indicating its in billions

#The solution to running into limits is usually converting the integer into numeric as shown below
Z <- as.numeric(1L) 
class(Z)

##----------II. MATRICES----------
##-----constructing a matrix-----
x = matrix(1:9, ncol= 3, nrow = 3)
x
#note that matrices are "column-first" so 1:9 will be filled column-first-wise

#matrix can have row names and column names(which do not have to be unique,like it wasn't for a vector, but unique names help you know why)
rownames(x) = letters [1:3]
#check dimension of matrix
dim(x) #indicates 3 rows and 3 columns: R,C
ncol(x)
nrow(x)

##-----Subsetting matrix-----
#subset a matrix using 2D subsetting
x[1:2,] #gives first two rows and all columns (2 R and all C)
x[,1:2] #gives all rows and 2 columns (only the columns i asked for)
x[1:2,1:2] #gives 2 rows and 2 columns

#subset to return my desired values only
x[x>5] #subset my x where x is greater than 5

##-----1D matrix-----
x["a",] # Nope! Dimension of matrix is lost, reading the returned line is hard/might create errors. Matrices are just REALLY long vectors. Having the dimension vector makes reading so much easier. lets try to not loose the dimension below)
#say I want to look at 1D matrix (single row but all columns returned i.e. dimensions not dropped entirely), below is how I'd subset
x["a",, drop=FALSE] #YEAH! Notice the difference

##-----byrow= TRUE-----
#Since matrices are "column-first" by default, 1:9 will be filled column-first-wise
#In order to fill the matrix row-wise use
x = matrix(1:9, ncol = 3, nrow = 3, byrow = TRUE)
x

##---------- III. LIST----------
##-----createList----- 

x <- list(a = rnorm(3), b = letters[1:5], matrix) 
#matrix prints as a function #position of 3rd object of this list [[3]] and this 3rd object is a matrix

#assign names to the three objects on my list(showing above a a,b and [[3]])
names(x) <- c("numbers", "letters", "matrix function description")
x

##-----Subsetting a list----- 
x <- list(a = rnorm(3), b = letters[1:5], matrix)
#Referring/calling out to components of matrix like you would with a vector
x$a
#calling out FIRST TWO elements of list (all entries of first 2 elements of list appear)
x[1:2]
#callout all entries of second element of list(b)
x[2]
#call out the FIRST COMPONENT<Existing within list> 
x[[2]]
#DOUBLE BRACKETS CALL OUT ELEMENTS OF ELEMENTS WITHIN THE LIST
##-----SHORTFORM NAME ISSUE-----
#AVOID indexing/looking for list components by using SHORT FORMS NAMES like let(short for letter)
#R is capable of working with short forms like this but suppose you need another vector letter2 later. Then, R will not be able to decide if it needs to index letter1 or letter 2 just by reading your short form "let")
#locating a few components: examples:

names(x) = (c("a", "letters", "c"))
#callout letters using the abbreviation 'let'
x$let #assumes let is letters and returns all elements of letters
#change the list component names to have 2 kinds of "letters"
names(x) = (c("a", "letters", "letters2"))
x$let #R gets confused about whether i want letter or letter2, so it returns NULL


##-----asList----- ?????????(i dont think i entirely get the use of it)???????????
#want to print some values/information/data as a list
as.list(1:3)

##-----lapply-----
#to use the same function on each element of the list. done by taking the list and applying the intended function after a comma)
#remember each element may or may not be of the same type in list

x = list(rnorm(3), 3:9)
x
lapply(x, mean)

##-----sapply-----
#returns a simplified view of the values obtained via lapply. same usage i.e applying one function to each element of list
sapply(x, mean)
#sapply can also be done by using unlist() function which also returns simplied(un-listed) view of lapply answers
unlist(lapply(x, mean))


##----------IV. DATAFRAMES----------
#allows to hold observations of different kind. Matrices do not let you hold two variables of different types (sex=character, age=number)

##-----Creating a Dataframe-----
x <- data.frame(sex=c("M", "M", "F"), age = c(32, 34, 29))
x

##-----Subsetting/indexing/dfColumns-----
#dataframes are column based. $sign calls out columns
x$sex #returns contents of this column
x[["sex"]] 
x[1,"sex"] #returns 1st row of the column "sex"

##-----dfApply-----
sapply(x,class)

##-----as.X, error=TRUE-----
##INTER-CONVERSION OF DATATYPES
#Converting datatype into another type
x #type dataframe
as.matrix(x) #converts df to matrix
as.list(x) #converts df to list

#sometimes the conversion you want to make with biological data is not these basic commands. Install "methods" package from library in the following way and then use 'as(dot)' followed by the conversion function
library(methods)
as(x, "matrix") #put your non-conventional conversion function command between "". so when using bioconductors, instead of matrix you'll type other object functions generally used for our data)


##-----------------------------------------------------------------------------------------END OF SECTION------------------------------------------------------------------------------------------


