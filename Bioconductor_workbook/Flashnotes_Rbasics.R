#Do i have bioconductor already installed?
BiocManager::available()

#I didn't have it installed already. So I used the code below
#Install Bioconductor using BiocManager instead of BiocLite since version>3.7
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

#Update all/some/none? [a/s/n]
a


