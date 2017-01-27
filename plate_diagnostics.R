# Script for making a pdf containing diagnostic plots of SORT-Seq data. It goes trough the following steps:
# optional: transform CS1 style data to CS2-style layout (from 4 96-cell libraries to 1 384-cell library)
# read in ALL .cout* files in one folder
# label cells with unique name taken from the input .cout* files and a number for each cell
# optional: rename part of specific plates in the case of two different experiments in one plate
# plot one pdf per plate containing diagnostic plots
# merge specified plates into one file that can be loaded into RaceID

## questions can be addressed to m.muraro@hubrecht.eu (and/or plijnzaad@gmail.com)

library(RColorBrewer)
library(oce)

script <- "~/git/single-cell-sequencing/plate_diagnostics_functions.R"

source(script)

### ---- Configurable stuff: ----
### defaults (override them in config.R if needed)

inputdir <- "."
outputdir <- "."

landscape.mode <- TRUE

## specify the location of your empty wells (follows primer number order):
emptywells <- c(357:360,381:384)

##if and where to write merged data:
output.RaceID.file <- NULL

config.dir <- '.'
config.file <- paste0(config.dir, "/plate_diagnostics_config.R")

source(config.file) # for inspiration, see config.R.example

### ---- No configurable stuff below this line ----

#variables for the script
names<-list() # for the shortened names of the libraries
split_files<-list() # for  full names of the plates without .cout* extension
rc<-list() # list of .coutc files (raw reads)
bc<-list() # list of .coutb files (umis)
tc<-list() # list of .coutt files (transcripts)
stats<-list() # list of stats (written as comment on top to the *.coutt.csv file)

#### OPTIONAL: Merge CS1-type libraries into one 384-column long CS2-type file####
# do this ONLY if you have CS1-style (4 libraries in one plate) data that you want to run this script on
# this assumes the order of A1=lib1, A2=lib2, B1=lib3, B2=lib4 when it reorders the columns
# setwd("/Users/mauro/AvO_lab/R/GK_GateID/raw_files/GK2E_GK2F/") # set directory containing the files
# CS1files <- read_files(dir = getwd()) # lists all the unique library names in the working directory
# input required is 4 unique library handles in the right order, followed by the name of the new, merged file, like so:
# reorder.cs1(c("EP10-12-1_AHVK72BGXX_S1","EP10-12-2_AHVK72BGXX_S2","EP10-12-3_AHVK72BGXX_S3","EP10-12-4_AHVK72BGXX_S4"),"EP10-CS2-live")

####read files and label cells ####
#read data into four lists

setwd(inputdir)
files <- read_files(dir = getwd())

# read in files
split_files <-sub(".*\\/","",files) 
for(i in 1:length(files)){
  names[[i]] <-  sub("\\_.*","",split_files[[i]]) # split lib name to keep only name supplied by you to cuppen group 
  cat("\n",split_files[[i]],"was renamed to",names[[i]],"\n",sep = " ") # fyi how the libraries will be named
  cat("reading .cout files for plate",i, "out of", length(files),"\n",sep = " ") # reports progress
  counts <- paste(files[i],".coutc.csv", sep="")
  rc[[i]] <- read.csv(file=counts, header = TRUE, sep = "\t",row.names =1, comment.char="#")
  stats[[i]] <- read.stats(file=counts)
  bc[[i]] <- read.csv(paste(files[i],".coutb.csv", sep=""), header = TRUE, sep = "\t",row.names =1, comment.char="#")
  tc[[i]] <- read.csv(paste(files[i],".coutt.csv", sep=""), header = TRUE, sep = "\t",row.names =1, comment.char="#")
  cat("library",names[[i]],"contains a total of",nrow(tc[[i]]),"genes\n")
}

ngenes <- nrow(rmspike(rc[[1]]))
nspikes <- nrow(keepspike(rc[[1]]))
valid <- sum(rc[[1]])
st <- apply(stats[[1]], 1, sum)

#label cells: all cells in library will get a _1 to _384 extension to the library name specified in names object
for(i in 1:length(tc)){ colnames(tc[[i]])<-paste(names[[i]],c(1:384),sep="_") }

# OPTIONAL: rename part of a specified plate if you have different experiments in one plate
# do this only for one of the .coutt files at the time by choosing the correct location in the .coutt list (eg:tc[[5]])
# and manually specify the names of the experiments (eg:"alpha"/"beta") and the cell labels (eg: 1:192), in total should be 384 cells.
if(FALSE)
  names(tc[[5]])<-c(paste("alpha",c(1:192),sep="_"),paste("beta",c(1:192),sep="_")) 

#### diagnostic plots####
dir.create(outputdir, showWarnings = TRUE) # directory will be made if it doesn't exist
setwd(outputdir)
for(i in 1:length(tc)){
  A4.width <- 8.3; A4.height <- 11.7
  file <- paste0(names[[i]],"_plate_diagnostics.pdf")
  if(landscape.mode) { 
    pdf(file=file, width=A4.height, height=A4.width)
    par(mfrow = c(3,4)) # specify grid for plots on the pdf
  } else {
    pdf(file=file, width=A4.width, height=A4.height)
    par(mfrow = c(4,3)) # specify grid for plots on the pdf
  }
  
  totals <- c(ngenes=ngenes, nspikes=nspikes, reads=valid, unmapped=unname(st["unmapped"]), umis=sum(bc[[i]]), txpts=sum(tc[[i]]))
  infobox(script=script,dir=inputdir,filename=split_files[i],totals=totals)
  totalreads(tc[[i]],plotmethod = "hist", emptywells=emptywells) # plots total UMI reads/cell, can choose 4 different plot methods
  cellgenes(tc[[i]],plotmethod= "cumulative") # plot number of detected genes/cell, can choose 4 different plot methods
  overseq2(rc[[i]],bc[[i]]) # plot oversequencing per molecule
  plate.plots(tc[[i]]) # 3 plots: total reads, ERCC reads and division between the two over a plate layout
  topgenes(tc[[i]])  # 2 plots: top expressed and most variable genes
  leakygenes(data=tc[[i]], emptywells=emptywells) # NB: this function can give errors if you don't have ERCCs or very few succesfully sequenced cells. comment out in case of errors
  # leakygenes plots (1): number of genes and ERCC reads in the empty corner. Will give warning if a sample has more than plate average genes/5
  # (2): top expressed genes in the empty corner and % of their total reads in empty corner compared to total reads in whole plate.
  # (3): any genes that are in the top 50 genes in empty corner, but not in the top 200 genes in rest of plate (likely artifacts)
  dev.off()
} #make pdf with diagnostic plots

if (!is.null(output.RaceID.file)) { 
####merge and write multiple dataframes into one .csv####
    cdata_all<-tc[[1]] # should be first position you want to start merging from
    for (i in 2:length(tc)){
        cdata_all <- intersectmatrix(cdata_all,tc[[i]]) #  JC's function intersectmatrix
    } # specify second position to last position you want to merge from-to
    cdata_all<- cdata_all[order(rownames(cdata_all)), ] #make row names alphabetical
    write.table(cdata_all,file=output.RaceID.file, sep="\t") # this file can be loaded into RaceID
}

# Local variables:
# mode: R
# ess-indent-level: 2
# End:
