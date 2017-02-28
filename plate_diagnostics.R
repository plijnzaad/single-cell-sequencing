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

rm(list=ls())                           #to avoid copypaste-errors
script <- "~/git/single-cell-sequencing/plate_diagnostics_functions.R"

source(script)

### ---- Configurable stuff: ----
### defaults. If file plate_diagnostics_config.R is found in config.dir, it is read, which can be used
### to override the defaults

inputdir <- "."
outputdir <- inputdir

## specify the location of your empty wells (follows primer number order):
## emptywells <- c(357:360,381:384)        # bottomright 8 wells
emptywells <- seq(24,384, by=24)        # rightmost column

landscape.mode <- TRUE

##if and where to write merged data:
output.RaceID.file <- NULL

config.file <- "plate_diagnostics_config.R"

if (file.exists(config.file)) {
  warning("Seeing file ", normalizePath(config.file), ", sourcing it\n")
  source(config.file)            # for inspiration, see config.R.example
} else
  warning("Not seeing file ", config.file, ", using default settings\n")

### ---- No configurable stuff below this line ----

#variables for the script
names<-list() # for the shortened names of the libraries
split_files<-list() # for  full names of the plates without .cout* extension
rc<-list() # list of .coutc files (raw reads)
bc<-list() # list of .coutb files (umis)
tc<-list() # list of .coutt files (transcripts)
genesseen <- list() # table with number of SAM lines, and number of different genes seen so far
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
  counts.file <- paste(files[i],".coutc.csv", sep="")
  umi.file <- paste(files[i],".coutb.csv", sep="")
  txpt.file <- paste(files[i],".coutt.csv", sep="")
  rc[[i]] <- read.counts(file=counts.file)
  bc[[i]] <- read.counts(file=umi.file)
  tc[[i]] <- read.counts(file=txpt.file)
  stats[[i]] <- read.stats(file=counts.file)
  genesseen.file <- paste0(names[[i]], "-occurrences.txt")
  if(file.exists(genesseen.file)) {
    tab <- read.table(file=genesseen.file, sep="\t",
                                 as.is=TRUE, quote="", header=TRUE,comment.char="", row.names=NULL)
    expected.cols <- c("reads", "genes", "umis") # first few columns, ignore the rest
    ## if(ncol(tab)!= length(expected.cols))
    ##  stop("File ", genesseen.file ," does not have 6 columns ", paste(expected.cols, sep=" "))
    colnames(tab)[1:length(expected.cols)] <- expected.cols
    genesseen[[i]] <- tab
  }
  cat("library",names[[i]],"contains a total of",nrow(tc[[i]]),"genes\n")
}

genes <- rmspike(rc[[1]])
ngenes <- nrow(genes)
spikes <- keepspike(rc[[1]])
nspikes <- nrow(spikes)
totvalid <- sum(rc[[1]])
st <- apply(stats[[1]], 1, sum)

# OPTIONAL: rename part of a specified plate if you have different experiments in one plate
# do this only for one of the .coutt files at the time by choosing the correct location in the .coutt list (eg:tc[[5]])
# and manually specify the names of the experiments (eg:"alpha"/"beta") and the cell labels (eg: 1:192), in total should be 384 cells.
if(FALSE)
  names(tc[[5]])<-c(paste("alpha",c(1:192),sep="_"),paste("beta",c(1:192),sep="_")) 

#### diagnostic plots####
dir.create(outputdir, showWarnings = TRUE) # directory will be made if it doesn't exist
setwd(outputdir)
for(i in 1:length(tc)) {
  A4.width <- 8.3; A4.height <- 11.7
  file <- paste0(names[[i]],"_plate_diagnostics.pdf")
  if(landscape.mode) { 
    pdf(file=file, width=A4.height, height=A4.width)
    par(mfrow = c(3,4)) # specify grid for plots on the pdf
  } else {
    pdf(file=file, width=A4.width, height=A4.height)
    par(mfrow = c(4,3)) # specify grid for plots on the pdf
  }
  
  totals <- c(ngenes=ngenes, nspikes=nspikes, reads=totvalid,
              unmapped=unname(st["unmapped"]),
              umis=sum(bc[[i]]), txpts=sum(tc[[i]]))
  infobox(script=script,dir=inputdir,filename=split_files[i],totals=totals)

  rawreads.total <- colSums(rmspike(rc[[i]]))
  genes <- rmspike(tc[[i]])
  gene.total <- colSums(genes) # sum of unique txpts after removing spike ins

  totaltxpts(gene.total, plotmethod = "hist", emptywells=emptywells) # total txpts/cell

  spike.total<-colSums(keepspike(tc[[i]]))
  complexity<-apply(genes,2,function(x)sum(x>=1))

  logdata <- log10(rawreads.total)
  f <- logdata[ logdata != -Inf ] 
  from <- floor(min(f))
  to <- ceiling(max(logdata))
  ticks <- from:to
  plate.plot(data=logdata,
              main='mapped raw gene reads',
              ticks=ticks,
              scale.name="log10",
              emptywells=emptywells)

  mapped <- apply(genes,2,sum)
  unmapped <- as.matrix(stats[[i]])['unmapped',] ## unmapped stats are combined for ERCC's and genes of course
  perc <- 100*(mapped/(mapped+unmapped))

  failed <- failedwells(spikes)         #names, not logical idx
  
  plate.plot(data=perc, main='% mapped raw gene reads',
             ticks=seq(0,100,10),
             scale.name="%",
             emptywells=emptywells
             ## failedwells=failed
             )

  plate.plot(data=log10(gene.total), main='gene txpts',
              ticks=0:5,
              scale.name="log10",
              emptywells=emptywells,
              hilite=(gene.total>=1000),
              mtext=">=1000 unique txpts")

  plate.plot(data=log10(complexity), main='complexity',
              ticks=0:4,
              scale.name="log10",
              emptywells=emptywells,
              hilite=complexity>=1000,
              mtext=">=1000 unique genes")

  plate.plot(data=log10(spike.total), main='ERCC txpts',
              ticks=-1:4,
              scale.name="log10",
              emptywells=emptywells,
              hilite=spike.total>100,
              mtext=">100 ERCCs")

  hilite <- (spike.total/gene.total>0.05)
  hilite[is.na(hilite)] <- FALSE
  plate.plot(data=log2(spike.total/gene.total), main='ratio ERCC/gene txpts',
              ticks=-4:4,
              scale.name="log2",
              emptywells=emptywells,
              hilite=hilite,
              mtext="ERCC/gene > 0.05")
  
  well.coverage(main="gene txpt coverage (non-empty wells)", gene.total[-emptywells])

  with(genesseen[[i]], 
       saturation.plot(main="gene saturation",x=reads, y=genes,
                       xlab="reads seen", ylab="unique genes seen"))

  with(genesseen[[i]], 
       saturation.plot(main="transcript saturation", x=reads, y=umis,
                       xlab="reads seen",ylab="unique transcripts seen"))

  with(genesseen[[i]],
       saturation.plot(main="saturation of mean transcripts/gene", x=reads, y=umis/genes,
                       xlab="reads seen",ylab="unique transcripts/unique gene seen"))

  cellgenes(complexity,plotmethod= "combo") # plot number of detected genes/cell, can choose 4 different plot methods
  topgenes(tc[[i]])  # 2 plots: top expressed and most variable genes
  leakygenes(data=tc[[i]], emptywells=emptywells)
  # leakygenes plots (1): number of genes and ERCC reads in the empty corner. Will give warning if a sample has more than plate average genes/5
  # (2): top expressed genes in the empty corner and % of their total reads in empty corner compared to total reads in whole plate.
  # (3): any genes that are in the top 50 genes in empty corner, but not in the top 200 genes in rest of plate (likely artifacts)
  overseq2(rc[[i]],bc[[i]]) # plot oversequencing per molecule
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
