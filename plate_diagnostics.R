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

regexp <- ".*wellsat.*"                           # file glob to select data set if > 1 in a directory

script <- "~/git/single-cell-sequencing/plate_diagnostics_functions.R"

source(script)

### ---- Configurable stuff: ----
### defaults. If file plate_diagnostics_config.R is found in config.dir, it is read, which can be used
### to override the defaults

inputdir <- "."
outputdir <- inputdir

## specify the location of your empty wells (follows primer number order):
## empties <- c(357:360,381:384)        # bottomright 8 wells
empties <- seq(24,384, by=24)        # rightmost column
names(empties) <- wellname(empties)
wells <- 1:384
names(wells) <- wellname(wells)
non.empties <- setdiff(wells, empties)
names(non.empties) <- wellname(non.empties)

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
paths<-list() # dir + basename, no extension
rc<-list() # list of .coutc files (raw reads)
bc<-list() # list of .coutb files (umis)
tc<-list() # list of .coutt files (transcripts)
saturations <- list() # table with number of SAM lines, and number of different genes seen so far
stats<-list() # list of stats (written as comment on top to the *.cout?.csv files)

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

dir <- getwd()

basenames <- get.file.names(dir = dir, regexp=regexp)
## i.e. filename without directory nor extension
if (is.null(basenames) || length(basenames)==0)
  stop("Found no files called ", sprintf("%s/%s.cout(t|b|c.csv", dir, prefix))

if (length(basenames) > 1)
  warning("Found more than one set of files", 
          paste(collapse=" ", sprintf("%s/%s.*", dir,names)))

for(i in 1:length(basenames)) {
  path <- paste(dir, basenames[i], sep="/")
  warning("Reading ", path, ".* ...")
  paths[[i]] <- path
  counts.file <- paste0(path,".coutc.csv")
  umi.file <- paste0(path,".coutb.csv")
  txpt.file <- paste0(path,".coutt.csv")
  rc[[i]] <- read.counts(file=counts.file)
  bc[[i]] <- read.counts(file=umi.file)
  tc[[i]] <- read.counts(file=txpt.file)
  stats[[i]] <- read.stats(file=counts.file)
  saturations[[i]] <- read.saturations(path)
  cat("library",path,"contains a total of",nrow(tc[[i]]),"genes\n")
}

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
  file <- paste0(paths[[i]],"_plate_diagnostics.pdf")
  if(landscape.mode) { 
    pdf(file=file, width=A4.height, height=A4.width)
    par(mfrow = c(3,4)) # specify grid for plots on the pdf
  } else {
    pdf(file=file, width=A4.width, height=A4.height)
    par(mfrow = c(4,3)) # specify grid for plots on the pdf
  }

  genes <- rmspike(rc[[i]])
  ngenes <- nrow(genes)
  spikes <- keepspike(rc[[i]])
  nspikes <- nrow(spikes)
  totvalid <- sum(rc[[i]])
  totumis <- sum(bc[[i]])
  st <- apply(stats[[i]], 1, sum)

  totals <- c(ngenes=ngenes, nspikes=nspikes, reads=totvalid,
              unmapped=unname(st["unmapped"]),
              umis=sum(bc[[i]]), txpts=sum(tc[[i]]))
  infobox(script=script,dir=inputdir,filename=paths[[i]],totals=totals)

  rawreads.total <- colSums(rmspike(rc[[i]]))
  genes <- rmspike(tc[[i]])
  gene.total <- colSums(genes) # sum of unique txpts after removing spike ins

  totaltxpts(gene.total, plotmethod = "hist", emptywells=empties) # total txpts/cell

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
              emptywells=empties)

  mapped <- apply(genes,2,sum)
  unmapped <- as.matrix(stats[[i]])['unmapped',] ## unmapped stats are combined for ERCC's and genes of course
  perc <- 100*(mapped/(mapped+unmapped))

  failed <- failedwells(spikes)         #names, not logical idx!
  
  plate.plot(data=perc, main='% mapped raw gene reads',
             ticks=seq(0,100,10),
             scale.name="%",
             emptywells=empties
             ## failedwells=failed
             )

  plate.plot(data=log10(gene.total), main='gene txpts',
              ticks=0:5,
              scale.name="log10",
              emptywells=empties,
              hilite=(gene.total>=1000),
              mtext=">=1000 unique txpts")

  plate.plot(data=log10(complexity), main='complexity',
              ticks=0:4,
              scale.name="log10",
              emptywells=empties,
              hilite=complexity>=1000,
              mtext=">=1000 unique genes")

  plate.plot(data=log10(spike.total), main='ERCC txpts',
              ticks=-1:4,
              scale.name="log10",
              emptywells=empties,
              hilite=spike.total>100,
              mtext=">100 ERCCs")

  hilite <- (spike.total/gene.total>0.05)
  hilite[is.na(hilite)] <- FALSE
  plate.plot(data=log2(spike.total/gene.total), main='ratio ERCC/gene txpts',
              ticks=-4:4,
              scale.name="log2",
              emptywells=empties,
              hilite=hilite,
              mtext="ERCC/gene > 0.05")
  
  well.coverage(main="gene txpt coverage (non-empty wells)", gene.total[-empties])

###   ## saturation plots, ignoring the wells
###   rug <- with(saturations[[i]]$all, nmapped[ reads %% 1e6 <= 1 ])
###   with(saturations[[i]]$all, 
###        saturation.plot(main="gene saturation",x=nmapped, y=genes, rug=rug, maxn=NULL,
###                        xlab="mapped reads seen", ylab="unique genes seen"))
### 
###   ### note: alternative is to use same but using 'txpt(s)' instead of 'umi(s)'
###   with(saturations[[i]]$all, 
###        saturation.plot(main="umi saturation", x=nmapped, y=umis, rug=rug, maxn=totals['umis'],
###                        xlab="mapped reads seen", ylab="unique umis seen"))
### 
###   with(saturations[[i]]$all,
###        saturation.plot(main="saturation of mean umis/gene", x=nmapped, y=umis/genes,rug=rug,
###                        xlab="mapped reads seen", ylab="unique umis/unique gene seen"))
### 
  ## ---- wellwise ----
  percentile <- 95
  lastrow <- nrow(saturations[[i]]$genes_perwell)
  maxg <- quantile(probs=percentile/100, unlist(saturations[[i]]$genes_perwell[lastrow, names(non.empties)]))
  maxu <- quantile(probs=percentile/100, unlist(saturations[[i]]$umis_perwell[lastrow, names(non.empties)]))

  ## genes, mean based:
  with(saturations[[i]]$wellwise_genes_mean, 
       { rug <- nmapped[ reads %% 1e6 <= 1 ]
         saturation.plot(main="mean gene saturation per well",
                         x=nmapped, y=genes_mean, rug=rug, maxn=maxg,
                         xlab="mapped reads seen", ylab="genes per well")
         })

  ## umis, mean based:
  with(saturations[[i]]$wellwise_umis_mean,
       { rug <- nmapped[ reads %% 1e6 <= 1 ]
         saturation.plot(main="mean umi saturation per well",
                         x=nmapped, y=umis_mean,rug=rug, maxn=maxu,
                         xlab="mapped reads seen", ylab="umis per well")
       })

  ## u/g, mean based:
  u.per.g <- with(saturations[[i]], wellwise_umis_mean$umis_mean/wellwise_genes_mean$genes_mean)
  with(saturations[[i]]$wellwise_umis_mean,
       { rug <- nmapped[ reads %% 1e6 <= 1 ]
         saturation.plot(main="saturation of mean umis/gene", 
                         x=nmapped, y=u.per.g, rug=rug, maxn= max(u.per.g),
                         xlab="mapped reads seen", ylab="unique umis/unique gene seen")
       })

  cellgenes(complexity,plotmethod= "combo") # plot number of detected genes/cell, can choose 4 different plot methods
  topgenes(tc[[i]])  # 2 plots: top expressed and most variable genes
  leakygenes(plate=basenames[1], data=tc[[i]], emptywells=empties)
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
