#Functions needed for reading, filtering and normalizing Cel-SEQ data

#remove or keep only spike ins from given data frame
rmspike<-function(x){
  ERCCs<-grep("ERCC-",row.names(x)) # gives vector with row # of spike ins
  data<-x[-ERCCs,,drop=FALSE] # make new data frame without the specified rows
  return(data) # output new data frame
}

keepspike<-function(x){
  ERCCs<-grep("ERCC-",row.names(x)) # gives vector with row # of spike ins
  data<-x[ERCCs,,drop=FALSE] # make new data frame with only the specified rows
  return(data) # output new data frame
}

# chop of chromosome lables (Abel)
chop_chr <- function  (name, splitcharacter = "__") {   strsplit(name, splitcharacter)[[1]][1]}

#plot expression of one gene as barplot
plotgene<-function(x,n){
  barplot(as.matrix(x[grep(n,rownames(x)),]),main=n)
}

# make GENEID the rownames and remove column GENEID
mvgeneid<-function(data){
  data <- as.data.frame(data)
  rownames(data) = data[,1]
  data= data[,-1]
  return(data)
}

#reorder cells from four CS1 primer libraries into one 384 column-long library
# libraries is a vector containing the four library names, in the order A1,A2,B1,B2
reorder.cs1<-function(libraries,name){
  stop("this function may have stopped working due to explicit column naming and automatic reordering")
  tc<-list()
  rc<-list()
  bc<-list()
  for(i in 1:4){
    tc[[i]]<-  read.csv(paste(libraries[i],".coutt.csv", sep=""), header = TRUE, sep = "\t",row.names =1)
    rc[[i]] <- read.csv(paste(libraries[i],".coutc.csv", sep=""), header = TRUE, sep = "\t",row.names =1) 
    bc[[i]] <- read.csv(paste(libraries[i],".coutb.csv", sep=""), header = TRUE, sep = "\t",row.names =1)
  }
  merge.tc<-intersectmatrix(tc[[1]],intersectmatrix(tc[[2]],intersectmatrix(tc[[3]],tc[[4]])))
  merge.bc<-intersectmatrix(bc[[1]],intersectmatrix(bc[[2]],intersectmatrix(bc[[3]],bc[[4]])))
  merge.rc<-intersectmatrix(rc[[1]],intersectmatrix(rc[[2]],intersectmatrix(rc[[3]],rc[[4]])))

  order<-c(matrix(c(96*0+seq(1,96), 96*1+seq(1,96)), 2, byrow = T)) 
  order2<-c(matrix(c(96*2+seq(1,96), 96*3+seq(1,96)), 2, byrow = T)) 
  all<-c()
  for(i in 0:7){
    all<-c(all,order[(1+i*24):((i+1)*24)],order2[(1+i*24):((i+1)*24)])
  }
  merge.order.tc<-merge.tc[all]
  merge.order.bc<-merge.bc[all]
  merge.order.rc<-merge.rc[all]
  
  merge.order.tc<- merge.order.tc[order(rownames( merge.order.tc)), ]
  merge.order.bc<- merge.order.bc[order(rownames( merge.order.bc)), ]
  merge.order.rc<- merge.order.rc[order(rownames( merge.order.rc)), ]
  
  write.table(merge.order.tc,paste(name,".coutt.csv",sep=""),sep="\t")
  write.table(merge.order.bc,paste(name,".coutb.csv",sep=""),sep="\t")
  write.table(merge.order.rc,paste(name,".coutc.csv",sep=""),sep="\t")
}

#JC's merge function (produces new rows on bottom of new dataframe, so reorder rows alphabetically afterwards)
intersectmatrix<-function(x,y){
  a<-setdiff(row.names(x),row.names(y))
  b<-setdiff(row.names(y),row.names(x))
  d<-matrix(data = 0,nrow = length(a),ncol = ncol(y))
  row.names(d)<-a
  colnames(d)<-colnames(y)
  c<-matrix(data = 0,nrow = length(b),ncol = ncol(x))
  row.names(c)<-b
  colnames(c)<-colnames(x)
  y<-rbind(y,d)
  x<-rbind(x,c)
  e <- match(rownames(x), rownames(y))
  f <- cbind( x, y[e,])
  return(f)
}


#overseq2, plot oversequencing per transcript
overseq2 <- function(x,y) {
  main=paste("oversequencing")   # mixes string + name of choice
#  xlab=bquote(log[10] ~ "read counts / barcode counts")    #  subscript in string
  xlab <- "log2(read count/UMI count)"
  rc.v<- unlist(x); rc.v <- rc.v[rc.v>0]
  bc.v<- unlist(y); bc.v <- bc.v[bc.v>0]
  results<-rc.v/bc.v
  hist(log2(results),breaks=75, col="red", border=NA, main=main,xlab=xlab)
  mtext(sprintf("median: %.1f",median(rc.v/bc.v)), col="red", cex=0.6)
}

commafy <- function(x, preserve.width="common") { 
  ## commafy(12345678) => "12,345,678"
  formatC(as.integer(x), format="d", big.mark=",", preserve.width=preserve.width)
}

getversion <- function(path) {
  dir <- dirname(normalizePath(path))
  cmd <- sprintf("cd %s && git describe --match 'v[0-9]*' --tags --dirty", dir)
  res <- system(cmd, intern=TRUE)
  if( length(res)==0)
      res <- 'UNKNOWN'
  res
}                                       #getversion

infobox <- function(script, dir, filename,totals) {
  ## print version filename and some overal statistics
  plot(type="n", x=c(0,1), y=c(0,1), axes=FALSE, xlab=NA,ylab=NA, main="Information")
  dir <- dirname(normalizePath(dir))
  len <- nchar(dir); maxlen <- 30
  if(len>maxlen)
    dir <- substr(dir,start=(len-maxlen)+1, stop=len)

  unmapped <- totals['unmapped']
  valid <- totals['reads']
  totreads <- valid + unmapped #note: ignores all those without a valid CBC, and more!
  text <- c(
    paste0("script version: ", getversion(script)),
    paste0("dir: ", dir),
    paste0("file: ", filename),
    sprintf("refgenes found: %d", totals['ngenes']),
    sprintf("ECCS found: %d", totals['nspikes']),
    sprintf("total reads with CBC*: %s\n\t mapped: %s (%.1f %%) ",
            commafy(totreads), commafy(valid), 100*valid/totreads),
    paste0("total umis: ",  commafy(totals['umis'])),
    paste0("total txpts: ",  commafy(totals['txpts']))
    )

  text(pos=4, x=0, y=seq(1, 0, length.out=length(text)), labels=text)
}                                       #infobox

#plot total number of reads per sample
totalreads <- function(data,plotmethod=c("barplot","hist","cumulative","combo"),
                       emptywells=NULL) {
  col.full <- rgb(0, 0.6, 0, 0.8)
  col.empty <- rgb(0, 0, 0.6, 0.6)
  magnif.empty <- 10
  cex <- 0.6
  ## if ( ! plotmethod %in% c("barplot","hist","cumulative","combo") ) stop("invalid method")
  if(plotmethod == "hist"){
    a<-hist(plot=FALSE,x=log10(colSums(rmspike(data))),breaks=100)
    ticks <- log10(as.vector(c(1,2,5) %o% 10^(0:9)))
    last <- which(ticks > max(a$breaks))[1]
    ticks <- ticks[1:last]
    rest.args <-list(xlab="counts",ylab="frequency",main="gene transcripts/well",
                     xaxt="n",col.sub="red", breaks=a$breaks,
                     xlim=range(a$breaks), ylim=c(0, max(a$counts))) 
    if(is.null(emptywells)) { 
      do.call(hist, args=c(list(x=log10(colSums(data)),
                      col=col.full, border=NA), rest.args))
    } else {
      do.call(hist, args=c(list(x=log10(colSums(data[,-emptywells])),
                      col=col.full, border=NA), rest.args))
      do.call(hist, args=c(list(x=rep(log10(colSums(data[,emptywells])),magnif.empty),
                      col=col.empty, border=NA, add=TRUE), rest.args))

      legend(x = "topleft", legend = c("full", sprintf("empty x %d", magnif.empty)),
                              fill = c(col.full, col.empty), bty="n", cex=cex)
    }
    axis(1, at=ticks,labels=sprintf("%s",commafy(round(10^ticks))),las=3, cex.axis=cex)
    mn <- mean(colSums(data))
    md <- median(colSums(data))
    abline(v=log10(c(mn/2,md,mn)),col=c("purple", "red", "brown"), lty=2, lwd=0.5)
    text(x=min(a$breaks)+(0.5 + 0.2*c(-1,1))*diff(range(a$breaks)), y=max(a$counts),
         labels=sprintf("%s: %.0f", c("median","mean"), c(md,mn)), col=c("red","brown"),
         cex=cex)
    text(log10(mn/2),max(a$counts)-2, 'half-mean', srt=0.2, col = "purple",pos=2, cex=cex*0.66)

    return()
  }
  stop("Invalid method ", plotmethod," (older versions also had barplot, cumulative and combo, untested bu may still work)")

  if(!is.null(emptywells))
    stop("emptywells argument currently only works with plotmethod='hist'")
  if(plotmethod == "barplot"){
    b<-barplot(colSums(data),xaxt="n",xlab="cells",sub=paste("mean total read:",round(mean(colSums(data)))),main="total unique reads",col="black",border=NA) 
    axis(1,at=b,labels=c(1:length(data))) # 1=horizontal at = position of marks
    abline(h=mean(colSums(data)),col="red")
  }
  if(plotmethod == "cumulative"){
    plot(ecdf(colSums(data)),xlab="total reads",ylab="fraction",main="total unique reads",col="red",tck=1,pch=19,cex=cex,cex.axis=0.8) 
    abline(v=mean(colSums(data)/2),col="red")
    mtext(paste("mean:",round(mean(colSums(data)))," median:",round(median(colSums(data)))),side=3,col="red",cex=cex)
  }
  
  if(plotmethod == "combo"){
    a<-hist(log10(colSums(data)),breaks=100,xlab="log10(counts)",ylab="frequency",main="total unique reads",col="grey",xaxt="n",col.sub="red") 
    mtext(paste("mean:",round(mean(colSums(data)))," median:",round(median(colSums(data)))),side=3,col="red",cex=cex)
    axis(1,at=a$breaks[which(a$breaks %in% c(0,1,2,3,4,5))],labels=a$breaks[which(a$breaks %in% c(0,1,2,3,4,5))])
    abline(v=log10(mean(colSums(data))/2),col="red")
    text(log10(mean(colSums(data))/2),max(a$counts)-2, round(mean(colSums(data))/2), srt=0.2, col = "red",pos=2)
    plotInset(log10(1),max(a$counts)/4,log10(250), max(a$counts),mar=c(1,1,1,1),
              plot(ecdf(colSums(data)),pch=".",col="red",cex=cex,ylab=NA,xlab=NA,main=NA,cex.axis=0.8,xaxt="n",las=3,mgp=c(2,0.1,0),tck=1,bty="n"),
              debug = getOption("oceDebug"))
  }
}                                       #totalreads

#plot complexity (number of unique genes detected) per cell
cellgenes<-function(data,plotmethod=c("hist","cumulative","combo")) {
  cex <- 0.6
  if ( ! plotmethod %in% c("hist","cumulative","combo") ) stop("invalid plotting method")

  genes<-apply(data,2,function(x) sum(x>=1))

  if(plotmethod == "hist"){
    a<-hist(genes,breaks=100,xlab="total genes",ylab="frequency",main="unique genes detected/well",col="steelblue1",xaxt="n") 
    mtext(paste("mean:",round(mean(genes))," median:",round(median(genes))),side=3,col="red",cex=cex)
    axis(1,at=a$breaks[which(a$breaks %in% seq(0,max(a$breaks),1000))],labels=a$breaks[which(a$breaks %in% seq(0,max(a$breaks),1000))])
  }

  if(plotmethod == "cumulative"){
    plot(ecdf(genes),pch=19,cex=0.2,col="red",ylab="frequency",xlab="detected genes/cell",main="cumul. complexity",cex.axis=1,las=1,tck=1)
    mtext(paste("mean:",round(mean(genes))," median:",round(median(genes))),side=3,col="red",cex=cex)
  }

  if(plotmethod == "combo"){
    a<-hist(genes,breaks=100,xlab="log10(counts)",ylab="frequency",main="detected genes/cell",col="steelblue1",xaxt="n") 
    mtext(paste("mean:",round(mean(genes))," median:",round(median(genes))),side=3,col="red",cex=cex)
    axis(1,at=a$breaks[which(a$breaks %in% seq(0,max(a$breaks),1000))],labels=a$breaks[which(a$breaks %in% seq(0,max(a$breaks),1000))])

    plotInset(max(genes)/3,max(a$counts)/3,max(genes), max(a$counts),mar=c(1,1,1,1),
              plot(ecdf(colSums(data)),pch=19,col="red",cex=cex,ylab=NA,xlab=NA,main=NA,cex.axis=0.6,las=3),
              debug = getOption("oceDebug"))
  }
}                                       #cellgenes

coverage.plot <- function(data, type='genes') {       
  if (is.null(data)) {
    .empty.plot(main="coverage", msg="no data")
    return()
  }
  max.ntxpts <- 25000                     # or new argument transcriptome.size?
  cex <- 0.6
  
  if (type=='genes') { 
    with(data,
         plot(main="gene coverage", x=reads, y=genes, xlab="reads seen",ylab="unique genes seen",
              ylim=c(0, max.ntxpts), type="l", lwd=2, col="red", cex.axis=cex, las=2,tck=1))
    return()
  }

  if (type=='umis') { 
    with(data,
         plot(main="transcript coverage", x=reads, y=umis, xlab="reads seen",ylab="unique transcripts seen",
              type="l", lwd=2, col="red", cex.axis=cex, las=2,tck=1))
    return()
  }

  if (type=='umis.per.gene') { 
    with(data,
         plot(main="mean transcripts/gene", x=reads, y=umis/genes, xlab="reads seen",ylab="unique transcripts/unique gene seen",
              type="l", lwd=2, col="red", cex.axis=cex, las=2,tck=1))
    return()
  }
  stop("Unknown coverage.plot type, only know 'genes', 'umis', 'umis.per.gene'")
}                                       #coverage.plot

#plot ERCC reads
plotspike<-function(data){
  erccs<-data[grep("ERCC-",rownames(data)),]
    b<-barplot(colSums(erccs),main="ERCC reads",ylab="total ERCC reads",xlab="cells",col="orange",xaxt="n",border=NA)
    axis(1,at=b,labels=c(1:length(data)))
}

#plot number of available transcripts vs cutoffs of median detected transcripts
testcutoff<-function(data,n,pdf=FALSE){
  main=paste("genes cutoff test",n)
  for(l in 1:15){
    z = apply(data,1,median) > l
    if(l==1){
      rc.cutoff = z
    } else {
       rc.cutoff = cbind(rc.cutoff,z)
    }
  }
  if (pdf){
    pdf(paste(getwd(),main,".pdf",sep="")) 
    plot(apply(rc.cutoff,2,sum),ylab = "number of transcripts",col="black",
    xlab = "cutoff (mean transcript no.)",main=main,type="b",lty=2,pch=19)
    dev.off()
  }
  else{
    plot(apply(rc.cutoff,2,sum),ylab = "number of transcripts",col="black",
    xlab = "cutoff (mean transcript no.)",main=main,type="b",lty=2,pch=19)
  }    
}

.empty.plot <- function(main, msg) {
  plot(type="n", x=c(0,1), y=c(0,1), main=main, axes=FALSE,xlab="",ylab="")
  text(pos=4, x=0, y=0.5, labels=msg, cex=1)
}
   
.plot.grid <- function(main, emptywells=NULL) {
  g <- expand.grid(x = c(1:24), y = c(1:16))
  plot(g,main=main, xlab=NA, ylab=NA, type="n",
       axes=FALSE, frame.plot=TRUE)
  axis(1, at=1:24, labels=1:24, cex.axis=0.4,las=3)
  axis(2, at=1:16, labels=rev(LETTERS[1:16]), cex.axis=0.4,las=1)
  if(is.null(emptywells))
    return
  x <- range(g[emptywells,'x'])
  x[1] <- x[1]-0.45
  x[2] <- x[2]+0.45
  y <- range(g[emptywells,'y'])
  y[1] <- 16 - (y[1]-0.45) + 1
  y[2] <- 16 - (y[2]+0.45) + 1
  rect(xleft=x[1],xright=x[2], ybottom=y[1],ytop=y[2], col=NA, border="grey")
}

#plot number of total reads, ERCC-reads and genes/cell over a 384-well plate layout
plate.plots<-function(data, welltotals=NULL, emptywells=NULL, rawreads=NULL) {
  cex <- 0.6
  cex.wells <- 1.1
  genes <- rmspike(data)
  gene.total<-colSums(genes) # sum of unique reads after removing spike ins
  spike.total<-colSums(keepspike(data))
  rawreads.total <- colSums(rmspike(rawreads))

  counts.palette <- colorRampPalette(c("white","blue4"))(91)

  mar <- par()$mar
  save.mar <- mar
  mar[1] <- mar[1]*1.3                  #to get right aspect ratio for the plates
  mar[3] <- mar[3]*1.3
  par(mar=mar, xpd=TRUE)

  scale.top <- -2.5; scale.bot <- -3.5

  coordinates <- expand.grid(seq(1,24),rev(seq(1,16)))

  if(is.null(rawreads))
    .empty.plot(main="% raw gene reads", msg="stats not given")
  else {
    .plot.grid(main="raw gene reads",emptywells=emptywells)
    l <- log10(rawreads.total)
    infinite <- is.infinite(l)
    l[infinite] <- NA
    col <- colorize(x=l, counts.palette, na.col='white')
    from <- floor(min(l, na.rm=TRUE))
    to <- ceiling(max(l, na.rm=TRUE))
    ticks <- from:to
    points(coordinates,pch=19,col=col, cex=cex.wells)
    points(coordinates[infinite,],pch=4, cex=0.5*cex.wells, col='black')
    draw.color.scale(xleft=1, xright=24,ybottom=scale.bot, ytop=scale.top, at=ticks, sep=0.2,cex.axis=0.5,
                     col=counts.palette, main="log10")
  }

  if(is.null(welltotals))
    .empty.plot(main="% mapped gene reads", msg="stats not given")
  else { 
    .plot.grid(main="% mapped gene reads", emptywells=emptywells)
    perc <- with(welltotals, 100*(mapped/(mapped+unmapped)))
    col <- colorize(perc, counts.palette, na.col='white')

    ticks <- seq(0,100,10)
    points(coordinates,pch=19,col=col, cex=cex.wells)
    ## mtext(sprintf(">1500 unique reads: %.0f%%",sum(colSums(data)>1500)/384*100),col="red",cex=cex)
    draw.color.scale(xleft=1, xright=24,ybottom=scale.bot, ytop=scale.top, at=ticks, sep=0.2,cex.axis=0.5,
                     col=counts.palette, main="%")
  }

  .plot.grid(main="gene txpts",emptywells=emptywells)
  l <- log10(gene.total); l[!is.finite(l)] <- NA
  col <- colorize(x=l, counts.palette, na.col='white')
  ticks <- 0:4 
  points(coordinates,pch=19,col=col, cex=cex.wells)
  mtext(sprintf(">1000 unique genes: %.0f%%",sum(colSums(genes)>1000)/384*100),col="red",cex=cex)
  draw.color.scale(xleft=1, xright=24,ybottom=scale.bot, ytop=scale.top, at=ticks, sep=0.2,cex.axis=0.5,
                   col=counts.palette, main="log10")

  .plot.grid(main="ERCC txpts",emptywells=emptywells)
  l <- log10(spike.total); l[!is.finite(l)] <- NA
  col <- colorize(x=l, counts.palette, na.col='white')
  ticks <- -1:4
  points(coordinates,pch=19,col=col, cex=cex.wells)
  mtext(sprintf(">100 ERCCs : %.0f%%",sum(spike.total>100)/384*100),col="red",cex=cex)
  draw.color.scale(xleft=1, xright=24, ybottom=scale.bot, ytop=scale.top, at=ticks, sep=0.2, cex.axis=0.5,
                   col=counts.palette, main="log10")
  
  .plot.grid(main="ratio ERCCs/genes", emptywells=emptywells)

  l <- log2(spike.total/gene.total); l[!is.finite(l)] <- NA
  col <- colorize(x=l, counts.palette, na.col='white')
  points(coordinates,pch=19,col=col, cex=cex.wells)
  mtext(sprintf(">ERCC/gene > 0.05: %.0f%% (non-empty only)",sum(na.rm=TRUE, (spike.total/gene.total)>0.05)/384*100),col="red",cex=cex)
  draw.color.scale(xleft=1, xright=24,ybottom=scale.bot, ytop=scale.top, at=ticks, sep=0.2, cex.axis=0.5,
                   col=counts.palette, main="log2")
  par(mar=save.mar, xpd=FALSE)
}                                       # plate.plots

# plot the top 20 genes with expression bar and then barplot their index of dispersion
topgenes<-function(data){
    data<-rmspike(data)
    means<-apply(data,1,mean)
    vars<-apply(data,1,var)
    idx.disp<-vars/means
    means<-means[order(means, decreasing = TRUE)]
    idx.disp<-idx.disp[order(idx.disp, decreasing = TRUE)]
    names(means)<-sapply(names(means),chop_chr)
    names(idx.disp)<-sapply(names(idx.disp),chop_chr)
    barplot(log2(rev(means[1:20])),las=1,cex.names = 0.5, main="top expressed genes", xlab="log2(mean expression)",horiz=TRUE)
    barplot(log2(rev(idx.disp[1:20])),las=1,cex.names = 0.5, main="top varying genes",xlab="log2(var/mean)",horiz=TRUE)
}


#Read files in specified directory automatically (based on Thoms script)
read_files <- function(dir = ""){
  
  #add "/" to dir
  if(substr(dir, start = nchar(dir), stop = nchar(dir)) != "/" && dir != ""){
    dir <- paste(dir, "/", sep = "")
  }
  
  #Read files
  files <- list.files(dir, ".cout(t|b|c).csv")
  split <- strsplit(files,split = ".cout")
  file_names <- unique(as.character(data.frame(split, stringsAsFactors = FALSE)[1,]))
  
  #This check if all necessary files are in the script
  error <- ""
  for(i in 1:length(file_names)){
    
    if(file.exists(paste(dir, file_names[i],".coutb.csv", sep="")) == FALSE){
      f <- paste(file_names[i], ".coutb.csv", " is not found!", sep = "")
      error <- paste(error, "\n", f)
    }
    
    if(file.exists(paste(dir, file_names[i],".coutc.csv", sep="")) == FALSE){
      f <- paste(file_names[i], ".coutc.csv", " is not found!", sep = "")
      error <- paste(error, "\n", f)
    }
    
    if(file.exists(paste(dir,file_names[i],".coutt.csv", sep="")) == FALSE){
      f <- paste(file_names[i], ".coutt.csv", " is not found!", sep = "")
      error <- paste(error, "\n", f)
    }
  }
  
  if(error != ""){
    stop(error)
  }
  cat("the following plates will be processed:\n")
  print(file_names)

  output <- paste(dir,file_names, sep="")
  return(output)
}

wellname <- function(i=NULL,j=NULL) {
    ## returns name of all wells (i==NULL), or given its index (J=NULL) or coordinates
    rows <- LETTERS[1:16]
    cols <- 1:24
    m <- matrix(kronecker(X=rows, Y=as.character(cols), FUN=paste0),byrow=TRUE,
                nrow=length(rows),ncol=length(cols), dimnames=list(rows,cols))
    if(is.null(i))
      return(t(m))
    if(is.null(j))
      t(m)[i]
    else
      m[row,col]
}


## check expression in empty wells of plate and calculate "leakyness" from highly expressed genes
leakygenes<-function(data, emptywells) {
  empty<-data[emptywells] # subset data to the wells specified as empty
  names(empty)<-sapply(emptywells, wellname)
  genes<-apply(data,2,function(x) sum(x>=1)) # check how many genes are detected
  genes.empty<-apply(rmspike(empty),2,function(x) sum(x>=1)) # remove ERCC reads
  spike.empty<-colSums(keepspike(empty))
  empties<-data.frame(genes=genes.empty,ERCC=spike.empty)
  n <- sum(genes.empty > median(genes)/5)
  if(n>0)
    warning(sprintf("plate %s: %d/%d empty wells contain >= median/5 reads", names[[i]], n, length(emptywells)))
  ## plot genes/cell and ERCC reads/cell for empty wells
  par(mar = c(5, 4, 6, 1))
  barplot(t(empties),main="transcripts in empty wells",
          col=c("blue","red"),space=rep(c(5/nrow(empties),0),nrow(empties)),cex.names = 0.5,las=3,beside=TRUE,
          legend=colnames(empties),
          args.legend = list(x = "topright", bty = "n",horiz=TRUE,inset=c(0,-0.25)))
  
  # determine top expressed genes in empty and compare to mean expressed genes in plate, but
  # only use empty wells that have at least 75 ERCC reads (otherwise that well did not work)
  use <- spike.empty > 75
  if( sum(use)== 0 ) {
    .empty.plot(main="top 10 of overlap between \n top50-empty and top200-all", msg="no empty wells with >= 75 reads")
    .empty.plot(main="genes from top50-empty\n not in top200-all", msg="no empty wells with >= 75 reads")
    return()
  }
  if( sum(use)== 1 )
    warning(paste("Just one empty well with at least 75 ERCC reads, in", names[[i]]))
  emptyz<-empty[use]  # take only wells that worked
  emptyz<-rmspike(emptyz)
  top.empty<-apply(emptyz,1,sum)[order(apply(emptyz,1,sum),decreasing=TRUE)][1:50] # pick top 50 in empty
  top.all<-apply(data,1,sum)[order(apply(data,1,sum),decreasing=TRUE)][1:200] # pick top 200 in plate
  names(top.empty)<-sapply(names(top.empty),chop_chr) # remove __chr* from name
  names(top.all)<-sapply(names(top.all),chop_chr) # remove __chr* from name
  o <- names(top.empty) %in% names(top.all)
  if(sum(o)==0) {
    .empty.plot(main="top 10 of overlap between \n top50-empty and top200-all", msg="no overlap!")
    .empty.plot(main="genes from top50-empty\n not in top200-all", msg="all genes!")
    return()
  } 
  overlap<-top.empty[o] # check overlap between top 50 empty and 200 in plate
  non.overlap<-top.empty[ !o ]
  b<-barplot(log2(rev(overlap[1:10])),las=1,cex.names = 0.6, main="top 10 of overlap between \n top50-empty and top200-all",
             sub="leakage to empty in %", xlab="log2(sum(reads in empty))",horiz=TRUE)
  text(0.5,b, round((top.empty[names(overlap)[1:10] ]/top.all[names(overlap)[1:10] ])*100,2))
  if (length(overlap)==50)
    warning(paste("there is complete overlap between empty genes and plate genes in ", names[[i]]))
  barplot(log2(rev(non.overlap)),las=1,cex.names = 0.6,
          main="genes from top50-empty\n not in top200-all", xlab="log2(mean expression)",horiz=TRUE)
}                                       #leakygenes

colorize <- function(x, palette, min=NULL, max=NULL, na.col='grey') {
## taken from svn/pub/tools/general/R/uuutils/R/uuutils.R rev 1240
  if(length(dim(x))>0)
    stop("colorize: need simple vector")
  stopifnot(is.numeric(x))
  names <- names(x)

  if(is.null(min)) min <- min(x, na.rm=TRUE)
  if(is.null(max)) max <- max(x, na.rm=TRUE)
  
  x[ x<min ] <- min
  x[ x>max ] <- max

  n.colors <- length(palette)
  breaks <- seq(from=min, to=max, length.out=n.colors+1)
  levels <- cut(x, breaks=breaks,right=FALSE, labels=FALSE, include=TRUE)
  col <- palette[levels]
  names(col) <- names
  col[is.na(col)] <- na.col
  col
}                                       #colorize

draw.color.scale <- function(xleft,xright, ybottom, ytop, # e.g. from locator
## taken from svn/pub/tools/general/R/uuutils/R/uuutils.R rev 1237
                        range=NULL, col,
                        main=NA, # main: the title (should call this label?)
                        cex.main=0.8, 
                        offset.main=0.5,
                        at=c(range[1],mean(range),range[2]),# tick label locations in scale coordinates
                        labels=as.character(at), # the tick labels
                        cex.axis=cex.main,
                        offset.axis=NA,
                        sep=0.5
                         ) {
  ## needs corners of the strip with the colors (get an initial
  ## position using locator()), rest is done automagically.
  ## See also plotrix::color.legend
  if(is.null(at) && is.null(range))
    stop("Not both the 'at' and the 'range' arguments can be NULL")
  if(is.null(at))
    at <- c(range[1],mean(range),range[2])
  if(is.null(range))
    range <- range(at)

  vertical <- (  abs(xleft-xright) < abs(ybottom-ytop)  )
  if(is.na(offset.axis))
    offset.axis <- if(vertical)0.2 else 0.4
  
  x <- range(c(xleft,xright))
  y <- range(c(ybottom,ytop))
  width <- x[2]-x[1]
  height <- y[2] - y[1]
  n.colors <- length(col)
  if(vertical) { 
    ys <- seq(from= y[1], to= y[2], length.out=n.colors+1)
    rect(xleft= x[1],
         xright = x[2],
         ybottom= ys[1:n.colors],
         ytop=ys[2:(n.colors+1)],
         col=col,
         border=NA)
    if(is.character(main)|| is.expression(main)) ## Stupid warning on is.na(some.expression)
      text(x=mean(c(xleft,xright)), y=ys[n.colors+1]+sep*width, offset=offset.main, pos=3,
           labels=main, cex=cex.main)

    ## axis:
    ys <- y[1]+height*(at-range[1])/(range[2]-range[1])
    segments(x0=xright+sep*width, x1=xright+sep*width, y0=ys[1],y1=ys[length(ys)]) #line
    segments(x0=xright+sep*width, x1=xright+(2*sep)*width, y0=ys,y1=ys) #ticks
    text(x=xright + 2*sep*width, y=ys, offset=offset.axis, pos=4,
         labels=labels, cex=cex.axis)
    
  } else {
    xs <- seq(from= x[1], to= x[2], length.out=n.colors+1)
    rect(ybottom= y[1],
         ytop = y[2],
         xleft= xs[1:n.colors],
         xright=xs[2:(n.colors+1)],
         col=col,
         border=NA
         )

    if(is.character(main)||is.expression(main)) ## title to the right of the bar
      text(y=ybottom + 0.5*height, x=xs[n.colors+1], pos=4,
           labels=main, cex=cex.main, offset=offset.main)
    xs <- x[1]+width*(at-range[1])/(range[2]-range[1])
    segments(y0=ybottom-sep*height, y1=ybottom-sep*height, x0=xs[1],x1=xs[length(xs)]) #line
    segments(y0=ybottom-sep*height, y1=ybottom-(2*sep)*height, x0=xs,x1=xs) #ticks
    text(y=ybottom-2*sep*height, x=xs, labels=labels, cex=cex.axis,
         offset=offset.axis, pos=1)
  }
}                                     #draw.color.scale

read.counts <- function(file) {
  wellnames <- as.vector(wellname())
  x <- read.csv(file=file, header = TRUE, sep = "\t",row.names =1, comment.char="#")
  if (! setequal(colnames(x), wellnames))
    stop("Wellnames should be A1-24, B1-24, ..., P1-P24")
  x <- x[,wellnames]
  x
}

read.stats <- function(file) {
  wellnames <- as.vector(wellname())
  x <- read.csv(file=file, nrows=4, header = TRUE, sep = "\t",row.names =1, comment.char="")
  rownames(x) <- gsub("^#","", rownames(x))
  if (! setequal(colnames(x), wellnames))
    stop("Wellnames should be A1-24, B1-24, ..., P1-P24")
  x <- x[,wellnames]
  x
}                                       #read.stats


# Local variables:
# mode: R
# ess-indent-level: 2
# End:
