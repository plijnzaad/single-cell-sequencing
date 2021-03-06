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
  mtext(sprintf("median: %.1f",median(rc.v/bc.v)), col="DarkGreen", cex=0.6)
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

.rightmost.nchar <- function(string,n)substring(string, first=(nchar(string)-n+1))

infobox <- function(script, dir, filename,totals) {
  cex <- 0.6
  ## print version filename and some overal statistics
  plot(type="n", x=c(0,1), y=c(0,1), axes=FALSE, xlab=NA,ylab=NA, main="Information", cex=cex)
  dir <- dirname(normalizePath(dir))
  len <- nchar(dir); maxlen <- 30
  if(len>maxlen)
    dir <- substr(dir,start=(len-maxlen)+1, stop=len)

  unmapped <- totals['unmapped']
  valid <- totals['reads']
  totreads <- valid + unmapped #note: ignores all those without a valid CBC, and more!
  text <- c(
    paste0("script version: ", getversion(script)),
    paste0("dir: ", .rightmost.nchar(dir,55)),
    paste0("file: ", .rightmost.nchar(filename,55)),
    sprintf("refgenes found: %d", totals['ngenes']),
    sprintf("ERCCs found: %d", totals['nspikes']),
    sprintf("total reads with CBC: %s\n\t mapped: %s (%.1f %%) ",
            commafy(totreads), commafy(valid), 100*valid/totreads),
    paste0("total umis (mapped, valid CBC): ",  commafy(totals['umis'])),
    paste0("total txpts(mapped, valid CBC): ",  commafy(totals['txpts']))
    )

  text(pos=4, x=0, y=seq(1, 0, length.out=length(text)), labels=text, cex=cex)
}                                       #infobox


failedwells.gen <- function(spikes, minmed=2, maxmissing=0.75) {
    use <- apply(spikes,1,median)>=minmed
    if (maxmissing %% 1 > 0)
      maxmissing <- sum(use)*maxmissing   # if more than this many 'good' ERCCs miss, reaction has failed 
    absent <- spikes[use,] == 0 
    absent <- apply(absent, 2, function(col)sum(col)>= maxmissing)
    names(absent[absent])
}

failedwells <- function(spikes) {
  use <- apply(spikes,1,median)>=4
  absent <- spikes[use,]==0
  absent <- apply(absent, 2, function(col)sum(col)>= 2)
  return(names(absent[absent]))
  ## additional criterion (too many!)
  surv <- spikes[use, !absent]
  mn <- apply(log2(surv), 2, mean)
  med <- apply(log2(surv), 2, median)
  ma <- apply(log2(surv), 2, mad)
  m <- abs(mn - med)> ma
  union(names(absent[absent]), names(m[m]))
}                                       #failedwells

#plot total number of transcripts per sample
totaltxpts <- function(txpts,plotmethod=c("barplot","hist","cumulative","combo"),
                       emptywells=NULL) {
  col.full <- rgb(0, 0.6, 0, 0.8)
  col.empty <- rgb(0, 0, 0.6, 0.6)
  magnif.empty <- 10
  cex <- 0.6
  ## if ( ! plotmethod %in% c("barplot","hist","cumulative","combo") ) stop("invalid method")
  if(plotmethod == "hist"){
    logdata <- log10(txpts)
    logdata[ logdata == -Inf ] <- log10(0.5)+0.01 # just to the right of to the zero-tick
    a<-hist(plot=FALSE,x=logdata,breaks=100, right=FALSE)
    ticks <- log10(as.vector(c(1,2,5) %o% 10^(-1:9)))
    first <- which(ticks < min(a$breaks)); first <- first[length(first)]
    last <- which(ticks > max(a$breaks))[1]
    ticks <- ticks[first:last]
    rest.args <-list(xlab="counts",ylab="frequency",main="gene transcripts/well",
                     xaxt="n",col.sub="red", breaks=a$breaks,
                     xlim=range(ticks), ylim=c(0, 1.2*max(a$counts))) 
    if(is.null(emptywells)) { 
      do.call(hist, args=c(list(x=logdata,
                      col=col.full, border=NA), rest.args))
    } else {
      do.call(hist, args=c(list(x=logdata[-emptywells],
                      col=col.full, border=NA), rest.args))
      do.call(hist, args=c(list(x=rep(logdata[emptywells],magnif.empty),
                      col=col.empty, border=NA, add=TRUE), rest.args))

      legend(x = "topleft", legend = c("full", sprintf("empty x %d", magnif.empty)),
             fill = c(col.full, col.empty), bty="n", cex=cex)
    }
    axis(1, at=ticks,labels=sprintf("%s",commafy(round(10^ticks))),las=3, cex.axis=cex)
    mn <- mean(txpts)
    md <- median(txpts)
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
    b<-barplot(txpts,xaxt="n",xlab="cells",sub=paste("mean total read:",round(mean(txpts))),main="total unique reads",col="black",border=NA) 
    axis(1,at=b,labels=c(1:length(data))) # 1=horizontal at = position of marks
    abline(h=mean(txpts),col="red")
  }
  if(plotmethod == "cumulative"){
    plot(ecdf(txpts),xlab="total reads",ylab="fraction",main="total unique reads",col="red",tck=1,pch=19,cex=cex,cex.axis=0.8) 
    abline(v=mean(txpts/2),col="red")
    mtext(paste("mean:",round(mean(txpts))," median:",round(median(txpts))),side=3,col="DarkGreen",cex=cex)
  }
  
  if(plotmethod == "combo"){
    a<-hist(log10(txpts),breaks=100,xlab="log10(counts)",ylab="frequency",main="total unique reads",col="grey",xaxt="n",col.sub="red") 
    mtext(paste("mean:",round(mean(txpts))," median:",round(median(txpts))),side=3,col="DarkGreen",cex=cex)
    axis(1,at=a$breaks[which(a$breaks %in% c(0,1,2,3,4,5))],labels=a$breaks[which(a$breaks %in% c(0,1,2,3,4,5))])
    abline(v=log10(mean(txpts)/2),col="red")
    text(log10(mean(txpts)/2),max(a$counts)-2, round(mean(txpts)/2), srt=0.2, col = "red",pos=2)
    plotInset(log10(1),max(a$counts)/4,log10(250), max(a$counts),mar=c(1,1,1,1),
              plot(ecdf(txpts),pch=".",col="red",cex=cex,ylab=NA,xlab=NA,main=NA,cex.axis=0.8,xaxt="n",las=3,mgp=c(2,0.1,0),tck=1,bty="n"),
              debug = getOption("oceDebug"))
  }
}                                       # totaltxpts

#plot complexity (number of unique genes detected) per cell
cellgenes<-function(complexity,plotmethod=c("hist","cumulative","combo")) {
  cex <- 0.6
  if ( ! plotmethod %in% c("hist","cumulative","combo") ) stop("invalid plotting method")

  if(plotmethod == "hist"){
    a<-hist(complexity,breaks=100,xlab="total genes",ylab="frequency",main="unique genes detected/well",col="steelblue1",xaxt="n") 
    mtext(paste("mean:",round(mean(complexity))," median:",round(median(complexity))),side=3,col="DarkGreen",cex=cex)
    axis(1,at=a$breaks[which(a$breaks %in% seq(0,max(a$breaks),1000))],labels=a$breaks[which(a$breaks %in% seq(0,max(a$breaks),1000))])
  }

  if(plotmethod == "cumulative"){
    plot(ecdf(complexity),pch=19,cex=0.2,col="red",ylab="frequency",xlab="detected genes/well",main="cumul. complexity",cex.axis=1,las=1,tck=1)
    mtext(paste("mean:",round(mean(complexity))," median:",round(median(complexity))),side=3,col="DarkGreen",cex=cex)
  }

  if(plotmethod == "combo"){
    ## a<-hist(complexity,breaks=100,ylab="frequency",main="complexity",xlab="unique genes/well",col="steelblue1",xaxt="n")
    d <- density(complexity, adjust=0.5)
    xlim <- c(0, max(complexity))
    plot(d,ylab="density",main="complexity",xlab="unique genes/well",xlim=xlim, col="blue", yaxt="n")
    rug(complexity, ticksize= 0.02)
    mtext(paste("mean:",round(mean(complexity))," median:",round(median(complexity))),side=3,col="DarkGreen",cex=cex)
    ### axis(1,at=a$breaks[which(a$breaks %in% seq(0,max(a$breaks),1000))],labels=a$breaks[which(a$breaks %in% seq(0,max(a$breaks),1000))])

    plotInset(xleft=xlim[2]/2.5, xright=xlim[2], ybottom=max(d$y)/2, ytop=max(d$y),mar=c(1,1,1,1),
              plot(ecdf(complexity),xlim=c(0, 1.1*max(complexity)), xlab=NA,
                   pch=NA,lwd=1,col="red",cex=cex,ylab=NA, main=NA, cex.axis=0.6, las=3),
              debug = getOption("oceDebug"))
  }
}                                       #cellgenes

well.coverage <- function(main, gene.total) {
  cex <- 0.6
  max <- ceiling(log10(max(gene.total)))
  ticks <- as.vector(c(1,2,5) %o% 10^(0:max))
  ticks <- ticks[ticks<max(gene.total)]
  freq <- sapply(ticks, function(cnt)sum(gene.total >= cnt))
  xlim <- log10(range(ticks))
  ylim <- c(0, 400)
  plot(cex=0.5,main=main, x=log10(ticks), y=freq, 
       type="n",tck=1, 
       lwd=2 ,axes=FALSE, xlab="N", ylab="wells with >= N transcripts",
       xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")
  axis(1, at=log10(ticks), labels=sprintf("%s",commafy(ticks)),las=3, 
       cex.axis=cex,tck=0, pos=min(ylim))
  axis(2, at=seq(0,400,by=50), labels=TRUE, 
       cex.axis=cex, tck=0, las=1, pos=min(xlim))
  abline(h=seq(0,400,by=50), v=log10(ticks), col="grey")
  axis(4, at=seq(0,384,length.out=11), labels=paste0(seq(0, 100,length.out=11), '%'),
       cex.axis=cex, tck= -0.01, las=1, pos=max(xlim))
  lines(pch=19, cex=0.5,x=log10(ticks), y=freq, 
       type="o",tck=1, col="red",
       lwd=2, xaxs="i", yaxs="i")
  mtext(sprintf("wells without txpts: %d", sum(gene.total==0) ), side=3, col="DarkGreen", cex=cex)
}                                       #well.coverage

## Support for use of Michaelis-Menten for saturation curves comes from e.g. 
## Pascal de Caprariis, Richard Lindemann, and Robert Haimes
## Mathematical Geology, Vol. 13, No. 4, 1981
## "A Relationship Between Sample Size and Accuracy of Species Richness Predictions"
nop <- function(...)invisible(NULL)

R.squared <- function(model, y) {
    RSS <- sum(residuals(model)^2)
    TSS <- sum((y - mean(y))^2)
    (1 - (RSS/TSS))
}

reads.required <- function(K, f)K*f/(1-f) 

format.kMG <- function(x, fmt="%.3e") {
  ## warning("format.kMG: fix for values < 1!")
  return(sprintf(fmt,x))
  ##   stopifnot(all( 1<x  & x < 1e+27))
  prefixes <- c(" ", "k","M", "G", "T", "P", "E", "Z", "Y")
  f <- floor(floor(log10(x))/3)
  prefix <- prefixes[f+1]
  om <- 10^(3*floor(floor(log10(x))/3))
  scaled <- x/om
  sprintf(paste0(fmt," %s"), scaled, prefix) 
}                                       #format.kMG

fit.details <- function(model,y, pred="") {
  if(is.null(model))
    return("(no fit")
  R2 <- R.squared(model,y)
  coef <- coef(model)                   #has names Vm and K

  extra <- ""
  if( names(coef)[1]== "Vm") { # Michaelis-Menten
    f <- 0.9
    nearing.at <- reads.required(K=coef['K'], f=f)
    ## extra <- sprintf("%.1f%% at %s reads", floor(f*100), format.kMG(nearing.at))
  }

  ## coef.string <- paste(collapse="  ", sep=": ", names(coef), format.kMG(coef))
  ## sprintf("%s: R2: %.3f\n%s\n%s", pred, R2, coef.string, extra)
  sprintf("plateau at %.1e (R2 = %.3f)\nreaching 90%% at %.1e reads", coef['Vm'], R2, nearing.at)
} # fit.details

mm.fit <- function(data) {
  ss <- NULL

  tryCatch({ ss <- getInitial(data=data,y ~ SSmicmen(input=reads, Vm=max(y), K=1))},
           error=nop,warning=nop,finally=nop)

  if(is.null(ss))
    return(NULL)
  
  Vm <- ss[1]
  K <- ss[2]

  model <- NULL
  tryCatch({ model <- nls(data=data, y~SSmicmen(reads, Vm, K))},
           error=nop,warning=nop,finally=nop)
  return(model)
}                                       #mm.fit

save.par <- function(par, names) {
  keep <- list()
  for(name in names)
    keep[[name]] <- par[[name]]
  keep
}

restore.par <- function(saved) {
  for(name in names(saved))
    do.call(par, args=saved[name])
}

saturation.plot <- function(main, data, xcol, ycol,
                            diffdata, ycoldiff, xlab, ylab, pred="", max=NULL, ...) {
  ## max=NULL: no relative scale on left; max=NA: use max(data[,ycol]) as max

  if (is.null(data) || is.null(data[,xcol]) || is.null(data[,ycol]) ) { 
    .empty.plot(main=main, msg="no data")
    return()
  }

  cex <- 0.6
  
  x <- data[,xcol]
  y <- data[,ycol]
  rug <- data[ data$reads %% 1e6 <= 1, 'nmapped']
  
  par <- par(no.readonly=TRUE)
  keep.par <- save.par(par,c("mar", "xpd"))
  on.exit(restore.par(keep.par))
  mar <- par$mar
  mar[2] <- mar[2]*1.2
  mar[4] <- mar[4]*1.2
  par(mar=mar, xpd=FALSE)
  
  plot(main=main, x=c(0, max(x)), y=c(0, max(y)), xlab=xlab,ylab=ylab, ...,
       type="n", lwd=2, col="black", cex.axis=cex, las=2, tck= -0.01, xaxs="i", yaxs="i")

  atx <- axTicks(1)
  aty <- axTicks(2)
  axis(1, at=atx, labels=FALSE, cex.axis=cex)
  axis(2, at=aty, labels=FALSE, cex.axis=cex)
  abline(h=aty, v=atx, col="grey")

  maxy <- max(y)
  if(! is.null(max)) {                  #relative scale on left
    if(is.na(max))
      max <- maxy
    max <- unname(max)               #confusing
    maxperc <- maxy/max*100          # exceeds 100 when using eg. 95%-ile as max
    at.perc <- axisTicks(c(0,maxperc), FALSE, nint=5)
    axis(2, at=at.perc/100*max, labels=sprintf("%10d%%",at.perc), hadj=0, # leftpad with spaces to set them on rhs of axis
         cex.axis=cex, tck= 0.01, las=1)
  }
  
  if(!is.null(rug))                     
    rug(x=rug, ticksize= 0.01, col="purple")
  
  data <- data.frame(y=y, reads=x)      #actually mapped reads

  model <- mm.fit(data=data)

  if(!is.null(model)) {
    fitted <- fitted(model)
    max <- coef(model)['Vm']
    lines(x=x, y=y, ..., lwd=2, col="black")
    abline(h=max,col="black", lty=3)
    lines(data=data, fitted~reads, col="black", lwd=2, lty=2)
    details <- fit.details(model, data$y, pred='full data')
    text(x=0.5*max(x), y=0.5*max(y), adj=c(0.5, NA), labels=details, cex=cex, col="black")
  }

  ### ---- derivative:
  col.deriv="blue"
  diffx <- diffdata[,xcol]
  d <- diffdata[,ycoldiff]
  if( any(d < 0) ) {
    warning("saturation curve is not monotonic (prolly umis/gene?)")
    legend(x="bottomright", legend="data is not not monotonic", text.col="red", bty="n", cex=1)
  }

  ndiv <- 10
  f <- maxy/max(d)                    # converting to same scale as left axis
  lines(type="o", pch=19, cex=0.3, y=f*d, x=diffx, lwd=2, col=col.deriv)

  ### scale for derivative
  at.diff <- axisTicks(c(0,max(d)), FALSE, nint=ndiv)
  axis(2, at=at.diff*f, labels=at.diff, pos=max(x), cex.axis=cex,
       las=1, col=col.deriv, col.ticks=col.deriv, col.axis=col.deriv)

  ### also relative to max of main curve
  maxperc <- max(d)/maxy*100
  at.perc <- axisTicks(c(0,maxperc), FALSE, nint=ndiv)
  axis(4, at=at.perc/max(maxperc)*maxy,
       labels=sprintf("%d%%",at.perc), hadj=0.5, # leftpad with spaces to set them on rhs of axis
       cex.axis=cex, tck= -0.01, las=1, col=col.deriv, col.ticks=col.deriv, col.axis=col.deriv)

  mtext("extra items per 1M mapped reads", side=4, line=1.4, col=col.deriv, cex=cex*0.8)
  return()
}                                       # saturation.plot()

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
}                                       #testcutoff

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
}                                       #.plot.grid

#plot number of total reads, ERCC-reads and genes/well over a 384-well plate layout
plate.plot<-function(data, main, ticks, scale.name, emptywells=NULL, hilite=NULL,
                     mtext=NULL, failedwells=NULL) {
  cex <- 0.6
  cex.wells <- 1.1
  if( grepl("X11",names(dev.cur())) )
    cex.wells <- 3.5
  pch.infinite <- 4 # small 'x' for Inf and NaN
  col.failed <- "black"

  if(!is.null(hilite)) {
    stopifnot(is.logical(hilite)&&length(hilite)==384)
    if(is.null(mtext))
      stop("Need mtext to explain/quantify the hilites")
    tot <- sum(hilite)
    mtext <- sprintf("%s: %d wells (%.0f%%)", mtext, tot, tot/384*100)
  }

  low <- "#313695";mid="#FFFFBF"; high="#A50026"   ## three-class RdYlBu, 11 levels
  RdYlBu.orig <- colorRampPalette(c(low,mid,high))(91)
  col.hilite <- "DarkGreen"

  low <- '#91BFDB';  mid <- '#FFFFBF'; high <- '#FC8D59';  ## three-class RdYlBu, pastel
  RdYlBu.pastel <- colorRampPalette(c(low,mid,high))(91)
  col.hilite <- "magenta"

  black <- "#000000"; red <- "#FF0000"; yellow <- "#FFFF00"
  redhotpoker <- colorRampPalette(c(black,red,yellow))(91)
  col.hilite <- "DarkGreen"

###  counts.palette <- RdYlBu.pastel
  counts.palette <- RdYlBu.orig
###  counts.palette <- redhotpoker

  par <- par(no.readonly=TRUE)
  keep.par <- save.par(par,c("mar", "xpd"))
  on.exit(restore.par(keep.par))
  mar <- par$mar
  mar[1] <- mar[1]*1.3                  #to get right aspect ratio for the plates
  mar[3] <- mar[3]*1.3
  par(mar=mar, xpd=TRUE)
  scale.top <- -2.5
  scale.bot <- -3.5

  if(is.null(data)) { 
    .empty.plot(main=main, msg="no data")
    return()
  }
  coordinates <- expand.grid(seq(1,24),rev(seq(1,16)))
  rownames(coordinates) <- as.vector(wellname())
  .plot.grid(main=main,emptywells=emptywells)
  infinite <- is.infinite(data) | is.nan(data)
  data[infinite] <- NA
  col <- colorize(x=data, counts.palette, na.col='white')
  points(coordinates,pch=19,col=col, cex=cex.wells)
  points(coordinates[infinite,],pch=pch.infinite, cex=0.5*cex.wells, col='black')

  if(!is.null(hilite))
    points(coordinates[hilite,],pch=1, cex=cex.wells*1.1,col=col.hilite, lwd=0.3)

  if(!is.null(failedwells))
    points(coordinates[failedwells,],pch=4, cex=cex.wells*1.1,col=col.failed, lwd=0.3)

  draw.color.scale(xleft=1, xright=24,ybottom=scale.bot, ytop=scale.top, at=ticks, sep=0.2, cex.axis=0.5,
                   col=counts.palette, main=scale.name)

  if(!is.null(mtext))
    mtext(text=mtext,col=col.hilite,cex=cex)
  return()
}                                       # plate.plot

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
get.file.names <- function(dir = "./", regexp=".*") {
  pat <- paste0(regexp,"\\.cout[cbt]\\.csv")
  files <- list.files(path=dir, pattern=pat)
  split <- strsplit(files,split = ".cout")
  names <- unique(unlist(lapply(split,function(l)l[1])))
  for (name in names) { 
    for (letter in c("c", "b", "t")) {
      expect <- sprintf("%s/%s.cout%s.csv", dir, name, letter)
      if(! file.exists(expect))
        stop(expect, " not found")
    }
  }
  return(names)                         #note: without the directory
}                                       #get.file.names

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
}                                       #well.name

## check expression in empty wells of plate and calculate "leakyness" from highly expressed genes
leakygenes<-function(plate, data, emptywells) {
  empty <- data[,emptywells] # subset data to the wells specified as empty
  names(empty)<-sapply(emptywells, wellname)
  genes<-apply(data,2,function(x) sum(x>=1)) # check how many genes are detected
  genes.empty<-apply(rmspike(empty),2,function(x) sum(x>=1)) # remove ERCC reads
  spike.empty<-colSums(keepspike(empty))
  empties.frame<-data.frame(genes=genes.empty,ERCC=spike.empty)
  n <- sum(genes.empty > median(genes)/5)
  if(n>0)
    warning(sprintf("plate %s: %d/%d empty wells contain >= median/5 reads", plate, n, length(emptywells)))

  par <- par(no.readonly=TRUE)
  keep.par <- save.par(par,c("mar"))
  on.exit(restore.par(keep.par))
  par(mar = c(5, 4, 6, 1))

  ## plot genes/well and ERCC reads/well for empty wells
  barplot(t(empties.frame),main="transcripts in empty wells",
          col=c("blue","red"),space=rep(c(5/nrow(empties.frame),0),nrow(empties.frame)),cex.names = 0.5,las=3,beside=TRUE,
          legend=colnames(empties.frame),
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
    warning(paste("Just one empty well with at least 75 ERCC reads, in", plate))
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
  b<-barplot(log2(rev(overlap[1:10])),las=1,cex.names = 0.6,
             main="top 10 of overlap between \n top50-empty and top200-all",
             sub="leakage to empty in %",
             xlab="log2(sum(txpts in empty))",horiz=TRUE)
  text(0.5,b, round((top.empty[names(overlap)[1:10] ]/top.all[names(overlap)[1:10] ])*100,2))
  if (length(overlap)==50)
    warning(paste("there is complete overlap between empty genes and plate genes in ", plate))
  barplot(log2(rev(non.overlap)),las=1,cex.names = 0.6,
          main="genes from top50-empty\n not in top200-all", xlab="log2(mean expression)",horiz=TRUE)
  return()
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

read.wellwise.saturation <- function(file) {
  if(! file.exists(file))
    stop("No wellwise saturation file found (expected ", file, ")")

  tab <- read.table(file=file, sep="\t",
                    as.is=TRUE, quote="", header=TRUE,comment.char="", row.names=NULL)
  names(tab)[1] <- gsub("^(X\\.)?n?", "",names(tab)[1])
  if( !setequal(colnames(tab), c("reads", as.vector(wellname()))))
    stop("wrong column names; expected 'nreads \\t A1 \\t A2 ... \\t P24")
  tab
}                                       #read.wellwise.saturation

read.saturations <- function(name, diff.every=100) {
  file <- paste0(name, "-saturation.txt")
  
  if(! file.exists(file))
    stop("No saturation files found (expected ", file, ")")
    
  all <- read.table(file=file, sep="\t",
                    as.is=TRUE, quote="", header=TRUE,comment.char="", row.names=NULL)

  if(all(all[1,]==1)) all <- all[-1,]
  
  names(all)[1] <- gsub("^X\\.n?", "",names(all)[1])
  expected.cols <- c("reads", "nmapped", "nvalid", "genes", "umis", "txpts")
  missing <- setdiff(expected.cols, colnames(all))
  if (length(missing)>0)
    stop("Columns missing from ", file, ": ", paste(missing,collapse=" "))

  sats <- list()

  sats$all <- all               #lumps everything together, discarding well info
  sats$perwell <- list()        #information for each well individually
  sats$wellwise_mean <- list()  #averages taken of over all wells

  sats$perwell_diff <- list()     #information for each well individually
  sats$wellwise_mean_diff <- list()  # diff of wellwise mean (identical to mean of wellwise diffs!!)

  ## wellwise files:
  for (type in c('genes', 'umis')) { 
    file <- sprintf("%s-wellsat-%s.txt", name, type)
    tab <- read.wellwise.saturation(file)
    if(is.null(tab))
      return(NULL)
    m <- merge(all, tab,by='reads')
    sats$perwell[[type]] <- m
    ## mean:
    aggr <- cbind(tab[,1], apply(tab[,-1], 1, mean))
    colnames(aggr) <- c("reads", "mean")
    sats$wellwise_mean[[ type ]] <- merge(all, aggr, by='reads')
  }


  ## differences between recording points:
  n <- nrow(all)
  x <- all[,'nmapped']
  every <- c(seq(1, n, by=diff.every))
  x.every <- x[ every[-1] ]           #skip first one
  all.diff <- all[every[-1],]

  for (type in c('genes', 'umis')) { 
    m <- sats$perwell[[type]]

    ## diff of overal wellwise mean (is same as mean of wellwise differences)
    wellwise_mean <- sats$wellwise_mean[[type]][, 'mean']
    d <- diff(wellwise_mean[every])
    tab <- cbind(all.diff, diffofmean=d)
    sats$wellwise_mean_diff[[type]] <- tab

    ## diffs per indiv. well:
    y <- m[names(wells)]
    y.every <- y[ every, ]

    d <- apply(y.every, 2, diff)        #this keeps the A1, A2, ... names
    tab <- cbind(all.diff, d)
    sats$perwell_diff[[type]] <- tab
  }

  ## umis per gene
  umispergene <- with(sats$wellwise_mean, umis$mean/genes$mean)
  sats$umispergene <- cbind(all, umispergene=umispergene)

  ## and their diffs
  d <- diff(umispergene[every])
  sats$umispergene_diff <- cbind(all.diff, diff=d)

  sats
}
                                        #read.saturations

# Local variables:
# mode: R
# ess-indent-level: 2
# End:
