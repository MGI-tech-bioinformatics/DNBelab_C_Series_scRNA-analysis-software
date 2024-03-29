#!/usr/bin/env Rscropt

suppressMessages({
    library(ggplot2)
    library(getopt)
    library(data.table)
    library(cowplot)
    library(DropletUtils)
})

arg<-matrix(c("input", "i","1","character","Path of input directory",
              "output","o","2","character","Path of output",
              "force","f","1","numeric","Force cells number to analysis. You must set -f or -e to process barcode list.",
              "expect","e","1","numeric","Number of expect cells",
              "umi","u","1","numeric","Number of umi count threshold",
              "low","l","1","numeric","lower bound on the total UMI count.",
              "help"  ,"h","2","logical",  "This information"),
            byrow=T,ncol=5
            )

opt = getopt(arg)

if( !is.null(opt$help) || is.null(opt$input)){
    cat(paste(getopt(arg, usage = T), "\n"))
    q()
}

if (is.null(opt$output)) {
    opt$output<-getwd()
}

bc <- fread(opt$input,header=TRUE)
bc <- as.data.frame(bc)
bc <- subset(bc, bc$UB>0)
bc <- bc[order(bc$UB, decreasing=T),]
len <- nrow(bc)
                                        #sor = sort(bc$UB, decreasing=T)
sor = bc$UB
a = log10(1:len)
b = log10(sor)
expect <- 0
cutoff <- 0
m <- 0
umi <- 0
low <-  50

if (!is.null(opt$expect)) {
    expect=as.numeric(opt$expect)
}

if (!is.null(opt$umi)) {
    umi=as.numeric(opt$umi)
}

if (!is.null(opt$low)) {
    low=as.numeric(opt$low)
}
if (!is.null(opt$force)) {
    force = as.numeric(opt$force)
}
if (umi >0) {
    cutoff<-length(which(sor>=umi))
    m <- umi
}else if (expect >0) {
    c = log10(expect*1.7)
    if(10^c>len){c=log10(expect)}
    lo <- loess(b~a,span = 0.004,degree = 2)
    print(c)
    xl <- seq(2,c, (c - 2)/10)
    out = predict(lo,xl)
    infl <- c(FALSE,abs(diff(out)/((c - 2)/10) - -1) == min(abs(diff(out)/((c - 2)/10)- -1)))
    m = 10 ^ out[infl] + 0.5
    m = round(m , digits =0 )
    cutoff<-length(which(sor>=m))
}else if (force > 0) {
    expect = force
    cutoff = expect
    m = sor[cutoff]
} else {
   # trace(barcodeRanks, quote(totals <- m[,1]), at=3)
   # test=bc[3]
   # colnames(test)="count"
   # test$count=as.numeric(test$count)
    use=bc[,c(1,3)]
    rownames(use)=use[,1]
    use=use[2]
    test=t(use)
    br <- barcodeRanks(test,lower = low)
    o <- order(br$rank)
    m <- br@metadata$inflection
    cutoff<-length(which(sor>=m))
    if(cutoff >= 10000){
        cutoff=10000
        m = bc[10000,3]
    }
}

tmp<-data.frame(barcodes=1:len,UMI=sor,cell=c(rep("true",cutoff),rep("noise",len-cutoff)))

cc <- bc[1:cutoff,]
called = sum(cc$UB)
in_cell <- round(called/sum(bc$UB),digits=4)

small<-bc[which(bc$UB>=m),]
CellNum = nrow(small)
UMI_mean = mean(small$UB)
cell_mean = mean(small$Raw)
UMI_median = median(small$UB)
Gene_mean = mean(small$GN)
Gene_median = median(small$GN)
UMI_Cell = sum(small$UB)
#Read_mean = mean(small$NUM_GENIC_READS)

cat(paste("Estimated Number of Cell,", CellNum, "\n",sep=""),file=paste(opt$output,"/../report/cell_report.csv",sep=""))
cat(paste("Reads in cell,", round(in_cell,3)*100, "%\n",sep=""),file=paste(opt$output,"/../report/cell_report.csv",sep=""), append=T)
cat(paste("Mean reads per cell,", round(cell_mean),"\n",sep=""),file=paste(opt$output,"/../report/cell_report.csv",sep=""), append=T)
cat(paste("Mean UMI counts per cell,", round(UMI_mean),"\n",sep=""),file=paste(opt$output,"/../report/cell_report.csv",sep=""), append=T)
cat(paste("Mean Genes per cell,", round(Gene_mean),"\n",sep=""),file=paste(opt$output,"/../report/cell_report.csv",sep=""), append=T)
cat(paste("Median UMI Counts per Cell,", round(UMI_median),"\n",sep=""),file=paste(opt$output,"/../report/cell_report.csv",sep=""), append=T)
cat(paste("Median Genes per Cell,", round(Gene_median),"\n",sep=""),file=paste(opt$output,"/../report/cell_report.csv",sep=""), append=T)
write.table(cc$BARCODE, file=paste(opt$output,"/cell_barcodes.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
write.csv(tmp,file=paste(opt$output,"/../report/cutoff.csv",sep=""),row.names=FALSE,quote=FALSE)
write.csv(small,file=paste(opt$output,"/../report/vln.csv",sep=""),row.names=FALSE,quote=FALSE)
                                        #png(file=paste(opt$output,"/cell_count_summary.png",sep=""), width=700,height=300,res=100,pointsize=10)
png(file=paste(opt$output,"/cell_count_summary.png",sep=""), width=1200,height=400,res=80)
p1 = ggplot(tmp,aes(x=barcodes,y=UMI)) + xlim(0,10) + ylim(0,9)
p1 = p1 +annotate("text",x=0.2,y=9,label="Estimated Number of Cell:",size=10,hjust=0)
p1 = p1 +annotate("text",x=3,y=6.5,label=CellNum,size=20,hjust=0, colour="red")
p1 = p1 +annotate("text",x=0.2,y=3.5,label="Reads in cell:",size=6,hjust=0)
p1 = p1 +annotate("text",x=9,y=3.5,label=round(in_cell,3),size=6,hjust=1)
p1 = p1 +annotate("segment",x = 0, xend = 10, y = 3, yend = 3,colour = "blue")
p1 = p1 +annotate("text",x=0.2,y=2.5,label="Mean UMI counts per cell:",size=6,hjust=0)
p1 = p1 +annotate("text",x=9,y=2.5,label=round(UMI_mean),size=6,hjust=1)
p1 = p1 +annotate("segment", x = 0, xend = 10, y = 2, yend = 2,colour = "blue")
p1 = p1 +annotate("text",x=0.2,y=1.5,label="Mean Genes per cell:",size=6,hjust=0)
p1 = p1 +annotate("text",x=9,y=1.5,label=round(Gene_mean),size=6,hjust=1)
p1 = p1 +annotate("segment",x  = 0, xend  = 10, y  = 1, yend = 1,colour = "blue")
p1 = p1 +theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(),panel.background = element_blank(), panel.grid.major=element_blank(),panel.grid.minor = element_blank())

p = ggplot(tmp,aes(x=barcodes,y=UMI))
p = p + geom_line(aes(color=cell),size=2) +scale_color_manual(values=c("#999999","blue"))
p = p + scale_x_log10(name="Barcodes",breaks=c(1,10,100,1000,10000,100000),labels=c(1,10,100,"1k","10K","100K"))
p = p + scale_y_log10(name="UB",breaks=c(1,10,100,1000,10000,100000),labels=c(1,10,100,"1k","10K","100K"))
p = p + theme_bw() + geom_vline(xintercept =cutoff)
#p = p + geom_text(aes(x=10,y=10,label = paste("cell=",cutoff)), color = 'blue',size=6)
#p = p + geom_text(aes(x=10,y=4,label = paste("UB=",m)), color = 'blue',size=6)
#p = p + geom_text(aes(x=10,y=1,label = paste("nUIM%=",in_cell)), color = 'blue',size=6)
p = p + theme(legend.position = "none")

                                        #p1 <- ggplot(bc) + geom_boxplot(aes(x=5,y=UB), outlier.shape = 8, width=10) + theme_classic(
#p3 <- p3 + geom_jitter(aes(x=sample(1:10,nrow(bc),replace = T),y=UB),alpha=0.2,color="blue") #+ scale_y_log10()
p3 <- ggplot(cc) + geom_violin(aes(x=5,y=UB),stat="ydensity") + theme_classic() +geom_jitter(aes(x=5,y=UB),alpha=0.2,color="blue")
p3 <- p3 + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
#p2 <- ggplot(bc) + geom_boxplot(aes(x=5,y=GN), outlier.shape = 8, width=10) + theme_classic()
                                        #p2 <- p2 + geom_jitter(aes(x=sample(1:10,nrow(bc),replace = T),y=GN),alpha=0.2,color="blue") #+ scale_y_log10()
p2 <- ggplot(cc) + geom_violin(aes(x=5,y=GN),stat="ydensity") + theme_classic() +geom_jitter(aes(x=5,y=GN),alpha=0.2,color="blue")
p2 <- p2 + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())

pict <- plot_grid(p1, p, p3,p2, rel_widths = c(4,3,1.5,1.5),ncol=4)
pict
dev.off()
