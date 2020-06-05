#!/usr/bin/env Rscript

suppressMessages({
    library(ggplot2)
    library(getopt)
    library(data.table)
    library(cowplot)
    library(DropletUtils)
})

arg<-matrix(c("input", "i","1","character","Path of input file",
              "output","o","2","character","Path of output directory",
              "force","f","1","numeric","Force cells number to analysis. You must set -f or -e to process barcode list.",
              "expect","e","1","numeric","Number of expect cells",
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
low = 50
if (!is.null(opt$expect)) {
    expect=as.numeric(opt$expect)
}

if (!is.null(opt$low)) {
    low=as.numeric(opt$low)
}

if (expect >0) {
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

if (!is.null(opt$force)) {
    force = as.numeric(opt$force)
    if (force > 0) {
        expect = force
        cutoff = expect
        m = sor[cutoff]
    }
}

tmp<-data.frame(barcodes=1:len,UMI=sor,cell=c(rep("true",cutoff),rep("noise",len-cutoff)))
cutoff=nrow(bc[which(bc$UB>=m),])

cc <- bc[1:cutoff,]
Species1=substr(colnames(cc)[5],1,nchar(colnames(cc)[5])-3)
Species2=substr(colnames(cc)[7],1,nchar(colnames(cc)[7])-3)
called = sum(cc$UB)
in_cell <- round(called/sum(bc$UB),digits=4)
in_cell_Human <- round(sum(bc[which(bc$UB>=m),][,paste(Species1,"_UB",sep="")])/sum(bc[,paste(Species1,"_UB",sep="")]),digits=4)
in_cell_Mouse <- round(sum(bc[which(bc$UB>=m),][,paste(Species2,"_UB",sep="")])/sum(bc[,paste(Species2,"_UB",sep="")]),digits=4)

small<-bc[which(bc$UB>=m),]
CellNum = nrow(small)
cell_mean = mean(small$Raw)
UMI_mean = mean(small$UB)
Gene_mean = mean(small$GN)

small$total = small[,paste(Species1,"_UB",sep="")] + small[,paste(Species2,"_UB",sep="")]
small = small[order(small$total, decreasing = T),]

small$Species = "Mix"
small[small[,paste(Species1,"_UB",sep="")]/small$total<0.1,]$Species=Species2
small[small[,paste(Species2,"_UB",sep="")]/small$total<0.1,]$Species=Species1

small$Species = as.factor(small$Species)
sum = summary(small$Species)
umi_mix_ratio = round(sum['Mix']/sum(sum),digits=4)
small1 <- subset(small, small$Species==Species2)
small2 <- subset(small, small$Species==Species1)
small1$UB = small1[,paste(Species2,"_UB",sep="")]
small1$GN = small1[,paste(Species2,"_GN",sep="")]
small2$UB = small2[,paste(Species1,"_UB",sep="")]
small2$GN = small2[,paste(Species1,"_GN",sep="")]
small.new = rbind(small1,small2)

CellNum_Human = nrow(small2)
CellNum_Mouse = nrow(small1)
cell_mean_Human = mean(small2$Raw)
cell_mean_Mouse = mean(small1$Raw)
UMI_median_Human = median(small2$Human_UB)
UMI_median_Mouse = median(small1$Mouse_UB)
Gene_median_Human = median(small2$Human_GN)
Gene_median_Mouse = median(small1$Mouse_GN)
cat(paste("Estimated Number of Cell,", CellNum, "\n",sep=""),file=paste(opt$output,"/cell_report.csv",sep=""))
cat(paste("Estimated Number of Cell (", Species1, "),", CellNum_Human, "\n",sep=""),file=paste(opt$output,"/cell_report.csv",sep=""), append=T)
cat(paste("Estimated Number of Cell (", Species2, "),", CellNum_Mouse, "\n",sep=""),file=paste(opt$output,"/cell_report.csv",sep=""),append=T)
cat(paste("Reads in cell,", round(in_cell,3)*100, "%\n",sep=""),file=paste(opt$output,"/cell_report.csv",sep=""), append=T)
cat(paste("Reads in cell (", Species1, "),", round(in_cell_Human,3)*100, "%\n",sep=""),file=paste(opt$output,"/cell_report.csv",sep=""), append=T)
cat(paste("Reads in cell (", Species2, "),", round(in_cell_Mouse,3)*100, "%\n",sep=""),file=paste(opt$output,"/cell_report.csv",sep=""), append=T)
cat(paste("Mean reads per cell,", round(cell_mean),"\n",sep=""),file=paste(opt$output,"/cell_report.csv",sep=""), append=T)
cat(paste("Mean reads per cell (", Species1, "),", round(cell_mean_Human),"\n",sep=""),file=paste(opt$output,"/cell_report.csv",sep=""), append=T)
cat(paste("Mean reads per cell (", Species2, "),", round(cell_mean_Mouse),"\n",sep=""),file=paste(opt$output,"/cell_report.csv",sep=""), append=T)
cat(paste("Median UMI Counts per Cell (",Species1,"),", round(UMI_median_Human),"\n",sep=""),file=paste(opt$output,"/cell_report.csv",sep=""), append=T)
cat(paste("Median UMI Counts per Cell (",Species2,"),", round(UMI_median_Mouse),"\n",sep=""),file=paste(opt$output,"/cell_report.csv",sep=""), append=T)
cat(paste("Median Genes per Cell (",Species1,"),", round(Gene_median_Human),"\n",sep=""),file=paste(opt$output,"/cell_report.csv",sep=""), append=T)
cat(paste("Median Genes per Cell (",Species2,"),", round(Gene_median_Mouse),"\n",sep=""),file=paste(opt$output,"/cell_report.csv",sep=""), append=T)
cat(paste("Single cell,", CellNum, "\n",sep=""),file=paste(opt$output,"/mix_report.csv",sep=""))
cat(paste("Mix Cell,", sum['Mix'], "\n",sep=""),file=paste(opt$output,"/mix_report.csv",sep=""), append=T)
cat(paste("Multiplet rate,", round(umi_mix_ratio,3)*100, "%\n",sep=""),file=paste(opt$output,"/mix_report.csv",sep=""), append=T)
write.csv(small[,c("BARCODE","Raw","Human_UB","Human_GN","Mouse_UB","Mouse_GN","Species")],file=paste(opt$output,"/vln.csv",sep=""),row.names=FALSE,quote=FALSE)
write.csv(tmp,file=paste(opt$output,"/cutoff.csv",sep=""),row.names=FALSE,quote=FALSE)
write.table(cc$BARCODE, file=paste(opt$output,"/cell_barcodes.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)

png(file=paste(opt$output,"/mixture_cells.png",sep=""), width=1200,height=400,res=80)
x = max(small$Human_UB)*.5
y = max(small$Mouse_UB)*.75
p1 <- ggplot(small, aes(x=Human_UB,y=Mouse_UB,color=Species)) + geom_point()  +scale_color_manual(values=c("Human" = "blue", "Mouse" =  "red", "Mix"="#999999")) + labs(title=paste("Multiplet rate = ", round(umi_mix_ratio,3)*100, "%", sep="")) + xlab(paste("Human cells : ",sum['Human'],sep="")) + ylab(paste("Mouse cells : ",sum['Mouse'])) + annotate("text",x=x,y=y, label=paste("Mix cells : ", sum['Mix']), size=10,color="#999999") + theme(plot.title = element_text(color="red", size=24, face="bold.italic"), panel.background = element_blank(),axis.line = element_line(colour = "black"), axis.title.x = element_text(color="blue", size=24, face="bold.italic"), axis.title.y = element_text(color="red", size=24, face="bold.italic"))
p2 <- ggplot(small.new, aes(y=UB,x=Species)) +geom_violin(stat="ydensity") + theme_classic() +geom_jitter(alpha=0.2,color="blue") + theme(axis.title.x=element_blank(),axis.title.y=element_text(size=20, face="bold.italic"),axis.text.x= element_text(size=14, face="bold.italic"),axis.ticks.x=element_blank()) + ylab("UMI")
p3 <- ggplot(small.new, aes(y=GN,x=Species)) +geom_violin(stat="ydensity") + theme_classic() +geom_jitter(alpha=0.2,color="blue") + theme(axis.title.x=element_blank(),axis.title.y=element_text(size=20, face="bold.italic"),axis.text.x= element_text(size=14, face="bold.italic"),axis.ticks.x=element_blank()) + ylab("Genes")
p4 <- ggplot(small.new, aes(y=Raw,x=Species)) +geom_violin(stat="ydensity") + theme_classic() +geom_jitter(alpha=0.2,color="blue") + theme(axis.title.x=element_blank(),axis.title.y=element_text(size=20, face="bold.italic"),axis.text.x= element_text(size=14, face="bold.italic"),axis.ticks.x=element_blank()) + ylab("Raw reads")
pict1 <- plot_grid(p1, p2, p3,p4, rel_widths = c(4,2,2,2),ncol=4)
pict1
dev.off()

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
p = p + theme(legend.position = "none")
p3 <- ggplot(cc) + geom_violin(aes(x=5,y=UB),stat="ydensity") + theme_classic() +geom_jitter(aes(x=5,y=UB),alpha=0.2,color="blue")
p3 <- p3 + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
p2 <- ggplot(cc) + geom_violin(aes(x=5,y=GN),stat="ydensity") + theme_classic() +geom_jitter(aes(x=5,y=GN),alpha=0.2,color="blue")
p2 <- p2 + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
pict2 <- plot_grid(p1, p, p3,p2, rel_widths = c(4,3,1.5,1.5),ncol=4)
pict2
dev.off()
