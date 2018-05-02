x <- read.table("Peaks_file_recalc.txt",as.is=T,header=T)

y <- data.frame(peak.size = x$Right.peak-x$Left.peak,
                Model.size = x$Right.model-x$Left.model)

peak_size<-161 # change this for lattice size

y$Model.size <- y$Model.size * peak_size
Model.peaks <- sum(!is.na(y$Model.size))
y$Model.size[y$Model.size<200] <- NA

png("histograms.png")

hist(y$peak.size,prob=TRUE,xlim=c(0,15000),breaks=seq(from=0,to=peak_size*1e4,by=peak_size),ylim=c(0,0.002),col="#FF000030",main="Peak size histogram",xlab="Peak size (bp)")
hist(y$Model.size,prob=TRUE,breaks=seq(from=0,to=peak_size*1e4,by=peak_size),col="#00FF0030",add=T)
#hist(y$X0.01.size,prob=TRUE,breaks=seq(from=0,to=peak_size*1e4,by=peak_size),col="#0000FF30",add=T)

legend("topright",
       c("Data peaks","Model peaks"),
       lty=0,
       fill=c("#FF000030","#00FF0030"))

dev.off()


png("histograms.zoom.png")

hist(y$peak.size,prob=TRUE,xlim=c(0,15000),breaks=seq(from=0,to=peak_size*1e4,by=peak_size),ylim=c(0,0.00025),col="#FF000030",main="Peak size histogram",xlab="Peak size (bp)")
hist(y$Model.size,prob=TRUE,breaks=seq(from=0,to=peak_size*1e4,by=peak_size),col="#00FF0030",add=T)
#hist(y$X0.01.size,prob=TRUE,breaks=seq(from=0,to=peak_size*1e4,by=peak_size),col="#0000FF30",add=T)

legend("topright",
       c("Data peaks","Model peaks"),
       lty=0,
       fill=c("#FF000030","#00FF0030"))

dev.off()

write.table(y$Model.size[!is.na(y$Model.size)],file="model_peaks.txt",col.names=F,row.names=F)
write.table(y$peak.size[!is.na(y$peak.size)],file="MACS2_peaks.txt",col.names=F,row.names=F)






