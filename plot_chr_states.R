# Basic output plotting script which reads in each txt output file and
# generates a PNG file of chromatin states per text file

in.files <- dir(pattern="^chr.*txt$")
 
for(i in in.files){
 
	print(i)
    
	y<-read.fwf(i,widths=c(4,20,20,20,20),as.is=T); x <- y[-1,]; colnames(x)<-as.character(y[1,])

	for(j in 1:5){
		x[,j]<-as.numeric(x[,j])
	}
	
	plot.name <- gsub("txt","png",i)
	
	png(plot.name,height=300,width=700)
	
	matplot(x[,1]-107,x[,3:5],xlab="Lattice unit",ylab="Chromatin state",col=c("black","red","blue"),type="l",lty=1)
	legend("topright",lty=1,col=c("black","red","blue"),legend=c("Het","Eu","CTCF"))
	
	dev.off()
	
}

