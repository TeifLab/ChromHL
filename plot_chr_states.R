# Basic output plotting script which reads in each txt output file and
# generates a PNG file of chromatin states per text file

in.files <- dir(pattern="^chr.*txt$")

# based on lattice size, calculate the centre of the region in terms of the lattice
lattice_size <- 161
no_nucs <- floor((ceiling(40201/lattice_size))/2)*2+1
reg_centre <-ceiling(no_nucs/2)

# for each input file
for(i in in.files){
 
	print(i)
    
	# read in the fixed-width txt file
	y<-read.fwf(i,widths=c(4,20,20,20,20),as.is=T); x <- y[-1,]; colnames(x)<-as.character(y[1,])

	# convert to numeric values rather than strings
	for(j in 1:5){
		x[,j]<-as.numeric(x[,j])
	}
	
	# Output PNG file is same as input file but with .png extension
	plot.name <- gsub("txt","png",i)
	
	# Plot chromatin state map
	png(plot.name,height=300,width=700)
	
	matplot(x[,1]-reg_centre,x[,3:5],xlab="Lattice unit",ylab="Chromatin state",col=c("black","red","blue"),type="l",lty=1)
	legend("topright",lty=1,col=c("black","red","blue"),legend=c("Het","Eu","CTCF"))
	
	dev.off()
	
}

