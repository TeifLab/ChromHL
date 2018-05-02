#
# Script for computing peak size of model and compare with data
#
#
#
#

# Read in data from bed file

peaks<-read.table("HET_6387_centred_20100_full.bed",header=F,as.is=T)

# col 1-3 extended intervals
# col 4 name
# col 5 ""
# col 6 ""
# col 7-9 peak positions
# col 10 same as col 4
# col 11,12 ""
# col 13,14

peaks$V13 <- NA # will hold peaks size
peaks$V14 <- NA

colnames(peaks) <- c("Chromosome","Left.ext","Right.ext","Name",".",".","Chromosome","Left.peak","Right.peak","Name",".",".","Left.model","Right.model")

peaks_threshold <- 0.5 # simple classifier takes the threshold to be 0.5

lattice_size <- 161
no_nucs <- floor((ceiling(40201/lattice_size))/2)*2+1
reg_centre <-ceiling(no_nucs/2)

for(i in 1:dim(peaks)[1]){
  filename<-paste0(peaks[i,1],".",peaks[i,2],"-",peaks[i,3],".txt");
  
  print(filename)
  
  x <- try(
    
    model <- read.table(filename,header=T,as.is=T) 
    
  )
  
  if(class(x)!="try-error"){
    
    # col 1 is site no
    # col 2 is HP1 binding
    # col 3 is State 1 (Heterochromatin) prob
    # col 4 is State 2 (Euchromatin) prob
    # col 5 is CTCF binding site
    
    centre<-reg_centre
    
    left <- NA  
    right <- NA
    if(model[centre,3]>peaks_threshold){
      
      left<-reg_centre
      right<-reg_centre
      
      while(model[left,3]>peaks_threshold){
        left <- left -1
        if(left==0) break;
      }
      while(model[right,3]>peaks_threshold){
        right <- right +1
        if(right>length(model[,2])) break;
      }
      
    } # if peak found left is nuc to left, right is nuc to right so need to add/subtract 1
    
    peaks[i,12:13]<-c(left+1,right-1)
    
  }
  
  
}

write.table(peaks,file="Peaks_file_recalc.txt",sep="\t",quote=F,row.names=F,col.names=T)


