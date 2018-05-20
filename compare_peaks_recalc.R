# This R script computes the peak size of the model and compares it with the
# data

# BED file has each region and the position of each peak within
# that region

peaks<-read.table("HET_6387_centred_20100_full.bed",header=F,as.is=T)

# This file has the following structure:
#
# Adjust this for whichever BED file you use
#
# col 1-3 extended intervals
# col 4 name
# col 5 ""
# col 6 ""
# col 7-9 peak positions
# col 10 same as col 4
# col 11,12 ""

# Add two columns to hold the model peak size
peaks$V13 <- NA
peaks$V14 <- NA

# Set column names for the output file, based on the BED file
colnames(peaks) <- c("Chromosome","Left.ext","Right.ext","Name",".",".","Chromosome","Left.peak","Right.peak","Name",".",".","Left.model","Right.model")

# Set a threshold for setting which p(heterochromatin) means a peak
peaks_threshold <- 0.5 

# Set the lattice size of the computation and the number of lattice units
lattice_size <- 161
no_nucs <- floor((ceiling(40201/lattice_size))/2)*2+1
reg_centre <-ceiling(no_nucs/2)

# Go through each peak one by one
for(i in 1:dim(peaks)[1]){
  # set the filename based on the input region
  filename<-paste0(peaks[i,1],".",peaks[i,2],"-",peaks[i,3],".txt");
  print(filename)
  
  # Try reading in the file
  x <- try(  
    model <- read.table(filename,header=T,as.is=T)    
  )
  
  if(class(x)!="try-error"){ 
    # if we do have an output file
    # it has the following structure:
    # col 1 is site no
    # col 2 is HP1 binding
    # col 3 is State 1 (Heterochromatin) prob
    # col 4 is State 2 (Euchromatin) prob
    # col 5 is CTCF binding site
    
    # Read outwards from the centre of the region
    # to get the peak size in terms of lattice
    # units 
    centre<-reg_centre
    
    left <- NA  
    right <- NA
    
    if(model[centre,3]>peaks_threshold){ 
      # if a heterchromatin region exists at the centre
      left<-reg_centre
      right<-reg_centre
      
      # start a search in either direction to get the maximum size of the region
      while(model[left,3]>peaks_threshold){
        left <- left -1
        if(left==0) break; # we have reached the end of the region
      }
      while(model[right,3]>peaks_threshold){
        right <- right +1
        if(right>length(model[,2])) break; # we have reached the end of the region
      }
      
    } 
    
    # If a peak is found left is unit to left, right is unit to right so need to add/subtract 1
    
    # add to data file
    # note that the output is in lattice units, not in bp - so MUST be scaled later
    peaks[i,12:13]<-c(left+1,right-1)
    
  }
  
  
}

# write file to text file, which can be further processed
write.table(peaks,file="Peaks_file_recalc.txt",sep="\t",quote=F,row.names=F,col.names=T)
