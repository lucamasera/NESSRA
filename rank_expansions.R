rm(list=ls()) #clear all vars

## MC4 ranker function as implemented last year
## MC4: da script del dott. Argentini
#data frame - the columns would be represented by a number (1, the highest rank, 
# to some number which is the lowerst rank), and each row is the ranking of a different network.
# the contents of each cell is a gene name.
# col_discarded is how many cols to discard in case you have other data in the begining
# of the columns before the ranks. K and alpha are the same k and alpha from the paper
# (I used k=10 and alpha=0.05, as was found to be the best in the paper)

mc4_ranker  <-  function(dataframe,col_discarded=0,k_max=ncol(dataframe) - col_discarded, alpha=0.05){
  dfl <- nrow(dataframe)
  x <- y <- vector("character")       
  for(j in 1:k_max){                 
    for(i in 1:dfl){                                     
      x[(i-1)*k_max + j] <- as.character(dataframe[i,col_discarded+j])
    }
  }

  y <- unique(x)  
  yy <- length(y)
  dfl <- nrow(dataframe)  
  
  GeneRank <- matrix(0,yy,dfl)
  rownames(GeneRank) <- y
  colnames(GeneRank) <- c(1:dfl)
  
  ind <- vector("numeric")
  for(j in 1:(dfl+1)){
    ind[j] <- (j-1)*k_max
  }
  
  for(i in 1:yy){
    for(j in 1:dfl){
      if (is.na(which(y[i] == x[(ind[j]+1):ind[j+1]])[1]%%k_max) == TRUE){
        GeneRank[i,j] <- k_max
      } else if (which(y[i] == x[(ind[j]+1):ind[j+1]])[1]%%k_max == 0){
        GeneRank[i,j] <- k_max - 1
      } else {
        GeneRank[i,j] <- which(y[i]==x[(ind[j]+1):ind[j+1]])[1]%%k_max -1
      }
    }
  }
  
  obj= as.character(y) 
  t_pro_mat<- matrix(0,ncol=yy,nrow=yy)     
  for (i in 1:yy){
    for (j in 1:yy){
      if (i != j ){
        kk=1
        good_count= 0
        while ( kk <= dfl){
          if (GeneRank[i,kk] < GeneRank[j,kk]){
            good_count= good_count+1
          }
          if (good_count >=round(dfl/2) ){ 
            t_pro_mat[i,j]= 1/yy
            break
          }
          kk=kk+1
          
        }      
      }
    }
    
  } 
  for (j in 1:yy){
    t_pro_mat[j,j]= 1 - sum(t_pro_mat[j,])
  }
  a=alpha
  t_p_mat_t = ((1-a)*t_pro_mat ) + a/yy
  MC4 <- new("markovchain", states = y, transitionMatrix = t_p_mat_t, name = "rank_aggr")
  result <- data.frame(key=y, rank=rank(round(steadyStates(MC4),digits=10),ties.method="average"))   
  result=result[order(result$rank),]
  row.names(result) <- NULL
  return(result)
}

# code for running the function over all expansions. No diffrence if you use it or go manually, 
#but if you do then have each network in it's own folder, and each experiment results
# in a file ending with .expansion, first column being ranks and second column being gene names.
args <- commandArgs(TRUE)

work_dir=args[1]; #put the path to the expansion folder here
require(markovchain)

setwd(work_dir)
d=list.dirs(recursive=FALSE) # get netowrk folders
expan_array <-
ranks <- list()

print(d)

max.len=0
min.len=strtoi(args[2])
expan_data <- list()

for (i in 1:(length(d))){ #loop over the network folders
  f=list.files(d[i]) #get all the files in network folder i
  exp_ind=grep('expansion',f)  #get index for all the files in network folder that includes 'expansion' in their name
  exp_file = paste(d[i],'/',f[exp_ind[1]],sep='')
  print(exp_file)
  temp=read.csv(exp_file,header=TRUE,skip=strtoi(args[3]))
  expan_data[[i]]<- data.frame(lapply(temp$node, as.character), stringsAsFactors=FALSE)
  max.len=max(max.len,length(expan_data[[i]]))
  min.len=min(min.len,length(expan_data[[i]]))
}

x <- list()
for (j in 1:length(d)){
  x[[j]]= c(expan_data[[j]], rep(NA, max.len - length(expan_data[[j]]))) #for the array to be even, path with NAs.
}
x=do.call(rbind, x)

ranks <- mc4_ranker(x,0,min.len) # call the rankers the third parameter is how many genes we want in the rank
rm(x,expan_data,f,temp,max.len)

write.csv(ranks, 'mc4_aggregation.csv', row.names = FALSE, quote = FALSE) # write results to file
