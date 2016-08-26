library("cluster")
library("factoextra")
library("optparse")
library("stringr")

# select the cancer type from the shell script
option_list = list(
  make_option(c("-c", "--cancer"), type = "character", default=NULL,
              help="cancer set abbreviation name", metavar = "character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
cancer_name = opt[[1]]

# The first part of the function makes the first row
# predicted = p
# residual = r
# what we need: 
load(paste0(cancer_name, "-rsqrs-", "1",".RData")) # get one of the runs together that way you can get the dimensions right
bottom_cloud.ol <- finalTree.ol
bottomcloud_array.ol <- array(rep(0, nrow(bottom_cloud.ol)*64*100), dim=c(100,nrow(bottom_cloud.ol),64)) #3d array of 0's soing to concatenate everything together

for(j in 1:100){ #for every run 
  
  load(paste0(cancer_name, "-rsqrs-", as.character(k),".RData")) #load the rsqrs file for the j'th run
  bottom_cloud.ol <- as.data.frame(bottom_cloud.ol) #be able to select colums
  pppppp_path <- c()
  pppprp_path <- c()
  ppprpp_path <- c()
  ppprrp_path <- c()
  pprppp_path <- c()
  pprprp_path <- c()
  pprrpp_path <- c()
  pprrrp_path <- c()
  prpppp_path <- c()
  prpprp_path <- c()
  prprpp_path <- c()
  prprrp_path <- c()
  prrppp_path <- c()
  prrprp_path <- c()
  prrrpp_path <- c()
  prrrrp_path <- c()
  rppppp_path <- c()
  rppprp_path <- c()
  rpprpp_path <- c()
  rpprrp_path <- c()
  rprppp_path <- c()
  rprprp_path <- c()
  rprrpp_path <- c()
  rprrrp_path <- c()
  rrpppp_path <- c()
  rrpprp_path <- c()
  rrprpp_path <- c()
  rrprrp_path <- c()
  rrrppp_path <- c()
  rrrprp_path <- c()
  rrrrpp_path <- c()
  rrrrrp_path <- c()
  #############################################
  pppppr_path <- c()
  pppprr_path <- c()
  ppprpr_path <- c()
  ppprrr_path <- c()
  pprppr_path <- c()
  pprprr_path <- c()
  pprrpr_path <- c()
  pprrrr_path <- c()
  prpppr_path <- c()
  prpprr_path <- c()
  prprpr_path <- c()
  prprrr_path <- c()
  prrppr_path <- c()
  prrprr_path <- c()
  prrrpr_path <- c()
  prrrrr_path <- c()
  rppppr_path <- c()
  rppprr_path <- c()
  rpprpr_path <- c()
  rpprrr_path <- c()
  rprppr_path <- c()
  rprprr_path <- c()
  rprrpr_path <- c()
  rprrrr_path <- c()
  rrpppr_path <- c()
  rrpprr_path <- c()
  rrprpr_path <- c()
  rrprrr_path <- c()
  rrrppr_path <- c()
  rrrprr_path <- c()
  rrrrpr_path <- c()
  rrrrrr_path <- c()
  print("a") 
  
  for(i in 1:nrow(bottom_cloud.ol)){ #for all of the genes 
	# below you just get the product of each element to get a combination. I you look at the tree plot
	# then you get every "a" element for each gene and concatenate it together at the end into a 3d matrix
	# dimension 1 = gene, dimension 2 = combination, dimension 3 = run
    pppppp_path[i] <- bottom_cloud.ol$p[i] * bottom_cloud.ol$pp[i] * bottom_cloud.ol$ppp[i] * bottom_cloud.ol$pppp[i] * bottom_cloud.ol$ppppp[i] * bottom_cloud.ol$pppppp[i]
    pppppr_path[i] <- bottom_cloud.ol$p[i] * bottom_cloud.ol$pp[i] * bottom_cloud.ol$ppp[i] * bottom_cloud.ol$pppp[i] * bottom_cloud.ol$ppppp[i] * bottom_cloud.ol$pppppr[i]
    
    pppprp_path[i] <- bottom_cloud.ol$p[i] * bottom_cloud.ol$pp[i] * bottom_cloud.ol$ppp[i] * bottom_cloud.ol$pppp[i] * bottom_cloud.ol$ppppr[i] * bottom_cloud.ol$pppprp[i]
    pppprr_path[i] <- bottom_cloud.ol$p[i] * bottom_cloud.ol$pp[i] * bottom_cloud.ol$ppp[i] * bottom_cloud.ol$pppp[i] * bottom_cloud.ol$ppppr[i] * bottom_cloud.ol$pppprr[i]
    
    ppprpp_path[i] <- bottom_cloud.ol$p[i] * bottom_cloud.ol$pp[i] * bottom_cloud.ol$ppp[i] * bottom_cloud.ol$pppr[i] * bottom_cloud.ol$ppprp[i] * bottom_cloud.ol$ppprpp[i]
    ppprpr_path[i] <- bottom_cloud.ol$p[i] * bottom_cloud.ol$pp[i] * bottom_cloud.ol$ppp[i] * bottom_cloud.ol$pppr[i] * bottom_cloud.ol$ppprp[i] * bottom_cloud.ol$ppprpr[i]
    
    ppprrp_path[i] <- bottom_cloud.ol$p[i] * bottom_cloud.ol$pp[i] * bottom_cloud.ol$ppp[i] * bottom_cloud.ol$pppr[i] * bottom_cloud.ol$ppprr[i] * bottom_cloud.ol$ppprrp[i]
    ppprrr_path[i] <- bottom_cloud.ol$p[i] * bottom_cloud.ol$pp[i] * bottom_cloud.ol$ppp[i] * bottom_cloud.ol$pppr[i] * bottom_cloud.ol$ppprr[i] * bottom_cloud.ol$ppprrr[i]
    
    pprppp_path[i] <- bottom_cloud.ol$p[i] * bottom_cloud.ol$pp[i] * bottom_cloud.ol$ppr[i] * bottom_cloud.ol$pprp[i] * bottom_cloud.ol$pprpp[i] * bottom_cloud.ol$pprppp[i]
    pprppr_path[i] <- bottom_cloud.ol$p[i] * bottom_cloud.ol$pp[i] * bottom_cloud.ol$ppr[i] * bottom_cloud.ol$pprp[i] * bottom_cloud.ol$pprpp[i] * bottom_cloud.ol$pprppr[i]
    
    pprprp_path[i] <- bottom_cloud.ol$p[i] * bottom_cloud.ol$pp[i] * bottom_cloud.ol$ppr[i] * bottom_cloud.ol$pprp[i] * bottom_cloud.ol$pprpr[i] * bottom_cloud.ol$pprprp[i]
    pprprr_path[i] <- bottom_cloud.ol$p[i] * bottom_cloud.ol$pp[i] * bottom_cloud.ol$ppr[i] * bottom_cloud.ol$pprp[i] * bottom_cloud.ol$pprpr[i] * bottom_cloud.ol$pprprr[i]
    
    pprrpp_path[i] <- bottom_cloud.ol$p[i] * bottom_cloud.ol$pp[i] * bottom_cloud.ol$ppr[i] * bottom_cloud.ol$pprr[i] * bottom_cloud.ol$pprrp[i] * bottom_cloud.ol$pprrpp[i]
    pprrpr_path[i] <- bottom_cloud.ol$p[i] * bottom_cloud.ol$pp[i] * bottom_cloud.ol$ppr[i] * bottom_cloud.ol$pprr[i] * bottom_cloud.ol$pprrp[i] * bottom_cloud.ol$pprrpr[i]
    
    pprrrp_path[i] <- bottom_cloud.ol$p[i] * bottom_cloud.ol$pp[i] * bottom_cloud.ol$ppr[i] * bottom_cloud.ol$pprr[i] * bottom_cloud.ol$pprrr[i] * bottom_cloud.ol$pprrrp[i]
    pprrrr_path[i] <- bottom_cloud.ol$p[i] * bottom_cloud.ol$pp[i] * bottom_cloud.ol$ppr[i] * bottom_cloud.ol$pprr[i] * bottom_cloud.ol$pprrr[i] * bottom_cloud.ol$pprrrr[i]
    
    prpppp_path[i] <- bottom_cloud.ol$p[i] * bottom_cloud.ol$pr[i] * bottom_cloud.ol$prp[i] * bottom_cloud.ol$prpp[i] * bottom_cloud.ol$prppp[i] * bottom_cloud.ol$prpppp[i]
    prpppr_path[i] <- bottom_cloud.ol$p[i] * bottom_cloud.ol$pr[i] * bottom_cloud.ol$prp[i] * bottom_cloud.ol$prpp[i] * bottom_cloud.ol$prppp[i] * bottom_cloud.ol$prpppr[i]
    
    prpprp_path[i] <- bottom_cloud.ol$p[i] * bottom_cloud.ol$pr[i] * bottom_cloud.ol$prp[i] * bottom_cloud.ol$prpp[i] * bottom_cloud.ol$prppr[i] * bottom_cloud.ol$prpprp[i]
    prpprr_path[i] <- bottom_cloud.ol$p[i] * bottom_cloud.ol$pr[i] * bottom_cloud.ol$prp[i] * bottom_cloud.ol$prpp[i] * bottom_cloud.ol$prppr[i] * bottom_cloud.ol$prpprr[i]
    
    prprpp_path[i] <- bottom_cloud.ol$p[i] * bottom_cloud.ol$pr[i] * bottom_cloud.ol$prp[i] * bottom_cloud.ol$prpr[i] * bottom_cloud.ol$prprp[i] * bottom_cloud.ol$prprpp[i]
    prprpr_path[i] <- bottom_cloud.ol$p[i] * bottom_cloud.ol$pr[i] * bottom_cloud.ol$prp[i] * bottom_cloud.ol$prpr[i] * bottom_cloud.ol$prprp[i] * bottom_cloud.ol$prprpr[i]
    
    prprrp_path[i] <- bottom_cloud.ol$p[i] * bottom_cloud.ol$pr[i] * bottom_cloud.ol$prp[i] * bottom_cloud.ol$prpr[i] * bottom_cloud.ol$prprr[i] * bottom_cloud.ol$prprrp[i]
    prprrr_path[i] <- bottom_cloud.ol$p[i] * bottom_cloud.ol$pr[i] * bottom_cloud.ol$prp[i] * bottom_cloud.ol$prpr[i] * bottom_cloud.ol$prprr[i] * bottom_cloud.ol$prprrr[i]
    
    prrppp_path[i] <- bottom_cloud.ol$p[i] * bottom_cloud.ol$pr[i] * bottom_cloud.ol$prr[i] * bottom_cloud.ol$prrp[i] * bottom_cloud.ol$prrpp[i] * bottom_cloud.ol$prrppp[i]
    prrppr_path[i] <- bottom_cloud.ol$p[i] * bottom_cloud.ol$pr[i] * bottom_cloud.ol$prr[i] * bottom_cloud.ol$prrp[i] * bottom_cloud.ol$prrpp[i] * bottom_cloud.ol$prrppr[i]
    
    prrprp_path[i] <- bottom_cloud.ol$p[i] * bottom_cloud.ol$pr[i] * bottom_cloud.ol$prr[i] * bottom_cloud.ol$prrp[i] * bottom_cloud.ol$prrpr[i] * bottom_cloud.ol$prrprp[i]
    prrprr_path[i] <- bottom_cloud.ol$p[i] * bottom_cloud.ol$pr[i] * bottom_cloud.ol$prr[i] * bottom_cloud.ol$prrp[i] * bottom_cloud.ol$prrpr[i] * bottom_cloud.ol$prrprr[i]
    
    prrrpp_path[i] <- bottom_cloud.ol$p[i] * bottom_cloud.ol$pr[i] * bottom_cloud.ol$prr[i] * bottom_cloud.ol$prrr[i] * bottom_cloud.ol$prrrp[i] * bottom_cloud.ol$prrrpp[i]
    prrrpr_path[i] <- bottom_cloud.ol$p[i] * bottom_cloud.ol$pr[i] * bottom_cloud.ol$prr[i] * bottom_cloud.ol$prrr[i] * bottom_cloud.ol$prrrp[i] * bottom_cloud.ol$prrrpr[i]
    
    prrrrp_path[i] <- bottom_cloud.ol$p[i] * bottom_cloud.ol$pr[i] * bottom_cloud.ol$prr[i] * bottom_cloud.ol$prrr[i] * bottom_cloud.ol$prrrr[i] * bottom_cloud.ol$prrrrp[i]
    prrrrr_path[i] <- bottom_cloud.ol$p[i] * bottom_cloud.ol$pr[i] * bottom_cloud.ol$prr[i] * bottom_cloud.ol$prrr[i] * bottom_cloud.ol$prrrr[i] * bottom_cloud.ol$prrrrr[i]
    
    rppppp_path[i] <- bottom_cloud.ol$r[i] * bottom_cloud.ol$rp[i] * bottom_cloud.ol$rpp[i] * bottom_cloud.ol$rppp[i] * bottom_cloud.ol$rpppp[i] * bottom_cloud.ol$rppppp[i]
    rppppr_path[i] <- bottom_cloud.ol$r[i] * bottom_cloud.ol$rp[i] * bottom_cloud.ol$rpp[i] * bottom_cloud.ol$rppp[i] * bottom_cloud.ol$rpppp[i] * bottom_cloud.ol$rppppr[i]
    
    rppprp_path[i] <- bottom_cloud.ol$r[i] * bottom_cloud.ol$rp[i] * bottom_cloud.ol$rpp[i] * bottom_cloud.ol$rppp[i] * bottom_cloud.ol$rpppr[i] * bottom_cloud.ol$rppprp[i]
    rppprr_path[i] <- bottom_cloud.ol$r[i] * bottom_cloud.ol$rp[i] * bottom_cloud.ol$rpp[i] * bottom_cloud.ol$rppp[i] * bottom_cloud.ol$rpppr[i] * bottom_cloud.ol$rppprr[i]
    
    rpprpp_path[i] <- bottom_cloud.ol$r[i] * bottom_cloud.ol$rp[i] * bottom_cloud.ol$rpp[i] * bottom_cloud.ol$rppr[i] * bottom_cloud.ol$rpprp[i] * bottom_cloud.ol$rpprpp[i]
    rpprpr_path[i] <- bottom_cloud.ol$r[i] * bottom_cloud.ol$rp[i] * bottom_cloud.ol$rpp[i] * bottom_cloud.ol$rppr[i] * bottom_cloud.ol$rpprp[i] * bottom_cloud.ol$rpprpr[i]
    
    
    rpprrp_path[i] <- bottom_cloud.ol$r[i] * bottom_cloud.ol$rp[i] * bottom_cloud.ol$rpp[i] * bottom_cloud.ol$rppr[i] * bottom_cloud.ol$rpprr[i] * bottom_cloud.ol$rpprrp[i]
    rpprrr_path[i] <- bottom_cloud.ol$r[i] * bottom_cloud.ol$rp[i] * bottom_cloud.ol$rpp[i] * bottom_cloud.ol$rppr[i] * bottom_cloud.ol$rpprr[i] * bottom_cloud.ol$rpprrr[i]  
    
    rprppp_path[i] <- bottom_cloud.ol$r[i] * bottom_cloud.ol$rp[i] * bottom_cloud.ol$rpr[i] * bottom_cloud.ol$rprp[i] * bottom_cloud.ol$rprpp[i] * bottom_cloud.ol$rprppp[i]
    rprppr_path[i] <- bottom_cloud.ol$r[i] * bottom_cloud.ol$rp[i] * bottom_cloud.ol$rpr[i] * bottom_cloud.ol$rprp[i] * bottom_cloud.ol$rprpp[i] * bottom_cloud.ol$rprppr[i]
    
    rprprp_path[i] <- bottom_cloud.ol$r[i] * bottom_cloud.ol$rp[i] * bottom_cloud.ol$rpr[i] * bottom_cloud.ol$rprp[i] * bottom_cloud.ol$rprpr[i] * bottom_cloud.ol$rprprp[i]
    rprprr_path[i] <- bottom_cloud.ol$r[i] * bottom_cloud.ol$rp[i] * bottom_cloud.ol$rpr[i] * bottom_cloud.ol$rprp[i] * bottom_cloud.ol$rprpr[i] * bottom_cloud.ol$rprprr[i]
    
    rprrpp_path[i] <- bottom_cloud.ol$r[i] * bottom_cloud.ol$rp[i] * bottom_cloud.ol$rpr[i] * bottom_cloud.ol$rprr[i] * bottom_cloud.ol$rprrp[i] * bottom_cloud.ol$rprrpp[i]
    rprrpr_path[i] <- bottom_cloud.ol$r[i] * bottom_cloud.ol$rp[i] * bottom_cloud.ol$rpr[i] * bottom_cloud.ol$rprr[i] * bottom_cloud.ol$rprrp[i] * bottom_cloud.ol$rprrpr[i]
    
    rprrrp_path[i] <- bottom_cloud.ol$r[i] * bottom_cloud.ol$rp[i] * bottom_cloud.ol$rpr[i] * bottom_cloud.ol$rprr[i] * bottom_cloud.ol$rprrr[i] * bottom_cloud.ol$rprrrp[i]
    rprrrr_path[i] <- bottom_cloud.ol$r[i] * bottom_cloud.ol$rp[i] * bottom_cloud.ol$rpr[i] * bottom_cloud.ol$rprr[i] * bottom_cloud.ol$rprrr[i] * bottom_cloud.ol$rprrrr[i]
    
    rrpppp_path[i] <- bottom_cloud.ol$r[i] * bottom_cloud.ol$rr[i] * bottom_cloud.ol$rrp[i] * bottom_cloud.ol$rrpp[i] * bottom_cloud.ol$rrppp[i] * bottom_cloud.ol$rrpppp[i]
    rrpppr_path[i] <- bottom_cloud.ol$r[i] * bottom_cloud.ol$rr[i] * bottom_cloud.ol$rrp[i] * bottom_cloud.ol$rrpp[i] * bottom_cloud.ol$rrppp[i] * bottom_cloud.ol$rrpppr[i]  
    
    rrpprp_path[i] <- bottom_cloud.ol$r[i] * bottom_cloud.ol$rr[i] * bottom_cloud.ol$rrp[i] * bottom_cloud.ol$rrpp[i] * bottom_cloud.ol$rrppr[i] * bottom_cloud.ol$rrpprp[i]
    rrpprr_path[i] <- bottom_cloud.ol$r[i] * bottom_cloud.ol$rr[i] * bottom_cloud.ol$rrp[i] * bottom_cloud.ol$rrpp[i] * bottom_cloud.ol$rrppr[i] * bottom_cloud.ol$rrpprr[i]  
    
    rrrppp_path[i] <- bottom_cloud.ol$r[i] * bottom_cloud.ol$rr[i] * bottom_cloud.ol$rrr[i] * bottom_cloud.ol$rrrp[i] * bottom_cloud.ol$rrrpp[i] * bottom_cloud.ol$rrrppp[i]
    rrrppr_path[i] <- bottom_cloud.ol$r[i] * bottom_cloud.ol$rr[i] * bottom_cloud.ol$rrr[i] * bottom_cloud.ol$rrrp[i] * bottom_cloud.ol$rrrpp[i] * bottom_cloud.ol$rrrppr[i]
    
    rrrprp_path[i] <- bottom_cloud.ol$r[i] * bottom_cloud.ol$rr[i] * bottom_cloud.ol$rrr[i] * bottom_cloud.ol$rrrp[i] * bottom_cloud.ol$rrrpr[i] * bottom_cloud.ol$rrrprp[i]
    rrrprr_path[i] <- bottom_cloud.ol$r[i] * bottom_cloud.ol$rr[i] * bottom_cloud.ol$rrr[i] * bottom_cloud.ol$rrrp[i] * bottom_cloud.ol$rrrpr[i] * bottom_cloud.ol$rrrprr[i]
    
    rrprpp_path[i] <- bottom_cloud.ol$r[i] * bottom_cloud.ol$rr[i] * bottom_cloud.ol$rrp[i] * bottom_cloud.ol$rrpr[i] * bottom_cloud.ol$rrprp[i] * bottom_cloud.ol$rrprpp[i]
    rrprpr_path[i] <- bottom_cloud.ol$r[i] * bottom_cloud.ol$rr[i] * bottom_cloud.ol$rrp[i] * bottom_cloud.ol$rrpr[i] * bottom_cloud.ol$rrprp[i] * bottom_cloud.ol$rrprpr[i]
    
    rrprrp_path[i] <- bottom_cloud.ol$r[i] * bottom_cloud.ol$rr[i] * bottom_cloud.ol$rrp[i] * bottom_cloud.ol$rrpr[i] * bottom_cloud.ol$rrprr[i] * bottom_cloud.ol$rrprrp[i]
    rrprrr_path[i] <- bottom_cloud.ol$r[i] * bottom_cloud.ol$rr[i] * bottom_cloud.ol$rrp[i] * bottom_cloud.ol$rrpr[i] * bottom_cloud.ol$rrprr[i] * bottom_cloud.ol$rrprrr[i]
    
    rrrrpp_path[i] <- bottom_cloud.ol$r[i] * bottom_cloud.ol$rr[i] * bottom_cloud.ol$rrr[i] * bottom_cloud.ol$rrrr[i] * bottom_cloud.ol$rrrrp[i] * bottom_cloud.ol$rrrrpp[i]
    rrrrpr_path[i] <- bottom_cloud.ol$r[i] * bottom_cloud.ol$rr[i] * bottom_cloud.ol$rrr[i] * bottom_cloud.ol$rrrr[i] * bottom_cloud.ol$rrrrp[i] * bottom_cloud.ol$rrrrpr[i]
    
    rrrrrp_path[i] <- bottom_cloud.ol$r[i] * bottom_cloud.ol$rr[i] * bottom_cloud.ol$rrr[i] * bottom_cloud.ol$rrrr[i] * bottom_cloud.ol$rrrrr[i] * bottom_cloud.ol$rrrrrp[i]
    rrrrrr_path[i] <- bottom_cloud.ol$r[i] * bottom_cloud.ol$rr[i] * bottom_cloud.ol$rrr[i] * bottom_cloud.ol$rrrr[i] * bottom_cloud.ol$rrrrr[i] * bottom_cloud.ol$rrrrrr[i]
    
  }
  
  
  bottom_cloud <- cbind(rrrrrr_path, rrrrrp_path, rrrrpr_path, rrrrpp_path, #concatenate all of the combinations together
                        rrrprr_path, rrrprp_path, rrrppr_path, rrrppp_path,
                        rrprrr_path, rrprrp_path, rrprpr_path, rrprpp_path,
                        rrpprr_path, rrpprp_path, rrpppr_path, rrpppp_path,
                        rprrrr_path, rprrrp_path, rprrpr_path, rprrpp_path,
                        rprprr_path, rprprp_path, rprppr_path, rprppp_path,
                        rpprrr_path, rpprrp_path, rpprpr_path, rpprpp_path,
                        rppprr_path, rppprp_path, rppppr_path, rppppp_path,
                        prrrrr_path, prrrrp_path, prrrpr_path, prrrpp_path,
                        prrprr_path, prrprp_path, prrppr_path, prrppp_path,
                        prprrr_path, prprrp_path, prprpr_path, prprpp_path,
                        prpprr_path, prpprp_path, prpppr_path, prpppp_path,
                        pprrrr_path, pprrrp_path, pprrpr_path, pprrpp_path,
                        pprprr_path, pprprp_path, pprppr_path, pprppp_path,
                        ppprrr_path, ppprrp_path, ppprpr_path, ppprpp_path,
                        pppprr_path, pppprp_path, pppppr_path, pppppp_path)
  rownames(bottom_cloud) <- rownames(bottom_cloud.ol) #set each gene or isoform name together

  bottom_cloud.ol <- bottom_cloud
  
  rownames(bottom_cloud.ol) <- rownames(bottom_cloud.ol)

  print("c")
  
  test.mat <- bottom_cloud.ol 

  bottomcloud_array.ol[j,,] <- matrix(unlist(test.mat), ncol=ncol(test.mat)) # add the run to the list of matricies

  print(dim(bottom_cloud.ol))

  print(j)
}

dimnames(bottomcloud_array.ol)[[2]] <- rownames(bottom_cloud.ol) # each gene/isoform
dimnames(bottomcloud_array.ol)[[3]] <- colnames(bottom_cloud.ol) # each combination
dimnames(bottomcloud_array.ol)[[1]] <- paste0("subsample.", as.character(1:100)) # each run 


datalist <- list() #this will store all of the information for the silhouette plots for each run
whichList <- c() #this will suggest the correct number of clusters for each run
for(i in 1:40){
  test <- fviz_nbclust(bottomcloud_array.ol[i,,], hcut, method = "silhouette", hc_method = "complete", k.max = 30)
  ### store the data from the silhouette plot 
  datalist[[i]] <- test$data #store silhouette lengths for clusters
  print(i)
  whichList[i] <- which.max(test$data[,2]) # retrieve optimal cluster number
}

numClust <- ceiling(mean(whichList)) #get the average optimal cluster number and ceiling it

AssMatThree.ol <- array(rep(0, nrow(AssMatThree.ol)*100), dim = c(ncol(AssMatThree.ol), 100))
for(i in 1:100){ 
  AssMatThree.ol[,i] <- cutree(hclust(dist(bottomcloud_array.ol[i,,])), k = numClust) # for each run, look at the clustering assignment for each gene
  print(i)
}
rownames(AssMatThree.ol) <- rownames(bottom_cloud.ol) # each row is a gene 
colnames(AssMatThree.ol) <- paste0("subsample.", as.character(1:100)) # each column is a subsample 

simBinary <- list()
for(i in 1:100){ # for each run 
  subsamp <- AssMatThree.ol[,i] # get the run under investigation
  runMat <- array(rep(0, nrow(AssMatThree.ol)*nrow(AssMatThree.ol)), dim = c(nrow(AssMatThree.ol),nrow(AssMatThree.ol)))
  # build a n_isoform x n_isoform matrix that will be stored with whether or not two genes cluster together in the same run
  for(j in 1:nrow(AssMatThree.ol)){ #for the jth gene
    for(k in 1:nrow(AssMatThree.ol)){#for the kth gene 
      if(subsamp[j] == subsamp[k]){ #if gene i and gene k are clustered the same
        runMat[j,k] <- 1 # denote that 
      } else { #otherwise 
        runMat[j, k] <- 0 #denote that they're not together 
      }
    }
    
  }
  simBinary[[i]] <- runMat #save the binary matrix for that run in a list
  print(i)
}

meansFromList.ol <- matrix(rep(0, nrow(AssMatThree.ol)*nrow(AssMatThree.ol)), nrow = nrow(AssMatThree.ol), ncol = nrow(AssMatThree.ol))
# avobe will make a n_isoforms by n_isoforms matrix that gets the frequency two genes clustered the same across all the runs
for(i in 1:nrow(meansFromList)){ #for the ith gene
  for(j in 1:nrow(meansFromList)){ #for the kth gene 
    binSum <- c() # this will be a vector of 1's and 0's, denoting the frequency gene i and gene j were clustered together
    for(k in 1:100) { #for each run 
      binSum[k] <- simBinary[[k]][i,j] #the kth element of binSum == whether or not gene i and gene j were clustered together in the kth run 
    }
	meansFromList.ol[i,j] <- mean(binSum) # the i,j spot is saved as the frequency gene/isoform i and j were together 
  }
  print(i)
}

save("meansFromList.ol", file = paste0(cancer_name, "-means-matrix-clustered.RData"))
#save the matrix of means 