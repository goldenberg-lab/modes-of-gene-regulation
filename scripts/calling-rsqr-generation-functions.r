library("optparse")
library("stringr") 
library("foreach") 
library("glmnet") 
library("parallel")
library("doParallel")
library("abind")

# select the cancer type from the shell script
option_list = list(
  make_option(c("-c", "--cancer"), type = "character", default="BRCA",
              help="cancer set abbreviation name", metavar = "character"),
  make_option(c("-k", "--sample"), type = "integer", default=1,
              help="set which subsampe you are making", metavar = "integer"),
  make_option(c("-i", "--indir"), type = "character", default="None",
              help="The directory to your independent variable", metavar = "character"),
  make_option(c("-p", "--parallel"), type = "character", default="FALSE",
              help="determines whether parallel (T) or genome-wide (F) is used.", metavar = "character"),
  make_option(c("-a", "--alpha"), type = "integer", default=1,
              help="Determines the alpha mixing parameter for the elastic net regression", metavar = "integer"),  
  make_option(c("-t", "--isTF"), type = "character", default="two",
              help="the level TFs are used at", metavar = "character"),  
  make_option(c("-l", "--genome-wide"), type = "character", default = "FALSE",
              help="Determines whether the genome-wide or singe-gene levels are used", metavar = "character"),
  make_option(c("-nc", "--cores"), type = "integer", default = 1,
              help ="Informs the number of cores required for the function to run", metavar = "integer"),
  make_option(c("-s", "--step"), type = "character", default = "one",
              help ="Informs the step of the regression you're on", metavar = "character"),
  make_option(c("-dir", "--homedir"), type = "character", default="/hpf/largeprojects/agoldenb/dustin/Data/TCGA/2016_01_28/",
              help="cancer set abbreviation name", metavar = "character")
);

opt_parser = OptionParser(option_list=option_list);

print("tothe")
opt = parse_args(opt_parser);

for(i in 1:length(opt)) {
  print(opt[[i]]) # returns what all of the options are in the output including default options
}
print(length(opt))

cancer_name = opt[[1]] # save the cancer name from the shell script
k = opt[[2]] # run number from shell script
directory <- opt[[3]] # the directory where your multiple file variables would be
isPar <- opt[[4]] # whether or not the program is run in parallel
isPar <- as.logical(isPar)
alpha <- opt[[5]] # the alpha mixing parameter for the ES regression
isTF <- opt[[6]] # whether or not the level is the transcription factor level
isGen <- opt[[7]] # whethe genome wide is used
isGen <- as.logical(isGen)
n.cores <- opt[[8]] # the number of cores required to run in parallel
stepNum <- opt[[9]] # which step in the tree plot is being looked at
homeDir <- opt[[10]] # which directory your cancer folders are in


if(stepNum == isTF){ # if you're at the TF level
  isTF <- TRUE # act like there are TFs
} else {
  isTF <- FALSE # act like there aren't TFs
}

print(isTF)

if(isGen == T){
  # make the directory the genome wide level, this will store the generating-rsqr-values scripts and data
  homeDir <- paste0(homeDir, cancer_name, "genome-wide/") 
} else {
  # make the directory the single gene level, this will store the generating-rsqr-values scripts and data
  
  homeDir <- paste0(homeDir, cancer_name, "singe-gene/")
}

set.seed(k) # set your seed according to the sampling scheme
print(k)

# Of course the code in here (and all of the other loaded data) will need to be changed
# In general, you may want to find a nice way to load data in that way people wont have to change the code when using this later
########
###############
###################
load(paste0(cancer_name, "cnv-rna-methyl-mirna-rem.RData")) # load the methylation data
load("HomoSapien-TFs.RData") # load in the H.sapien TF names
#####################
###############
########


sixtyp <- floor(ncol(rna.ol.all)*0.6) # get 60% of the individuals
cols <- sample(1:ncol(clinicalsame.ol), sixtyp) # get the collumn names for 60% of the individuals
######################################
######################################
#########
### THIS WILL NEED TO CHANGE ONCE THE FILES HAVE BEEN ADJUSTED
#########
######################################
######################################
if(isGen == T) { #if you're looking at the genome wide level
  load(paste0(cancer_name,"-pruned-snps.RData")) # load pruned snps, may also need to be adjusted
  rna.ol.all <- rna.ol.all[,cols] # get the 60% sample for mRNA expression
  cnv.ol.all <- cnv.ol.all[,cols] # get the 60% sample for cnvs
  clinicalsame.ol <- clinicalsame.ol[,cols] # get the 60% sample for clinica data
  mirna <- mirna.ol[,cols] # get the 60% sample for miRNAs...
  fullmRNA <- rna.ol.all
  miRNAgene <- t(mirna) # load the data so individuals are the collumns
  TF.indep <- t(fullmRNA[row.names(fullmRNA) %in% TF.symbol,]) # get the mRNA expression for the genes that transcribe transcription factors
  ############# This'll be need to be changed in particular 
  methyl_5utr.ol.all <- methyl_5utr.ol.all[,cols]
  methyl_3utr.ol.all <- methyl_3utr.ol.all[,cols]
  methyl_body.ol.all <- methyl_body.ol.all[,cols]
  methyl_tss1500.ol.all<- methyl_tss1500.ol.all[,cols]
  futr <- t(as.matrix(methyl_5utr.ol.all)) # load the 5' UTR methylation for the gene
  tutr <- t(as.matrix(methyl_3utr.ol.all)) # load the 3' UTR methylation for the gene
  body <- t(as.matrix(methyl_body.ol.all)) # load the body methylation for the gene
  tss <- t(as.matrix(methyl_tss1500.ol.all)) # load the tss methylation for the gene
  methyl <- abind(futr, tss, body, tutr) #bind the methyl groups together, will still need to be canged]
  ############
  snpLevel <- snpLevel.ol[cols,]
  CNV <- t(as.matrix(cnv.ol.all))
} else {
  
  cnv.ol.all <- cnv.ol.all[,cols]
  clinicalsame.ol <- clinicalsame.ol[,cols]
  # load the TF 3d matrix
  # load the miRNA 3d matrix
  # load the first snp file (or 3d snp matrix)
  # load the first methylation file (or 3d methylation matrix)
  CNV <- t(as.matrix(cnv.ol.all))
  
}

source("generating-rsqr-values.r") # call the three functions requried to make the path vectors
print(dim(rna.ol.all))
if(isPar == T) { # if you're running the pipeline in parallel
  chunk <- nrow(rna.ol.all) %/% (n.cores-1) # get the number of genes going into each core
  n <- nrow(rna.ol.all) # get the number of cores 
  r  <- rep(1:ceiling(n/chunk),each=chunk)[1:n] # generate vector for each pipeline
  d <- split(rna.ol.all,r) # split the mRNA expression into "chunk" cuts with "r" rows
  nameD <- split(rownames(rna.ol.all), r) #get the gene/isoform names into each split
  cuts <- cut(1:nrow(rna.ol.all), n.cores) 
  cl <- makeCluster(n.cores) # get a cluster for each core 
  registerDoParallel(cl) #get them ready
  
  thirdD.ol <- foreach(i = 1:n.cores, .combine = 'rbind', .packages = c("glmnet", "abind")) %dopar% treeGenFullFunct(d[i], rna_names = nameD[i], one = miRNAgene, two = TF.indep, three = methyl, four = CNV, five = snpLevel)
  # for each of the 30 chunks and generate the R^2 values for the level that you're looking at 
  
  saveParallels(thirdD.ol) # combines each of the chunks back together again to be able to get all of the R^2 values in the right spot
  print("Leve is saved in your directory :)")
  stopCluster(cl) # get rid of the clusters
  
} else {
  treeGenFullFunct(rna.ol.all, rna_names = rownames(rna.ol.all), one = miRNAgene, two = TF.indep, three = methyl, four = CNV, five = snpLevel)
  # get all of the R^2 value for the level that you're looking at as well as all of the path vectors. Look at the generating-rsqr-values.r 
  # files for more details about this function
  print("Level saved in your directory :)")
}
