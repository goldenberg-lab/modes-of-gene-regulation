get_rownamed_funct <- function(x, y){ 
  test <- matrix(unlist(x), ncol = ncol(cnv.ol.all))
  rownames(test) <- unlist(y)
  return(test)
}

glmStep <- function(independent, dependent) {
	# independent = the mode of gene regulation that you are using
		# format: Four options are allowed
			# 1) 3d Array:
				# a) dim 1 = genes
				# b) dim 2 = individuals
				# 3) dim 3 = variables
			# 2) multiple files: this is similar to the 3d array, just each gene is now it's own file
				# this is good if you'd have a very large 3d array and it would be more convenient to load 
				# one gene at a time. Here, each gene is saved as a different file in a directory that 
				# you specify as an optional argument. Each gene will have the rows as individuals and collumns are
				# variables attached to the gene of interest
			# 3) single file (genome wide level): this is a single speadsheet where the rows are all of the individuals
				# and the collumns are all of the values for all of the variables in your regulatory modality
			# 4) single file (single gene level): this is a single spreadsheet where the rows are genes and the collumns are individuals
				# it is used when there is only one variable modelled for each gene at the regulatory modality 
	# dependent = matrix or list with either your RNA, Predicted values, or Residual values
		# format: rows are the genes and the columns are the individuals
	# Returns: A list containing two matricies and one vector (all unnamed)
		# first element = $predicted. These are the predicted values for the model
		# second element = $residual. These are the residual values from the model
		# third element = $Rsqr. This is the R^2 value of the regression model
  
	store <- independent
	if(typeof(dependent) == "list") {
		dependent <- matrix(unlist(dependent), ncol = ncol(rna.ol.all))
		rownames(dependent) <- rownames(rna.ol.all)[1:nrow(dependent)]
		# Aain, you want to have both the independent and dependent variables as matricies (dependent can be a vector)
		# where the number of columns of the dependent variable (or the length) is equivalent to the number of rows of
		# the independent variable 
	}
	if(length(dim(independent)) == 3) { 
		pp <- c()
		predicted <- c()
		residual <- c()
		for(i in 1:nrow(dependent)){ # For all the genes
			indepGene <- independent[rownames(independent) == rownames(dependent[i,]),,] # make sure that the data for the independent variable is the same isoform/gene as the dependent variable selected
			if(isTF == T){
					indepGene <- indepGene[,colnames(indepGene) != rownames(dependent)[i]] # remove the values of your dependent variable from your independent variable
				}
			if(ncol(indepGene) + 104 <= nrow(indepGene)) { # If you pass the heuristic
				if(var(dependent[i,]) != 0 & max(range(var(indepGene))) != 0){ # if both dependent and independent variables have variance
					lm.indep <- as.data.frame(indepGene) # make the independent variable a data frame with the gene of interest
					linmod <- lm(dependent[i,] ~., lm.indep) # multiple linear regression with each variable as a cofactor
					predicted <- rbind(predicted, predict(linmod)) # find y_hat 
					residual <-  rbind(residual, rna.ol.all[i,] - predicted[i,]) # residual = y (expression) - y_hat
					pp[i] <- unlist(summary(linmod)[8]) # get the R^2 value
				} else {
					predicted <- rbind(predicted, rep(0, ncol(rna.ol.all))) # make the predicted value = 0
					residual <-  rbind(residual,rep(0, ncol(rna.ol.all))) # residual values will be the dependent variable, but the whole path will = 0 so to save running time we make this 0
					pp[i] <- 0 # generate the R^2 value for the gene
					
				}            
			} else { #if the heuristic is not passed
				fit.pred <- try(cv.glmnet(indepGene, dependent[i,], alpha = alpha), silent = T) # apply k-fold cross validation to optimize lambda value
				if(var(dependent[i,]) != 0 & max(range(var(indepGene))) != 0 & class(fit.pred) != "try-error") { # if both variables have variance and cross validation worked
					bestlab.pred <- fit.pred$lambda.min # minimize lambda value 
					test.pred <- glmnet(indepGene, dependent[i,], alpha = alpha, lambda = bestlab.pred) # run an ES regression with the optimal lambda value from predicted
					predicted <- rbind(predicted, t(predict(test.pred, newx=indepGene, s = bestlab.pred))) # find y_hat based on ES coefficients
					residual <- rbind(residual, rna.ol.all[i,] - predicted[i,]) # residual = y - y_hat
					pp[i] <- unlist(test.pred$dev.ratio) # generate the R^2 value for the gene
				} else {
					predicted <- rbind(predicted, rep(0, ncol(rna.ol.all))) # make the predicted value = 0
					residual <-  rbind(residual,rep(0, ncol(rna.ol.all))) # residual values will be the dependent variable, but the whole path will = 0 so to save running time we make this 0
					pp[i] <- 0 # generate the R^2 value for the gene
				}
			}
		}
	}
	if(length(table(!(rownames(independent) %in% rownames(dependent))))==1) { # If option 2 or 3 is chosen for independent variables
		pp <- c()
		predicted <- c()
		residual <- c()    
		for(i in 1:nrow(dependent)){ # apply the previous bit to every gene being investigated
			if(directory != "None") { # if we're looking at multiple files (optio 2)
				setwd(directory) # choose the multiplle file directory
				independent <- load(paste0(cancer_name,"-",rownames(rna.ol.all)[i], ".RData")) # load the file for your gene
				independent <- independent[cols,] # make sure that the same individuals are in there
				setwd(homeDir) # go back to the original diectory
			}
			if(isTF == T){ # If the level is a transcription factor
				if(directory == "None") { # if we're not using multiple files
					independent <- store # reset the independent variable
				}
				independent <- independent[,colnames(independent) != rownames(dependent)[i]] # remove the variable identical to the dependent vairables from the dataset
				#
				##
				####
				##### When we're at the isoform level, you'll probably have to remoove every isoform where the collumn of the independent variable == the depent variable
				##### Maybe look for a naming convention with the isoform where you an remove these rows more easily? Since I don't have the isoform names this'll be up to you
				###
				##
				#
			}
			if(ncol(independent) + 104 <= nrow(independent)) { # If the heuristic is passed  
				if(var(dependent[i,]) != 0 & max(range(var(independent))) != 0){ #and dependent/independent have no variance
					lm.indep <- as.data.frame(independent) # make the independent variable a dataframe for the linear model
					linmod <- lm(dependent[i,] ~., lm.indep) # multiple linear regression with each variable relating from gene
					predicted <- rbind(predicted, predict(linmod))# find y_hat 
					residual <-  rbind(residual, rna.ol.all[i,] - predicted[i,]) # residual = y - y_hat,
					pp[i] <- unlist(summary(linmod)[8]) # get the R^2 value
				} else {
					predicted <- rbind(predicted, rep(0, ncol(rna.ol.all))) # make the predicted value = 0
					residual <-  rbind(residual,rep(0, ncol(rna.ol.all))) # residual values will be the dependent variable, but the whole path will = 0 so to save running time we make this 0
					pp[i] <- 0 # generate the R^2 value for the gene
				}
			} else{
				fit.pred <- try(cv.glmnet(independent, dependent[i,], alpha = alpha), silent = T) # run k-fold cross validation 
				if(var(dependent[i,]) != 0 & max(range(var(independent))) != 0 & class(fit.pred) != "try-error") { #if the variance is not 0 and cross validation worked
					bestlab.pred <- fit.pred$lambda.min # find the minimized lambda value
					test.pred <- glmnet(independent, dependent[i,], alpha = alpha, lambda = bestlab.pred) # run an ES regression with the optimal lambda value from predicted
					predicted <- rbind(predicted, t(predict(test.pred, newx=independent, s = bestlab.pred))) # find y_hat
					residual <- rbind(residual, rna.ol.all[i,] - predicted[i,]) #residual = y - y_hat
					pp[i] <- unlist(test.pred$dev.ratio) # generate the R^2 value for the gene
				} else {
					predicted <- rbind(predicted, rep(0, ncol(rna.ol.all))) # make the predicted value = 0
					residual <-  rbind(residual,rep(0, ncol(rna.ol.all))) # residual values will be the dependent variable, but the whole path will = 0 so to save running time we make this 0
					pp[i] <- 0 # generate the R^2 value for the gene
				}
		  }
		}
	}
	else{ #if there is only one vairable being modelled (regular linear regression)
		pp <- c()
		predicted <- c()
		residual <- c()    
		for(i in 1:nrow(dependent)) { # repeat the same process for all the genes being investigated
			indepSingGene <- independent[rownames(independent) == rownames(dependent[i,]),] # assign the proper gene for the independent variable
			if(var(dependent[i,]) != 0 & max(range(var(indepSingGene))) != 0){ #both independent/dependent has variance
				linmod <- lm(dependent[i,] ~ indepSingGene) # do a regular linear regression
				predicted <- rbind(predicted, predict(linmod)) # find y_hat 
				residual <-  rbind(residual, rna.ol.all[i,] - predicted[i,]) # residual = y - y_hat
				pp[i] <- unlist(summary(linmod)[8]) # get the R^2 value
			} else {
				predicted <- rbind(predicted, rep(0, ncol(rna.ol.all))) # make the predicted value = 0
				residual <-  rbind(residual,rep(0, ncol(rna.ol.all))) # residual values will be the dependent variable, but the whole path will = 0 so to save running time we make this 0
				pp[i] <- 0 # generate the R^2 value for the gene
			}
		}
	}
  return(list(predicted = predicted, residual = residual, Rsqr = pp))
}


treeGenFullFunct <- function(rna.ol.all, rna_names, one = miRNAgene, two = TF.indep, three = methyl, four = CNV, five = snpLevel) {
	# Variables
	#rna.ol.all: a matrix (or list) containing the mRNA expression of all of your genes and all of the individuals you're looking at
	#rna_names: a vector containing all the gene names being observed in this function run
	#one-five: each level in the path the typical path shown below
	#one = miRNA, two = TF, three = methylation, four = CNV, five = eQTL
	# Returns
	# A matrix where the rows are different genes and the collumns are each R^2 value generated. 
	# these will be used later to make the path vectors
	if(typeof(rna.ol.all) == "list"){
		#converts the mRNA list into a matrix
		rna.ol.all <- get_rownamed_funct(rna.ol.all, rna_names)
	}
	rownames(rna.ol.all) <- rownames(rna.ol.all)
	print("start!")
	if(stepNum == "one") { #If we're looking at the first step
		mRNA <- unlist(rna.ol.all[1,]) # dependent variable for clinical level
		
		age <- clinicalsame.ol[1,] # age of diagnosis
		tis <- clinicalsame.ol[2,] # tissue of tumour
		gen <- clinicalsame.ol[3,] # gender of individuals
		sT <- clinicalsame.ol[4,] # tumour size
		sM <- clinicalsame.ol[5,] # if metasticis/size of secondary tumour
		sN <- clinicalsame.ol[6,] # size in lymph nodes/if able to make it
		print("A")
		mRNA.env <- lm(mRNA ~ age + tis + gen + sT + sM + sN)
		# multiple linear regression to look at role of clinical variables
		# on mRNA expression
		
		mRNA.env.resi <- residuals(mRNA.env) # residual from expression
		mRNA.env.pred <- predict(mRNA.env) # predicted values from expression
		p <- c()
		p[1] <- unlist(summary(mRNA.env)[8]) # get the R^2 value or that regression
		for(i in 2:nrow(rna.ol.all)){ #repeat the above process for every gene being observed
			mRNA <- unlist(rna.ol.all[i,])
			mRNA.env <- lm(mRNA ~ tis + gen + sT + sM + sN)
			mRNA.env.resi <- rbind(mRNA.env.resi, residuals(mRNA.env))
			mRNA.env.pred <- rbind(mRNA.env.pred, predict(mRNA.env))
			p[i] <- unlist(summary(mRNA.env)[8])
		}
		print("b")
		row.names(mRNA.env.resi) <- rownames(rna.ol.all)
		row.names(mRNA.env.pred) <- rownames(rna.ol.all)
		r <- 1-p # get the % variation that is not accounted for in this model
		env.ol <- cbind(r, p) # combine the R^2 values for the mRNA~Clinical level
		row.names(env.ol) <- row.names(rna.ol.all) # assign the gene names to those R^2 levels
		### Env to miRNA level  
		
		
		# After this, the function gets very very repetitive
		# The glmStep function is applied where the predicted or residual
		# values of the previous regression are the dependent variables
		# the level that you've specified in the parameters of the function are 
		# your independent variables. This process is repeated for the entire tree.
		miRNA_p <- glmStep(independent = one, mRNA.env.pred) #get y_hat, residuals (y-y_hat), and R^2 values for each gene and save it at each level
		env_miRNA_pp.ol <- miRNA_p$predicted
		env_miRNA_pr.ol <- miRNA_p$residual
		pp<- miRNA_p$Rsqr
		miRNA_r <- glmStep(independent = one, mRNA.env.resi)
		env_miRNA_rp.ol <- miRNA_r$predicted
		env_miRNA_rr.ol <- miRNA_r$residual
		rp <- miRNA_r$Rsqr
		print("c")
		row.names(env_miRNA_pp.ol) <- rownames(rna.ol.all) # get the rownames for the predicted/residual values
		row.names(env_miRNA_pr.ol) <- rownames(rna.ol.all)
		row.names(env_miRNA_rp.ol) <- rownames(rna.ol.all)
		row.names(env_miRNA_rr.ol) <- rownames(rna.ol.all)
		pr <- 1 - pp
		print("d")
		rr <- 1-rp
		miRNA.ol <- cbind(rr, rp, pr, pp) # get all the R^2 values for the second level
		row.names(miRNA.ol) <- row.names(rna.ol.all)
		saveLvlOne.ol <- c("env_miRNA_pp.ol", "env_miRNA_pr.ol", "env_miRNA_rp.ol", "env_miRNA_rr.ol")
		# this is a vector with the names of all of the files that need to be saved in the RData file format.
		if(isPar==T){ #if this is being run in prallel
			###############
			return(list(env_miRNA_rr.ol, env_miRNA_rp.ol, env_miRNA_pr.ol, env_miRNA_pp.ol, miRNA.ol, env.ol))
			#############
		} else {
			save(list = c("env.ol", "miRNA.ol"), file = paste0(cancer_name, "-one-rsqrs-", as.character(k), ".RData"))
			save(list = saveLvlOne.ol, file = paste0(cancer_name, "-one-", as.character(k), ".RData"))
			return(env.ol)
		}
	}
	if(stepNum == "two"){ #If you're looking at the second level (likely the transcription factor level)
		load(paste0(cancer_name, "-one-", as.character(k), ".RData")) # load the y_hat and residual values from the previous level
		if(isPar == T){
			env_miRNA_pp.ol <- env_miRNA_pp.ol[rownames(env_miRNA_pp.ol) %in% rownames(rna.ol.all),] 
			# since each level is done in parallel and they are combined at the end you need to select the slice that you need 
			env_miRNA_pr.ol <- env_miRNA_pr.ol[rownames(env_miRNA_pr.ol) %in% rownames(rna.ol.all),]
			env_miRNA_rp.ol <- env_miRNA_rp.ol[rownames(env_miRNA_rp.ol) %in% rownames(rna.ol.all),]
			env_miRNA_rr.ol <- env_miRNA_rr.ol[rownames(env_miRNA_rr.ol) %in% rownames(rna.ol.all),]
		}
		TF_pp <- glmStep(two, env_miRNA_pp.ol)
		TF_ppp <- TF_pp$predicted
		TF_ppr <- TF_pp$residual
		ppp <- TF_pp$Rsqr
		
		TF_pr <- glmStep(two, env_miRNA_pr.ol)
		TF_prp <- TF_pr$predicted
		TF_prr <- TF_pr$residual
		prp <- TF_pr$Rsqr
		
		TF_rp <- glmStep(two, env_miRNA_rp.ol)
		TF_rpp <- TF_rp$predicted
		TF_rpr <- TF_rp$residual
		rpp <- TF_rp$Rsqr
		
		TF_rr <- glmStep(two, env_miRNA_rr.ol)
		TF_rrp <- TF_rr$predicted
		TF_rrr <- TF_rr$residual
		rrp <- TF_rr$Rsqr
		print("e")
		row.names(TF_prp) <- row.names(rna.ol.all) # match rownames of matrix with rownames from previous step
		row.names(TF_prr) <- row.names(rna.ol.all)
		row.names(TF_ppp) <- row.names(rna.ol.all)
		row.names(TF_ppr) <- row.names(rna.ol.all)
		row.names(TF_rrp) <- row.names(rna.ol.all)
		row.names(TF_rrr) <- row.names(rna.ol.all)
		row.names(TF_rpp) <- row.names(rna.ol.all)
		row.names(TF_rpr) <- row.names(rna.ol.all)

		rpr <- 1-rpp
		rrr <- 1-rrp
		ppr <- 1-ppp
		prr <- 1-prp
		print("f")
		# this is done to insure that the format of the output 
		# for the previous level is in proper shape for the next 
		# level
		TF_ppp.ol <- data.matrix(TF_ppp)
		TF_ppr.ol <- data.matrix(TF_ppr)
		TF_prp.ol <- data.matrix(TF_prp)
		TF_prr.ol <- data.matrix(TF_prr)
		TF_rpp.ol <- data.matrix(TF_rpp)
		TF_rpr.ol <- data.matrix(TF_rpr)
		TF_rrp.ol <- data.matrix(TF_rrp)
		TF_rrr.ol <- data.matrix(TF_rrr)

		TF.ol <- cbind(rrr, rrp, rpr, rpp, prr, prp, ppr, ppp) #combine all of the R^2 files at the third level of the path
		row.names(TF.ol) <- row.names(rna.ol.all)

		TF_files.ol <- c("TF_ppr.ol", "TF_prp.ol", "TF_prr.ol", "TF_ppp.ol", 
						"TF_rpp.ol", "TF_rpr.ol", "TF_rrp.ol", "TF_rrr.ol")
		if(isPar == T){
		
		return(list(TF_rrr.ol, TF_rrp.ol, TF_rpr.ol, TF_rpp.ol, TF_prr.ol, TF_prp.ol, TF_ppr.ol, TF_ppp.ol, TF.ol))
		} else {
			save(list = TF_files.ol, file = paste0(cancer_name, "-two-", as.character(k), ".RData"))
			save(list = "TF.ol",  file = paste0(cancer_name, "-two-rsqrs-", as.character(k), ".RData"))
			return(TF.ol)
		}
	}
	if(stepNum == "three") {
		load(paste0(cancer_name, "-two-", as.character(k), ".RData")) # save the predicted and residual levels from the prior step
		
		if(isPar == T){
			TF_ppp.ol <- TF_ppp.ol[rownames(TF_ppp.ol) %in% rownames(rna.ol.all),]
			TF_ppr.ol <- TF_ppr.ol[rownames(TF_ppr.ol) %in% rownames(rna.ol.all),]
			TF_prp.ol <- TF_prp.ol[rownames(TF_prp.ol) %in% rownames(rna.ol.all),]
			TF_prr.ol <- TF_prr.ol[rownames(TF_prr.ol) %in% rownames(rna.ol.all),]
			TF_rpp.ol <- TF_rpp.ol[rownames(TF_rpp.ol) %in% rownames(rna.ol.all),]
			TF_rpr.ol <- TF_rpr.ol[rownames(TF_rpr.ol) %in% rownames(rna.ol.all),]
			TF_rrp.ol <- TF_rrp.ol[rownames(TF_rrp.ol) %in% rownames(rna.ol.all),]
			TF_rrr.ol <- TF_rrr.ol[rownames(TF_rrr.ol) %in% rownames(rna.ol.all),]
		}
		
		meth_ppp <- glmStep(three, TF_ppp.ol)
		TF_meth_pppp <- meth_ppp$predicted
		TF_meth_pppr <- meth_ppp$residual
		pppp <- meth_ppp$Rsqr
		print("is")
		meth_ppr <- glmStep(three, TF_ppr.ol)
		TF_meth_pprp <- meth_ppr$predicted
		TF_meth_pprr <- meth_ppr$residual
		pprp <- meth_ppr$Rsqr
		print("this")
		meth_prp <- glmStep(three, TF_prp.ol)
		TF_meth_prpp <- meth_prp$predicted
		TF_meth_prpr <- meth_prp$residual
		prpp <- meth_prp$Rsqr
		print("probz")
		meth_prr <- glmStep(three, TF_prr.ol)
		TF_meth_prrp <- meth_prr$predicted
		TF_meth_prrr <- meth_prr$residual
		prrp <- meth_prr$Rsqr
	  
		meth_rpp <- glmStep(three, TF_rpp.ol)
		TF_meth_rppp <- meth_rpp$predicted
		TF_meth_rppr <- meth_rpp$residual
		rppp <- meth_rpp$Rsqr
	  
		meth_rpr <- glmStep(three, TF_rpr.ol)
		TF_meth_rprp <- meth_rpr$predicted
		TF_meth_rprr <- meth_rpr$residual
		rprp <- meth_rpr$Rsqr
	  
		meth_rrp <- glmStep(three, TF_rrp.ol)
		TF_meth_rrpp <- meth_rrp$predicted
		TF_meth_rrpr <- meth_rrp$residual
		rrpp <- meth_rrp$Rsqr
	  
		meth_rrr <- glmStep(three, TF_rrr.ol)
		TF_meth_rrrp <- meth_rrr$predicted
		TF_meth_rrrr <- meth_rrr$residual
		rrrp <- meth_rrr$Rsqr
	  

	  
		rrrr <- 1 - rrrp
		print("h")
		rrpr <- 1 - rrpp
		rprr <- 1 - rprp
		rppr <- 1 - rppp
		prrr <- 1 - prrp
		prpr <- 1 - prpp
		pprr <- 1 - pprp
		pppr <- 1 - pppp
		print("i")
		rownames(TF_meth_pppp) <- rownames(rna.ol.all)
		rownames(TF_meth_pppr) <- rownames(rna.ol.all)
		rownames(TF_meth_pprp) <- rownames(rna.ol.all)
		rownames(TF_meth_pprr) <- rownames(rna.ol.all)
		
		rownames(TF_meth_prpp) <- rownames(rna.ol.all)
		rownames(TF_meth_prpr) <- rownames(rna.ol.all)
		rownames(TF_meth_prrp) <- rownames(rna.ol.all)
		rownames(TF_meth_prrr) <- rownames(rna.ol.all)
		
		rownames(TF_meth_rppp) <- rownames(rna.ol.all)
		rownames(TF_meth_rppr) <- rownames(rna.ol.all)
		rownames(TF_meth_rprp) <- rownames(rna.ol.all)
		rownames(TF_meth_rprr) <- rownames(rna.ol.all)
		
		rownames(TF_meth_rrpp) <- rownames(rna.ol.all)
		rownames(TF_meth_rrpr) <- rownames(rna.ol.all)
		rownames(TF_meth_rrrp) <- rownames(rna.ol.all)
		rownames(TF_meth_rrrr) <- rownames(rna.ol.all)
		
		meth_pppp.ol <- data.matrix(TF_meth_pppp)
		meth_pppr.ol <- data.matrix(TF_meth_pppr)
		meth_pprp.ol <- data.matrix(TF_meth_pprp)
		meth_pprr.ol <- data.matrix(TF_meth_pprr)
		
		meth_prpp.ol <- data.matrix(TF_meth_prpp)
		meth_prpr.ol <- data.matrix(TF_meth_prpr)
		meth_prrp.ol <- data.matrix(TF_meth_prrp)
		meth_prrr.ol <- data.matrix(TF_meth_prrr)
		
		meth_rppp.ol <- data.matrix(TF_meth_rppp)
		meth_rppr.ol <- data.matrix(TF_meth_rppr)
		meth_rprp.ol <- data.matrix(TF_meth_rprp)
		meth_rprr.ol <- data.matrix(TF_meth_rprr)
		
		meth_rrpp.ol <- data.matrix(TF_meth_rrpp)
		meth_rrpr.ol <- data.matrix(TF_meth_rrpr)
		meth_rrrp.ol <- data.matrix(TF_meth_rrrp)
		meth_rrrr.ol <- data.matrix(TF_meth_rrrr)
		
		
		#### Removing novar meth level here
		
		
		print(dim(meth_pppp.ol))
		print(dim(rrrr))
		methnames.ol <- c("meth_pppp.ol", "meth_pppr.ol","meth_pprp.ol", "meth_pprr.ol",
						"meth_prpp.ol", "meth_prpr.ol","meth_prrp.ol", "meth_prrr.ol",
						"meth_rppp.ol", "meth_rppr.ol","meth_rprp.ol", "meth_rprr.ol",
						"meth_rrpp.ol", "meth_rrpr.ol","meth_rrrp.ol", "meth_rrrr.ol")
		
		methyl.ol <- cbind(rrrr, rrrp, rrpr, rrpp,
							rprr, rprp, rppr, rppp,
							prrr, prrp, prpr, prpp,
							pprr, pprp, pppr, pppp)
		row.names(methyl.ol) <- row.names(rna.ol.all)
		if(isPar == T) {
			return(list(meth_rrrr.ol, meth_rrrp.ol, meth_rrpr.ol, meth_rrpp.ol,
						meth_rprr.ol, meth_rprp.ol, meth_rppr.ol, meth_rppp.ol,
						meth_prrr.ol, meth_prrp.ol, meth_prpr.ol, meth_prpp.ol,
						meth_pprr.ol, meth_pprp.ol, meth_pppr.ol, meth_pppp.ol, methyl.ol))
		} else {
			save(list = methnames.ol, file = paste0(cancer_name, "-three-", as.character(k), ".RData"))
			save(list = "methyl.ol",  file = paste0(cancer_name, "-three-rsqrs-", as.character(k), ".RData"))
			return(methyl.ol)
		}
	}
	if(stepNum == "four"){
		#### meth to CNV level
		load(paste0(cancer_name, "-three-", as.character(k), ".RData"))
		if(isPar == T){
			meth_pppp.ol <- meth_pppp.ol[rownames(meth_pppp.ol) %in% rownames(rna.ol.all),]
			meth_pppr.ol <- meth_pppr.ol[rownames(meth_pppr.ol) %in% rownames(rna.ol.all),]
			meth_pprp.ol <- meth_pprp.ol[rownames(meth_pprp.ol) %in% rownames(rna.ol.all),]
			meth_pprr.ol <- meth_pprr.ol[rownames(meth_pprr.ol) %in% rownames(rna.ol.all),]
			
			meth_prpp.ol <- meth_prpp.ol[rownames(meth_prpp.ol) %in% rownames(rna.ol.all),]
			meth_prpr.ol <- meth_prpr.ol[rownames(meth_prpr.ol) %in% rownames(rna.ol.all),]
			meth_prrp.ol <- meth_prrp.ol[rownames(meth_prrp.ol) %in% rownames(rna.ol.all),]
			meth_prrr.ol <- meth_prrr.ol[rownames(meth_prrr.ol) %in% rownames(rna.ol.all),]
			
			meth_rppp.ol <- meth_rppp.ol[rownames(meth_rppp.ol) %in% rownames(rna.ol.all),]
			meth_rppr.ol <- meth_rppr.ol[rownames(meth_rppr.ol) %in% rownames(rna.ol.all),]
			meth_rprp.ol <- meth_rprp.ol[rownames(meth_rprp.ol) %in% rownames(rna.ol.all),]
			meth_rprr.ol <- meth_rprr.ol[rownames(meth_rprr.ol) %in% rownames(rna.ol.all),]
			
			meth_rrpp.ol <- meth_rrpp.ol[rownames(meth_rrpp.ol) %in% rownames(rna.ol.all),]
			meth_rrpr.ol <- meth_rrpr.ol[rownames(meth_rrpr.ol) %in% rownames(rna.ol.all),]
			meth_rrrp.ol <- meth_rrrp.ol[rownames(meth_rrrp.ol) %in% rownames(rna.ol.all),]
			meth_rrrr.ol <- meth_rrrr.ol[rownames(meth_rrrr.ol) %in% rownames(rna.ol.all),]
		}
		
		cnv_pppp <- glmStep(four, meth_pppp.ol)
		meth_CNV_ppppp <- cnv_pppp$predicted
		meth_CNV_ppppr <- cnv_pppp$residual
		ppppp <- cnv_pppp$Rsqr
		
		cnv_pppr <- glmStep(four, meth_pppr.ol)
		meth_CNV_ppprp <- cnv_pppr$predicted
		meth_CNV_ppprr <- cnv_pppr$residual
		ppprp <- cnv_pppr$Rsqr
		
		cnv_pprp <- glmStep(four, meth_pprp.ol)
		meth_CNV_pprpp <- cnv_pprp$predicted
		meth_CNV_pprpr <- cnv_pprp$residual
		pprpp <- cnv_pprp$Rsqr
		
		cnv_pprr <- glmStep(four, meth_pprr.ol)
		meth_CNV_pprrp <- cnv_pprr$predicted
		meth_CNV_pprrr <- cnv_pprr$residual
		pprrp <- cnv_pprr$Rsqr
		##########################################
		cnv_prpp <- glmStep(four, meth_prpp.ol)
		meth_CNV_prppp <- cnv_prpp$predicted
		meth_CNV_prppr <- cnv_prpp$residual
		prppp <- cnv_prpp$Rsqr
		
		cnv_prpr <- glmStep(four, meth_prpr.ol)
		meth_CNV_prprp <- cnv_prpr$predicted
		meth_CNV_prprr <- cnv_prpr$residual
		prprp <- cnv_prpr$Rsqr
		
		cnv_prrp <- glmStep(four, meth_prrp.ol)
		meth_CNV_prrpp <- cnv_prrp$predicted
		meth_CNV_prrpr <- cnv_prrp$residual
		prrpp <- cnv_prrp$Rsqr
		
		cnv_prrr <- glmStep(four, meth_prrr.ol)
		meth_CNV_prrrp <- cnv_prrr$predicted
		meth_CNV_prrrr <- cnv_prrr$residual
		prrrp <- cnv_prrr$Rsqr
		############################################
		############################################
		cnv_rppp <- glmStep(four, meth_rppp.ol)
		meth_CNV_rpppp <- cnv_rppp$predicted
		meth_CNV_rpppr <- cnv_rppp$residual
		rpppp <- cnv_rppp$Rsqr
		
		cnv_rppr <- glmStep(four, meth_rppr.ol)
		meth_CNV_rpprp <- cnv_rppr$predicted
		meth_CNV_rpprr <- cnv_rppr$residual
		rpprp <- cnv_rppr$Rsqr
		
		cnv_rprp <- glmStep(four, meth_rprp.ol)
		meth_CNV_rprpp <- cnv_rprp$predicted
		meth_CNV_rprpr <- cnv_rprp$residual
		rprpp <- cnv_rprp$Rsqr
		
		cnv_rprr <- glmStep(four, meth_rprr.ol)
		meth_CNV_rprrp <- cnv_rprr$predicted
		meth_CNV_rprrr <- cnv_rprr$residual
		rprrp <- cnv_rprr$Rsqr
		##########################################
		cnv_rrpp <- glmStep(four, meth_rrpp.ol)
		meth_CNV_rrppp <- cnv_rrpp$predicted
		meth_CNV_rrppr <- cnv_rrpp$residual
		rrppp <- cnv_rrpp$Rsqr
		
		cnv_rrpr <- glmStep(four, meth_rrpr.ol)
		meth_CNV_rrprp <- cnv_rrpr$predicted
		meth_CNV_rrprr <- cnv_rrpr$residual
		rrprp <- cnv_rrpr$Rsqr
		
		cnv_rrrp <- glmStep(four, meth_rrrp.ol)
		meth_CNV_rrrpp <- cnv_rrrp$predicted
		meth_CNV_rrrpr <- cnv_rrrp$residual
		rrrpp <- cnv_rrrp$Rsqr
		
		cnv_rrrr <- glmStep(four, meth_rrrr.ol)
		meth_CNV_rrrrp <- cnv_rrrr$predicted
		meth_CNV_rrrrr <- cnv_rrrr$residual
		rrrrp <- cnv_rrrr$Rsqr

		ppppr <- 1 - ppppp
		print("j")
		ppprr <- 1 - ppprp
		pprpr <- 1 - pprpp
		pprrr <- 1 - pprrp
		prppr <- 1 - prppp
		prprr <- 1 - prprp
		prrpr <- 1 - prrpp
		prrrr <- 1 - prrrp
		rpppr <- 1 - rpppp
		rpprr <- 1 - rpprp
		rprpr <- 1 - rprpp
		rprrr <- 1 - rprrp
		rrppr <- 1 - rrppp
		rrprr <- 1 - rrprp
		rrrpr <- 1 - rrrpp
		rrrrr <- 1 - rrrrp
		print("k")

		rownames(meth_CNV_ppppp) <- rownames(rna.ol.all)
		print(dim(meth_CNV_ppppp))
		print(dim(rna.ol.all))
		
		rownames(meth_CNV_ppppr) <- rownames(rna.ol.all)
		print(dim(meth_CNV_ppppr))

		rownames(meth_CNV_ppprp) <- rownames(rna.ol.all)
		print(dim(meth_CNV_ppprp))
	  
		rownames(meth_CNV_ppprr) <- rownames(rna.ol.all)
		print(dim(meth_CNV_ppprr))
	  
		print("meth_CNV_pprpp")
		rownames(meth_CNV_pprpp) <- rownames(rna.ol.all)
	  
		print("meth_CNV_pprpr")
		print(dim(meth_CNV_pprpr))
		#return(meth_CNV_pprpr)
		rownames(meth_CNV_pprpr) <- rownames(rna.ol.all)
		
		print("meth_CNV_pprrp")
		rownames(meth_CNV_pprrp) <- rownames(rna.ol.all)
		rownames(meth_CNV_pprrr) <- rownames(rna.ol.all)
		rownames(meth_CNV_prppp) <- rownames(rna.ol.all)
		rownames(meth_CNV_prppr) <- rownames(rna.ol.all)
		rownames(meth_CNV_prprp) <- rownames(rna.ol.all)
		rownames(meth_CNV_prprr) <- rownames(rna.ol.all)
		rownames(meth_CNV_prrpp) <- rownames(rna.ol.all)
		rownames(meth_CNV_prrpr) <- rownames(rna.ol.all)
		print("prrrp")
		rownames(meth_CNV_prrrp) <- rownames(rna.ol.all)
		rownames(meth_CNV_prrrr) <- rownames(rna.ol.all)
		rownames(meth_CNV_rpppp) <- rownames(rna.ol.all)
		rownames(meth_CNV_rpppr) <- rownames(rna.ol.all)
		print("rpprp")
		rownames(meth_CNV_rpprp) <- rownames(rna.ol.all)
		rownames(meth_CNV_rpprr) <- rownames(rna.ol.all)
		rownames(meth_CNV_rprpp) <- rownames(rna.ol.all)
		rownames(meth_CNV_rprpr) <- rownames(rna.ol.all)
		rownames(meth_CNV_rprrp) <- rownames(rna.ol.all)
		rownames(meth_CNV_rprrr) <- rownames(rna.ol.all)
		
		rownames(meth_CNV_rrppp) <- rownames(rna.ol.all)
		print("rrppr")
		rownames(meth_CNV_rrppr) <- rownames(rna.ol.all)
		rownames(meth_CNV_rrprp) <- rownames(rna.ol.all)
		rownames(meth_CNV_rrprr) <- rownames(rna.ol.all)
		
		rownames(meth_CNV_rrrpp) <- rownames(rna.ol.all)
		rownames(meth_CNV_rrrpr) <- rownames(rna.ol.all)
		print("rrrpr")
		rownames(meth_CNV_rrrrp) <- rownames(rna.ol.all)
		rownames(meth_CNV_rrrrr) <- rownames(rna.ol.all)
		print("f")
		#print(dim(rsqr_until_cnv.ol))
		meth_CNV_rrrrr.ol <- meth_CNV_rrrrr
		meth_CNV_rrrrp.ol <- meth_CNV_rrrrp
		meth_CNV_rrrpr.ol <- meth_CNV_rrrpr
		meth_CNV_rrrpp.ol <- meth_CNV_rrrpp
		meth_CNV_rrprr.ol <- meth_CNV_rrprr
		meth_CNV_rrprp.ol <- meth_CNV_rrprp
		meth_CNV_rrppr.ol <- meth_CNV_rrppr
		meth_CNV_rrppp.ol <- meth_CNV_rrppp
		meth_CNV_rprrr.ol <- meth_CNV_rprrr
		meth_CNV_rprrp.ol <- meth_CNV_rprrp
		meth_CNV_rprpr.ol <- meth_CNV_rprpr
		meth_CNV_rprpp.ol <- meth_CNV_rprpp
		meth_CNV_rpprr.ol <- meth_CNV_rpprr
		meth_CNV_rpprp.ol <- meth_CNV_rpprp
		meth_CNV_rpppr.ol <- meth_CNV_rpppr
		meth_CNV_rpppp.ol <- meth_CNV_rpppp
		meth_CNV_prrrr.ol <- meth_CNV_prrrr
		meth_CNV_prrrp.ol <- meth_CNV_prrrp
		meth_CNV_prrpr.ol <- meth_CNV_prrpr
		meth_CNV_prrpp.ol <- meth_CNV_prrpp
		meth_CNV_prprr.ol <- meth_CNV_prprr
		meth_CNV_prprp.ol <- meth_CNV_prprp
		meth_CNV_prppr.ol <- meth_CNV_prppr
		meth_CNV_prppp.ol <- meth_CNV_prppp
		meth_CNV_pprrr.ol <- meth_CNV_pprrr
		meth_CNV_pprrp.ol <- meth_CNV_pprrp
		meth_CNV_pprpr.ol <- meth_CNV_pprpr
		meth_CNV_pprpp.ol <- meth_CNV_pprpp
		meth_CNV_ppprr.ol <- meth_CNV_ppprr
		meth_CNV_ppprp.ol <- meth_CNV_ppprp
		meth_CNV_ppppr.ol <- meth_CNV_ppppr
		meth_CNV_ppppp.ol <- meth_CNV_ppppp

		### CNV level
	  
		cnv.ol <- cbind(rrrrr, rrrrp, rrrpr, rrrpp,
						rrprr, rrprp, rrppr, rrppp,
						rprrr, rprrp, rprpr, rprpp,
						rpprr, rpprp, rpppr, rpppp,
						prrrr, prrrp, prrpr, prrpp,
						prprr, prprp, prppr, prppp,
						pprrr, pprrp, pprpr, pprpp,
						ppprr, ppprp, ppppr, ppppp)
		row.names(cnv) <- row.names(rna.ol.all)
		cnv_files.ol <- c("meth_CNV_rrrrr.ol", "meth_CNV_rrrrp.ol", "meth_CNV_rrrpr.ol", "meth_CNV_rrrpp.ol",
						"meth_CNV_rrprr.ol", "meth_CNV_rrprp.ol", "meth_CNV_rrppr.ol", "meth_CNV_rrppp.ol",
						"meth_CNV_rprrr.ol", "meth_CNV_rprrp.ol", "meth_CNV_rprpr.ol", "meth_CNV_rprpp.ol",
						"meth_CNV_rpprr.ol", "meth_CNV_rpprp.ol", "meth_CNV_rpppr.ol", "meth_CNV_rpppp.ol",
						"meth_CNV_prrrr.ol", "meth_CNV_prrrp.ol", "meth_CNV_prrpr.ol", "meth_CNV_prrpp.ol",
						"meth_CNV_prprr.ol", "meth_CNV_prprp.ol", "meth_CNV_prppr.ol", "meth_CNV_prppp.ol",
						"meth_CNV_pprrr.ol", "meth_CNV_pprrp.ol", "meth_CNV_pprpr.ol", "meth_CNV_pprpp.ol",
						"meth_CNV_ppprr.ol", "meth_CNV_ppprp.ol", "meth_CNV_ppppr.ol", "meth_CNV_ppppp.ol")
		if(isPar == T){
		##########
			return(list(meth_CNV_rrrrr.ol, meth_CNV_rrrrp.ol, meth_CNV_rrrpr.ol, meth_CNV_rrrpp.ol,
				meth_CNV_rrprr.ol, meth_CNV_rrprp.ol, meth_CNV_rrppr.ol, meth_CNV_rrppp.ol,
				meth_CNV_rprrr.ol, meth_CNV_rprrp.ol, meth_CNV_rprpr.ol, meth_CNV_rprpp.ol,
				meth_CNV_rpprr.ol, meth_CNV_rpprp.ol, meth_CNV_rpppr.ol, meth_CNV_rpppp.ol,
				meth_CNV_prrrr.ol, meth_CNV_prrrp.ol, meth_CNV_prrpr.ol, meth_CNV_prrpp.ol,
				meth_CNV_prprr.ol, meth_CNV_prprp.ol, meth_CNV_prppr.ol, meth_CNV_prppp.ol,
				meth_CNV_pprrr.ol, meth_CNV_pprrp.ol, meth_CNV_pprpr.ol, meth_CNV_pprpp.ol,
				meth_CNV_ppprr.ol, meth_CNV_ppprp.ol, meth_CNV_ppppr.ol, meth_CNV_ppppp.ol,
				cnv.ol))
		#
		} else {
			save(list = cnv_files.ol, file = paste0(cancer_name, "-four-", as.character(k), ".RData"))
			save(list = "cnv.ol",  file = paste0(cancer_name, "-four-rsqrs-", as.character(k), ".RData"))
			return(cnv.ol)
		}
	}
	if(stepNum == "five") {
		load(paste0(cancer_name, "-four-", as.character(k), ".RData"))
				meth_CNV_rrrrr.ol <- meth_CNV_rrrrr.ol[rownames(meth_CNV_rrrrr.ol) %in% rownames(rna.ol.all),]
		meth_CNV_rrrrp.ol <- meth_CNV_rrrrp.ol[rownames(meth_CNV_rrrrp.ol) %in% rownames(rna.ol.all),]
		meth_CNV_rrrpr.ol <- meth_CNV_rrrpr.ol[rownames(meth_CNV_rrrpr.ol) %in% rownames(rna.ol.all),]
		meth_CNV_rrrpp.ol <- meth_CNV_rrrpp.ol[rownames(meth_CNV_rrrpp.ol) %in% rownames(rna.ol.all),]
		meth_CNV_rrprr.ol <- meth_CNV_rrprr.ol[rownames(meth_CNV_rrprr.ol) %in% rownames(rna.ol.all),]
		meth_CNV_rrprp.ol <- meth_CNV_rrprp.ol[rownames(meth_CNV_rrprp.ol) %in% rownames(rna.ol.all),]
		meth_CNV_rrppr.ol <- meth_CNV_rrppr.ol[rownames(meth_CNV_rrppr.ol) %in% rownames(rna.ol.all),]
		meth_CNV_rrppp.ol <- meth_CNV_rrppp.ol[rownames(meth_CNV_rrppp.ol) %in% rownames(rna.ol.all),]
		meth_CNV_rprrr.ol <- meth_CNV_rprrr.ol[rownames(meth_CNV_rprrr.ol) %in% rownames(rna.ol.all),]
		meth_CNV_rprrp.ol <- meth_CNV_rprrp.ol[rownames(meth_CNV_rprrp.ol) %in% rownames(rna.ol.all),]
		meth_CNV_rprpr.ol <- meth_CNV_rprpr.ol[rownames(meth_CNV_rprpr.ol) %in% rownames(rna.ol.all),]
		meth_CNV_rprpp.ol <- meth_CNV_rprpp.ol[rownames(meth_CNV_rprpp.ol) %in% rownames(rna.ol.all),]
		meth_CNV_rpprr.ol <- meth_CNV_rpprr.ol[rownames(meth_CNV_rpprr.ol) %in% rownames(rna.ol.all),]
		meth_CNV_rpprp.ol <- meth_CNV_rpprp.ol[rownames(meth_CNV_rpprp.ol) %in% rownames(rna.ol.all),]
		meth_CNV_rpppr.ol <- meth_CNV_rpppr.ol[rownames(meth_CNV_rpppr.ol) %in% rownames(rna.ol.all),]
		meth_CNV_rpppp.ol <- meth_CNV_rpppp.ol[rownames(meth_CNV_rpppp.ol) %in% rownames(rna.ol.all),]
		meth_CNV_prrrr.ol <- meth_CNV_prrrr.ol[rownames(meth_CNV_prrrr.ol) %in% rownames(rna.ol.all),]
		meth_CNV_prrrp.ol <- meth_CNV_prrrp.ol[rownames(meth_CNV_prrrp.ol) %in% rownames(rna.ol.all),]
		meth_CNV_prrpr.ol <- meth_CNV_prrpr.ol[rownames(meth_CNV_prrpr.ol) %in% rownames(rna.ol.all),]
		meth_CNV_prrpp.ol <- meth_CNV_prrpp.ol[rownames(meth_CNV_prrpp.ol) %in% rownames(rna.ol.all),]
		meth_CNV_prprr.ol <- meth_CNV_prprr.ol[rownames(meth_CNV_prprr.ol) %in% rownames(rna.ol.all),]
		meth_CNV_prprp.ol <- meth_CNV_prprp.ol[rownames(meth_CNV_prprp.ol) %in% rownames(rna.ol.all),]
		meth_CNV_prppr.ol <- meth_CNV_prppr.ol[rownames(meth_CNV_prppr.ol) %in% rownames(rna.ol.all),]
		meth_CNV_prppp.ol <- meth_CNV_prppp.ol[rownames(meth_CNV_prppp.ol) %in% rownames(rna.ol.all),]
		meth_CNV_pprrr.ol <- meth_CNV_pprrr.ol[rownames(meth_CNV_pprrr.ol) %in% rownames(rna.ol.all),]
		meth_CNV_pprrp.ol <- meth_CNV_pprrp.ol[rownames(meth_CNV_pprrp.ol) %in% rownames(rna.ol.all),]
		meth_CNV_pprpr.ol <- meth_CNV_pprpr.ol[rownames(meth_CNV_pprpr.ol) %in% rownames(rna.ol.all),]
		meth_CNV_pprpp.ol <- meth_CNV_pprpp.ol[rownames(meth_CNV_pprpp.ol) %in% rownames(rna.ol.all),]
		meth_CNV_ppprr.ol <- meth_CNV_ppprr.ol[rownames(meth_CNV_ppprr.ol) %in% rownames(rna.ol.all),]
		meth_CNV_ppprp.ol <- meth_CNV_ppprp.ol[rownames(meth_CNV_ppprp.ol) %in% rownames(rna.ol.all),]
		meth_CNV_ppppr.ol <- meth_CNV_ppppr.ol[rownames(meth_CNV_ppppr.ol) %in% rownames(rna.ol.all),]
		meth_CNV_ppppp.ol <- meth_CNV_ppppp.ol[rownames(meth_CNV_ppppp.ol) %in% rownames(rna.ol.all),]
		snp_ppppp <- glmStep(five, meth_CNV_ppppp.ol)
		snp_pppppp.ol <- snp_ppppp$predicted
		snp_pppppr.ol <- snp_ppppp$residual
		pppppp <- snp_ppppp$Rsqr
	  
		snp_ppppr <- glmStep(five, meth_CNV_ppppr.ol)
		snp_pppprp.ol <- snp_ppppr$predicted
		snp_pppprr.ol <- snp_ppppr$residual
		pppprp <- snp_ppppr$Rsqr
		
		snp_ppprp <- glmStep(five, meth_CNV_ppprp.ol)
		snp_ppprpp.ol <- snp_ppprp$predicted
		snp_ppprpr.ol <- snp_ppprp$residual
		ppprpp <- snp_ppprp$Rsqr
	  
		snp_ppprr <- glmStep(five, meth_CNV_ppprr.ol)
		snp_ppprrp.ol <- snp_ppprr$predicted
		snp_ppprrr.ol <- snp_ppprr$residual
		ppprrp <- snp_ppprr$Rsqr
		###########################################
		snp_pprpp <- glmStep(five, meth_CNV_pprpp.ol)
		snp_pprppp.ol <- snp_pprpp$predicted
		snp_pprppr.ol <- snp_pprpp$residual
		pprppp <- snp_pprpp$Rsqr
		
		snp_pprpr <- glmStep(five, meth_CNV_pprpr.ol)
		snp_pprprp.ol <- snp_pprpr$predicted
		snp_pprprr.ol <- snp_pprpr$residual
		pprprp <- snp_pprpr$Rsqr
		
		snp_pprrp <- glmStep(five, meth_CNV_pprrp.ol)
		snp_pprrpp.ol <- snp_pprrp$predicted
		snp_pprrpr.ol <- snp_pprrp$residual
		pprrpp <- snp_pprrp$Rsqr
		
		snp_pprrr <- glmStep(five, meth_CNV_pprrr.ol)
		snp_pprrrp.ol <- snp_pprrr$predicted
		snp_pprrrr.ol <- snp_pprrr$residual
		pprrrp <- snp_pprrr$Rsqr
		#####################################
		#####################################
		snp_prppp <- glmStep(five, meth_CNV_prppp.ol)
		snp_prpppp.ol <- snp_prppp$predicted
		snp_prpppr.ol <- snp_prppp$residual
		prpppp <- snp_prppp$Rsqr
		
		snp_prppr <- glmStep(five, meth_CNV_prppr.ol)
		snp_prpprp.ol <- snp_prppr$predicted
		snp_prpprr.ol <- snp_prppr$residual
		prpprp <- snp_prppr$Rsqr
		
		snp_prprp <- glmStep(five, meth_CNV_prprp.ol)
		snp_prprpp.ol <- snp_prprp$predicted
		snp_prprpr.ol <- snp_prprp$residual
		prprpp <- snp_prprp$Rsqr
		
		snp_prprr <- glmStep(five, meth_CNV_prprr.ol)
		snp_prprrp.ol <- snp_prprr$predicted
		snp_prprrr.ol <- snp_prprr$residual
		prprrp <- snp_prprr$Rsqr
		###########################################
		snp_prrpp <- glmStep(five, meth_CNV_prrpp.ol)
		snp_prrppp.ol <- snp_prrpp$predicted
		snp_prrppr.ol <- snp_prrpp$residual
		prrppp <- snp_prrpp$Rsqr
		
		snp_prrpr <- glmStep(five, meth_CNV_prrpr.ol)
		snp_prrprp.ol <- snp_prrpr$predicted
		snp_prrprr.ol <- snp_prrpr$residual
		prrprp <- snp_prrpr$Rsqr
		
		snp_prrrp <- glmStep(five, meth_CNV_prrrp.ol)
		snp_prrrpp.ol <- snp_prrrp$predicted
		snp_prrrpr.ol <- snp_prrrp$residual
		prrrpp <- snp_prrrp$Rsqr
		
		snp_prrrr <- glmStep(five, meth_CNV_prrrr.ol)
		snp_prrrrp.ol <- snp_prrrr$predicted
		snp_prrrrr.ol <- snp_prrrr$residual
		prrrrp <- snp_prrrr$Rsqr
		
		
		snp_rpppp <- glmStep(five, meth_CNV_rpppp.ol)
		snp_rppppp.ol <- snp_rpppp$predicted
		snp_rppppr.ol <- snp_rpppp$residual
		rppppp <- snp_rpppp$Rsqr
		
		snp_rpppr <- glmStep(five, meth_CNV_rpppr.ol)
		snp_rppprp.ol <- snp_rpppr$predicted
		snp_rppprr.ol <- snp_rpppr$residual
		rppprp <- snp_rpppr$Rsqr
		
		snp_rpprp <- glmStep(five, meth_CNV_rpprp.ol)
		snp_rpprpp.ol <- snp_rpprp$predicted
		snp_rpprpr.ol <- snp_rpprp$residual
		rpprpp <- snp_rpprp$Rsqr
		
		snp_rpprr <- glmStep(five, meth_CNV_rpprr.ol)
		snp_rpprrp.ol <- snp_rpprr$predicted
		snp_rpprrr.ol <- snp_rpprr$residual
		rpprrp <- snp_rpprr$Rsqr
		###########################################
		snp_rprpp <- glmStep(five, meth_CNV_rprpp.ol)
		snp_rprppp.ol <- snp_rprpp$predicted
		snp_rprppr.ol <- snp_rprpp$residual
		rprppp <- snp_rprpp$Rsqr
		
		snp_rprpr <- glmStep(five, meth_CNV_rprpr.ol)
		snp_rprprp.ol <- snp_rprpr$predicted
		snp_rprprr.ol <- snp_rprpr$residual
		rprprp <- snp_rprpr$Rsqr
		
		snp_rprrp <- glmStep(five, meth_CNV_rprrp.ol)
		snp_rprrpp.ol <- snp_rprrp$predicted
		snp_rprrpr.ol <- snp_rprrp$residual
		rprrpp <- snp_rprrp$Rsqr
		
		snp_rprrr <- glmStep(five, meth_CNV_rprrr.ol)
		snp_rprrrp.ol <- snp_rprrr$predicted
		snp_rprrrr.ol <- snp_rprrr$residual
		rprrrp <- snp_rprrr$Rsqr
		#####################################
		#####################################
		snp_rrppp <- glmStep(five, meth_CNV_rrppp.ol)
		snp_rrpppp.ol <- snp_rrppp$predicted
		snp_rrpppr.ol <- snp_rrppp$residual
		rrpppp <- snp_rrppp$Rsqr
		
		snp_rrppr <- glmStep(five, meth_CNV_rrppr.ol)
		snp_rrpprp.ol <- snp_rrppr$predicted
		snp_rrpprr.ol <- snp_rrppr$residual
		rrpprp <- snp_rrppr$Rsqr
		
		snp_rrprp <- glmStep(five, meth_CNV_rrprp.ol)
		snp_rrprpp.ol <- snp_rrprp$predicted
		snp_rrprpr.ol <- snp_rrprp$residual
		rrprpp <- snp_rrprp$Rsqr
		
		snp_rrprr <- glmStep(five, meth_CNV_rrprr.ol)
		snp_rrprrp.ol <- snp_rrprr$predicted
		snp_rrprrr.ol <- snp_rrprr$residual
		rrprrp <- snp_rrprr$Rsqr
		###########################################
		snp_rrrpp <- glmStep(five, meth_CNV_rrrpp.ol)
		snp_rrrppp.ol <- snp_rrrpp$predicted
		snp_rrrppr.ol <- snp_rrrpp$residual
		rrrppp <- snp_rrrpp$Rsqr
		
		snp_rrrpr <- glmStep(five, meth_CNV_rrrpr.ol)
		snp_rrrprp.ol <- snp_rrrpr$predicted
		snp_rrrprr.ol <- snp_rrrpr$residual
		rrrprp <- snp_rrrpr$Rsqr
		
		snp_rrrrp <- glmStep(five, meth_CNV_rrrrp.ol)
		snp_rrrrpp.ol <- snp_rrrrp$predicted
		snp_rrrrpr.ol <- snp_rrrrp$residual
		rrrrpp <- snp_rrrrp$Rsqr
		
		snp_rrrrr <- glmStep(five, meth_CNV_rrrrr.ol)
		snp_rrrrrp.ol <- snp_rrrrr$predicted
		snp_rrrrrr.ol <- snp_rrrrr$residual
		rrrrrp <- snp_rrrrr$Rsqr
		
		
		#######################################
		#######################################
		pppppr <- 1 - pppppp
		pppprr <- 1 - pppprp
		ppprpr <- 1 - ppprpp
		ppprrr <- 1 - ppprrp
		pprppr <- 1 - pprppp
		pprprr <- 1 - pprprp
		pprrpr <- 1 - pprrpp
		pprrrr <- 1 - pprrrp
		prpppr <- 1 - prpppp
		prpprr <- 1 - prpprp
		prprpr <- 1 - prprpp
		prprrr <- 1 - prprrp
		prrppr <- 1 - prrppp
		prrprr <- 1 - prrprp
		prrrpr <- 1 - prrrpp
		prrrrr <- 1 - prrrrp
		################################################
		rppppr <- 1 - rppppp
		rppprr <- 1 - rppprp
		rpprpr <- 1 - rpprpp
		rpprrr <- 1 - rpprrp
		rprppr <- 1 - rprppp
		rprprr <- 1 - rprprp
		rprrpr <- 1 - rprrpp
		rprrrr <- 1 - rprrrp
		rrpppr <- 1 - rrpppp
		rrpprr <- 1 - rrpprp
		rrprpr <- 1 - rrprpp
		rrprrr <- 1 - rrprrp
		rrrppr <- 1 - rrrppp
		rrrprr <- 1 - rrrprp
		rrrrpr <- 1 - rrrrpp
		rrrrrr <- 1 - rrrrrp

		snp_rem.ol <- cbind(rrrrrr, rrrrrp, rrrrpr, rrrrpp,
						rrrprr, rrrprp, rrrppr, rrrppp,
						rrprrr, rrprrp, rrprpr, rrprpp,
						rrpprr, rrpprp, rrpppr, rrpppp,
						rprrrr, rprrrp, rprrpr, rprrpp,
						rprprr, rprprp, rprppr, rprppp,
						rpprrr, rpprrp, rpprpr, rpprpp,
						rppprr, rppprp, rppppr, rppppp,
						prrrrr, prrrrp, prrrpr, prrrpp,
						prrprr, prrprp, prrppr, prrppp,
						prprrr, prprrp, prprpr, prprpp,
						prpprr, prpprp, prpppr, prpppp,
						pprrrr, pprrrp, pprrpr, pprrpp,
						pprprr, pprprp, pprppr, pprppp,
						ppprrr, ppprrp, ppprpr, ppprpp,
						pppprr, pppprp, pppppr, pppppp)
		rownames(snp_rem.ol) <- rownames(rna.ol.all)
		if(isPar == T){
			return(snp_rem.ol)
		} else {
			# load all of the R^2 values from the previous steps and combine them together
			
			load(paste0(cancer_name, "-one-rsqrs-", as.character(k), ".RData"))
			load(paste0(cancer_name, "-two-rsqrs-", as.character(k), ".RData"))
			load(paste0(cancer_name, "-three-rsqrs-", as.character(k), ".RData"))
			load(paste0(cancer_name, "-four-rsqrs-", as.character(k), ".RData"))
			
			#combine them into one final data fram and save it
			finalTree.ol <- cbind(snp_rem.ol, cnv.ol, methyl.ol, TF.ol, miRNA.ol, env.ol)
			save(list = "finalTree.ol", file = paste0(cancer_name, "-rsqrs-", as.character(k),".RData"))
		}

	}
}


saveParallels <- function(x){
	# This function takes a list that did not recombine through the foreach function
	# then, it combines the list based on the level of modality
	# this only needs to be done if the function is run in parallel
	if(stepNum == "one"){  #for the first step
		env_miRNA_rr.ol <- c() #make the names that fit in with the second step
		
		for(i in 1:n.cores){ #for all of the cores that you're using 
			env_miRNA_rr.ol <- rbind(env_miRNA_rr.ol, x[[i]]) #rbind those splits to rebuild the matrix at that level
		}
		x <- x[n.cores+1:length(x)] # take the first matriex (i.e. env_miRNA_rr.ol) from the list oof split matricies
		
		##what is above is repeated for every other step
		
		env_miRNA_rp.ol <- c()
		for(i in 1:n.cores){
			env_miRNA_rp.ol <- rbind(env_miRNA_rp.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		env_miRNA_pr.ol <- c()
		
		for(i in 1:n.cores){
			env_miRNA_pr.ol <- rbind(env_miRNA_pr.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		env_miRNA_pp.ol <- c()
		for(i in 1:n.cores){
			env_miRNA_pp.ol <- rbind(env_miRNA_pp.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]
		
		miRNA.ol <- c()
		for(i in 1:n.cores){
			miRNA.ol <- rbind(miRNA.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		env.ol <- c()
		for(i in 1:n.cores){
			env.ol <- rbind(env.ol, x[[i]])
		}
		#now that all of the files have been combined properly, you want to save them as your datafiles
		#this will be repeated in the exact same manner for every step
		saveLvlOne.ol <- c("env_miRNA_pp.ol", "env_miRNA_pr.ol", "env_miRNA_rp.ol", "env_miRNA_rr.ol")
		rsqrs <- c("env.ol", "miRNA.ol")
		save(list = c("env.ol", "miRNA.ol"), file = paste0(cancer_name, "-one-rsqrs-", as.character(k), ".RData"))
		save(list = saveLvlOne.ol, file = paste0(cancer_name, "-one-", as.character(k), ".RData"))
		return("The first level is in your saved directory :)")
	}
	if(stepNum == "two"){
		TF_rrr.ol <- c()
		for(i in 1:n.cores){
			TF_rrr.ol <- rbind(TF_rrr.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		TF_rrp.ol <- c()
		for(i in 1:n.cores){
			TF_rrp.ol <- rbind(TF_rrp.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		TF_rpr.ol <- c()
		for(i in 1:n.cores){
			TF_rpr.ol <- rbind(TF_rpr.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		TF_rpp.ol <- c()
		for(i in 1:n.cores){
			TF_rpp.ol <- rbind(TF_rpp.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		TF_prr.ol <- c()
		for(i in 1:n.cores){
			TF_prr.ol <- rbind(TF_prr.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		TF_prp.ol <- c()
		for(i in 1:n.cores){
			TF_prp.ol <- rbind(TF_prp.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		TF_ppr.ol <- c()
		for(i in 1:n.cores){
			TF_ppr.ol <- rbind(TF_ppr.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		TF_ppp.ol <- c()
		for(i in 1:n.cores){
			TF_ppp.ol <- rbind(TF_ppp.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		TF.ol <- c()
		for(i in 1:n.cores){
			TF.ol <- rbind(TF.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		TF_files.ol <- c("TF_ppr.ol", "TF_prp.ol", "TF_prr.ol", "TF_ppp.ol", 
						"TF_rpp.ol", "TF_rpr.ol", "TF_rrp.ol", "TF_rrr.ol")
		rsqrs <- "TF.ol"
		save(list = TF_files.ol, file = paste0(cancer_name, "-two-", as.character(k), ".RData"))
		save(list = "TF.ol",  file = paste0(cancer_name, "-two-rsqrs-", as.character(k), ".RData"))
		return("The level is in your saved directoy :)")
	}
	if(stepNum == "three"){
		meth_rrrr.ol <- c()
		for(i in 1:n.cores){
			meth_rrrr.ol <- rbind(meth_rrrr.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_rrrp.ol <- c()
		for(i in 1:n.cores){
			meth_rrrp.ol <- rbind(meth_rrrp.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_rrpr.ol <- c()
		for(i in 1:n.cores){
			meth_rrpr.ol <- rbind(meth_rrpr.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_rrpp.ol <- c()
		for(i in 1:n.cores){
			meth_rrpp.ol <- rbind(meth_rrpp.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_rprr.ol <- c()
		for(i in 1:n.cores){
			meth_rprr.ol <- rbind(meth_rprr.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_rprp.ol <- c()
		for(i in 1:n.cores){
			meth_rprp.ol <- rbind(meth_rprp.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_rppr.ol <- c()
		for(i in 1:n.cores){
			meth_rppr.ol <- rbind(meth_rppr.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_rppp.ol <- c()
		for(i in 1:n.cores){
			meth_rppp.ol <- rbind(meth_rppp.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_prrr.ol <- c()
		for(i in 1:n.cores){
			meth_prrr.ol <- rbind(meth_prrr.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_prrp.ol <- c()
		for(i in 1:n.cores){
			meth_prrp.ol <- rbind(meth_prrp.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_prpr.ol <- c()
		for(i in 1:n.cores){
			meth_prpr.ol <- rbind(meth_prpr.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_prpp.ol <- c()
		for(i in 1:n.cores){
			meth_prpp.ol <- rbind(meth_prpp.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_pprr.ol <- c()
		for(i in 1:n.cores){
			meth_pprr.ol <- rbind(meth_pprr.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_pprp.ol <- c()
		for(i in 1:n.cores){
			meth_pprp.ol <- rbind(meth_pprp.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_pppr.ol <- c()
		for(i in 1:n.cores){
			meth_pppr.ol <- rbind(meth_pppr.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_pppp.ol <- c()
		for(i in 1:n.cores){
			meth_pppp.ol <- rbind(meth_pppp.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		methyl.ol <- c()
		for(i in 1:n.cores){
			methyl.ol <- rbind(methyl.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

	
		methnames.ol <- c("meth_pppp.ol", "meth_pppr.ol","meth_pprp.ol", "meth_pprr.ol",
						"meth_prpp.ol", "meth_prpr.ol","meth_prrp.ol", "meth_prrr.ol",
						"meth_rppp.ol", "meth_rppr.ol","meth_rprp.ol", "meth_rprr.ol",
						"meth_rrpp.ol", "meth_rrpr.ol","meth_rrrp.ol", "meth_rrrr.ol")
		rsqrs <- "methyl.ol"
		save(list = methnames.ol, file = paste0(cancer_name, "-three-", as.character(k), ".RData"))
		save(list = "methyl.ol",  file = paste0(cancer_name, "-three-rsqrs-", as.character(k), ".RData"))
		return("The level is in your saved directory :)")
	}
	if(stepNum == "four"){
		meth_CNV_rrrrr.ol <- c()
		for(i in 1:n.cores){
			meth_CNV_rrrrr.ol <- rbind(meth_CNV_rrrrr.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_CNV_rrrrp.ol <- c()
		for(i in 1:n.cores){
			meth_CNV_rrrrp.ol <- rbind(meth_CNV_rrrrp.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_CNV_rrrpr.ol <- c()
		for(i in 1:n.cores){
			meth_CNV_rrrpr.ol <- rbind(meth_CNV_rrrpr.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_CNV_rrrpp.ol <- c()
		for(i in 1:n.cores){
			meth_CNV_rrrpp.ol <- rbind(meth_CNV_rrrpp.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_CNV_rrprr.ol <- c()
		for(i in 1:n.cores){
			meth_CNV_rrprr.ol <- rbind(meth_CNV_rrprr.ol, x[[i]])
		}
		meth_CNV_rrprp.ol <- c()
		for(i in 1:n.cores){
			meth_CNV_rrprp.ol <- rbind(meth_CNV_rrprp.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_CNV_rrppr.ol <- c()
		for(i in 1:n.cores){
			meth_CNV_rrppr.ol <- rbind(meth_CNV_rrppr.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_CNV_rrppp.ol <- c()
		for(i in 1:n.cores){
			meth_CNV_rrppp.ol <- rbind(meth_CNV_rrppp.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_CNV_rprrr.ol <- c()
		for(i in 1:n.cores){
			meth_CNV_rprrr.ol <- rbind(meth_CNV_rprrr.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_CNV_rprrp.ol <- c()
		for(i in 1:n.cores){
			meth_CNV_rprrp.ol <- rbind(meth_CNV_rprrp.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_CNV_rprpr.ol <- c()
		for(i in 1:n.cores){
			meth_CNV_rprpr.ol <- rbind(meth_CNV_rprpr.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_CNV_rprpp.ol <- c()
		for(i in 1:n.cores){
			meth_CNV_rprpp.ol <- rbind(meth_CNV_rprpp.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_CNV_rpprr.ol <- c()
		for(i in 1:n.cores){
			meth_CNV_rpprr.ol <- rbind(meth_CNV_rpprr.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_CNV_rpprp.ol <- c()
		for(i in 1:n.cores){
			meth_CNV_rpprp.ol <- rbind(meth_CNV_rpprp.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_CNV_rpppr.ol <- c()
		for(i in 1:n.cores){
			meth_CNV_rpppr.ol <- rbind(meth_CNV_rpppr.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_CNV_rpppp.ol <- c()
		for(i in 1:n.cores){
			meth_CNV_rpppp.ol <- rbind(meth_CNV_rpppp.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_CNV_prrrr.ol <- c()
		for(i in 1:n.cores){
			meth_CNV_prrrr.ol <- rbind(meth_CNV_prrrr.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_CNV_prrrp.ol <- c()
		for(i in 1:n.cores){
			meth_CNV_prrrp.ol <- rbind(meth_CNV_prrrp.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_CNV_prrpr.ol <- c()
		for(i in 1:n.cores){
			meth_CNV_prrpr.ol <- rbind(meth_CNV_prrpr.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_CNV_prrpp.ol <- c()
		for(i in 1:n.cores){
			meth_CNV_prrpp.ol <- rbind(meth_CNV_prrpp.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_CNV_prprr.ol <- c()
		for(i in 1:n.cores){
			meth_CNV_prprr.ol <- rbind(meth_CNV_prprr.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_CNV_prprp.ol <- c()
		for(i in 1:n.cores){
			meth_CNV_prprp.ol <- rbind(meth_CNV_prprp.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_CNV_prppr.ol <- c()
		for(i in 1:n.cores){
			meth_CNV_prppr.ol <- rbind(meth_CNV_prppr.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_CNV_prppp.ol <- c()
		for(i in 1:n.cores){
			meth_CNV_prppp.ol <- rbind(meth_CNV_prppp.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_CNV_pprrr.ol <- c()
		for(i in 1:n.cores){
			meth_CNV_pprrr.ol <- rbind(meth_CNV_pprrr.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_CNV_pprrp.ol <- c()
		for(i in 1:n.cores){
			meth_CNV_pprrp.ol <- rbind(meth_CNV_pprrp.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_CNV_pprpr.ol <- c()
		for(i in 1:n.cores){
			meth_CNV_pprpr.ol <- rbind(meth_CNV_pprpr.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_CNV_pprpp.ol <- c()
		for(i in 1:n.cores){
			meth_CNV_pprpp.ol <- rbind(meth_CNV_pprpp.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_CNV_ppprr.ol <- c()
		for(i in 1:n.cores){
			meth_CNV_ppprr.ol <- rbind(meth_CNV_ppprr.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_CNV_ppprp.ol <- c()
		for(i in 1:n.cores){
			meth_CNV_ppprp.ol <- rbind(meth_CNV_ppprp.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_CNV_ppppr.ol <- c()
		for(i in 1:n.cores){
			meth_CNV_ppppr.ol <- rbind(meth_CNV_ppppr.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		meth_CNV_ppppp.ol <- c()
		for(i in 1:n.cores){
			meth_CNV_ppppp.ol <- rbind(meth_CNV_ppppp.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		cnv.ol <- c()
		for(i in 1:n.cores){
			cnv.ol <- rbind(cnv.ol, x[[i]])
		}
		x <- x[n.cores+1:length(x)]

		cnv_files.ol <- c("meth_CNV_rrrrr.ol", "meth_CNV_rrrrp.ol", "meth_CNV_rrrpr.ol", "meth_CNV_rrrpp.ol",
						"meth_CNV_rrprr.ol", "meth_CNV_rrprp.ol", "meth_CNV_rrppr.ol", "meth_CNV_rrppp.ol",
						"meth_CNV_rprrr.ol", "meth_CNV_rprrp.ol", "meth_CNV_rprpr.ol", "meth_CNV_rprpp.ol",
						"meth_CNV_rpprr.ol", "meth_CNV_rpprp.ol", "meth_CNV_rpppr.ol", "meth_CNV_rpppp.ol",
						"meth_CNV_prrrr.ol", "meth_CNV_prrrp.ol", "meth_CNV_prrpr.ol", "meth_CNV_prrpp.ol",
						"meth_CNV_prprr.ol", "meth_CNV_prprp.ol", "meth_CNV_prppr.ol", "meth_CNV_prppp.ol",
						"meth_CNV_pprrr.ol", "meth_CNV_pprrp.ol", "meth_CNV_pprpr.ol", "meth_CNV_pprpp.ol",
						"meth_CNV_ppprr.ol", "meth_CNV_ppprp.ol", "meth_CNV_ppppr.ol", "meth_CNV_ppppp.ol")
		rsqrs <- "cnv.ol"
		save(list = cnv_files.ol, file = paste0(cancer_name, "-four-", as.character(k), ".RData"))
		save(list = "cnv.ol",  file = paste0(cancer_name, "-four-rsqrs-", as.character(k), ".RData"))
		return("The files have been saved in your directoy :)")
	}
	if(stepNum == "five"){
		#for the final step you want to combine the R^2 values from the last step, then combine the 
		#R^2 values from the previous steps in the function, cbind them together, and then save the
		#final functions. 
		snp_rem.ol <- c()
		for(i in 1:n.cores){
			snp_rem.ol <- rbind(snp_rem.ol, x[[i]])
		
		}
		load(paste0(cancer_name, "-one-rsqrs-", as.character(k), ".RData"))
		load(paste0(cancer_name, "-two-rsqrs-", as.character(k), ".RData"))
		load(paste0(cancer_name, "-three-rsqrs-", as.character(k), ".RData"))
		load(paste0(cancer_name, "-four-rsqrs-", as.character(k), ".RData"))

		# combine all of the R^2 values that will be required to 
		# build the tree of predicted and residual values
	  
		finalTree.ol <- cbind(snp_rem.ol, cnv.ol, methyl.ol, TF.ol, miRNA.ol, env.ol)
		save(list = "finalTree.ol", file = paste0(cancer_name, "-rsqrs-", as.character(k),".RData"))
	}

}


