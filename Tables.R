## Cold Case Kaszubinski et al. 
# Tables

#Packages
library(car)
library(exactRankTests)
library(nlme)
library(GGally)
library(ggforce)
library(ggplot2)
library(ggpubr)
library(ggRandomForests)
library(ggthemes)
library(lme4)
library(phyloseq)
library(plyr)
library(PMCMR)
library(randomForest)
library(rsample)
library(tidyverse)
library(vegan)
set.seed(1234)

# Load Microbial Data
otu <- read.csv("otutable_fp_otuid.csv")
tax <- read.csv("taxtable_fp_otuid.csv")
tree <- read_tree('tree.nwk')

metadata=(read.csv("ForensicShowPigMetadata.csv",header=TRUE))
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$SampleID
sampdat$Day <- factor(sampdat$Day, levels = c('1', '5', '9', '13', '17', '22'))

metadata_complete <- metadata[complete.cases(metadata),]
metadata_complete %>% group_by(Day) %>% summarise_at(c('Temp', 'ADH'), funs(mean, sd))


rownames(otu) <- otu$OTUID
otu <- otu[,-1]
OTU=otu_table(otu, taxa_are_rows=TRUE)
physeq_otu.tree=phyloseq(OTU,tree, sampdat)

rownames(tax) <- tax$OTUID
tax_reduc <- merge(tax, otu, by = "row.names")
rownames(tax_reduc) <- tax_reduc$Row.names
tax_reduc <- tax_reduc[,-1]
tax_reduc <- tax_reduc[,-1]

tax_f <- tax_reduc[,1:7]
tax_f <- as.matrix(tax_f)
TAX=tax_table(tax_f)
taxa_names(TAX)=row.names(OTU)

physeq_otu.tax=phyloseq(OTU,TAX, sampdat)
#physeq object with tree
physeq_beta <- rarefy_even_depth(physeq_otu.tree, sample.size = 5000)
#physeq object without tree
physeq <- rarefy_even_depth(physeq_otu.tax, sample.size = 5000)

## Table 1

# This table was created and calculated in excel

## Table 2
 
# This table was created and calculated in excel

## Table 3
# Differentially abundant families among PMSI days identified by ANCOM. 
# Pairwise comparison among days are included as well as the number of significant 
# families identified with the corresponding phyla. 
# The average W stat among differentially abundant taxa is also reported. 

# ANCOM function 
# code from https://sites.google.com/site/siddharthamandal1985/research 
ancom.W = function(otu_data,var_data,
                   adjusted,repeated,
                   main.var,adj.formula,
                   repeat.var,long,rand.formula,
                   multcorr,sig){
  
  n_otu=dim(otu_data)[2]-1
  
  otu_ids=colnames(otu_data)[-1]
  
  if(repeated==F){
    data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID",all.y=T),row.names=NULL)
    #data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var)],by="Sample.ID",all.y=T),row.names=NULL)
  }else if(repeated==T){
    data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID"),row.names=NULL)
    # data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var,repeat.var)],by="Sample.ID"),row.names=NULL)
  }
  
  base.formula = paste0("lr ~ ",main.var)
  if(repeated==T){
    repeat.formula = paste0(base.formula," | ", repeat.var)
  }
  if(adjusted==T){
    adjusted.formula = paste0(base.formula," + ", adj.formula)
  }
  
  if( adjusted == F & repeated == F ){
    fformula  <- formula(base.formula)
  } else if( adjusted == F & repeated == T & long == T ){
    fformula  <- formula(base.formula)   
  }else if( adjusted == F & repeated == T & long == F ){
    fformula  <- formula(repeat.formula)   
  }else if( adjusted == T & repeated == F  ){
    fformula  <- formula(adjusted.formula)   
  }else if( adjusted == T & repeated == T  ){
    fformula  <- formula(adjusted.formula)   
  }else{
    stop("Problem with data. Dataset should contain OTU abundances, groups, 
         and optionally an ID for repeated measures.")
  }
  
  
  
  if( repeated==FALSE & adjusted == FALSE){
    if( length(unique(data_comp[,which(colnames(data_comp)==main.var)]))==2 ){
      tfun <- exactRankTests::wilcox.exact
    } else{
      tfun <- stats::kruskal.test
    }
  }else if( repeated==FALSE & adjusted == TRUE){
    tfun <- stats::aov
  }else if( repeated== TRUE & adjusted == FALSE & long == FALSE){
    tfun <- stats::friedman.test
  }else if( repeated== TRUE & adjusted == FALSE & long == TRUE){
    tfun <- nlme::lme
  }else if( repeated== TRUE & adjusted == TRUE){
    tfun <- nlme::lme
  }
  
  logratio.mat <- matrix(NA, nrow=n_otu, ncol=n_otu)
  for(ii in 1:(n_otu-1)){
    for(jj in (ii+1):n_otu){
      data.pair <- data_comp[,which(colnames(data_comp)%in%otu_ids[c(ii,jj)])]
      lr <- log((1+as.numeric(data.pair[,1]))/(1+as.numeric(data.pair[,2])))
      
      lr_dat <- data.frame( lr=lr, data_comp,row.names=NULL )
      
      if(adjusted==FALSE&repeated==FALSE){  ## Wilcox, Kruskal Wallis
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }else if(adjusted==FALSE&repeated==TRUE&long==FALSE){ ## Friedman's 
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }else if(adjusted==TRUE&repeated==FALSE){ ## ANOVA
        model=tfun(formula=fformula, data = lr_dat,na.action=na.omit)   
        picker=which(gsub(" ","",row.names(summary(model)[[1]]))==main.var)  
        logratio.mat[ii,jj] <- summary(model)[[1]][["Pr(>F)"]][picker]
      }else if(repeated==TRUE&long==TRUE){ ## GEE
        model=tfun(fixed=fformula,data = lr_dat,
                   random = formula(rand.formula),
                   correlation=corAR1(),
                   na.action=na.omit)   
        picker=which(gsub(" ","",row.names(anova(model)))==main.var)
        logratio.mat[ii,jj] <- anova(model)[["p-value"]][picker]
      }
      
    }
  } 
  
  ind <- lower.tri(logratio.mat)
  logratio.mat[ind] <- t(logratio.mat)[ind]
  
  
  logratio.mat[which(is.finite(logratio.mat)==FALSE)] <- 1
  
  mc.pval <- t(apply(logratio.mat,1,function(x){
    s <- p.adjust(x, method = "BH")
    return(s)
  }))
  
  a <- logratio.mat[upper.tri(logratio.mat,diag=FALSE)==TRUE]
  
  b <- matrix(0,ncol=n_otu,nrow=n_otu)
  b[upper.tri(b)==T] <- p.adjust(a, method = "BH")
  diag(b)  <- NA
  ind.1    <- lower.tri(b)
  b[ind.1] <- t(b)[ind.1]
  
  #########################################
  ### Code to extract surrogate p-value
  surr.pval <- apply(mc.pval,1,function(x){
    s0=quantile(x[which(as.numeric(as.character(x))<sig)],0.95)
    # s0=max(x[which(as.numeric(as.character(x))<alpha)])
    return(s0)
  })
  #########################################
  ### Conservative
  if(multcorr==1){
    W <- apply(b,1,function(x){
      subp <- length(which(x<sig))
    })
    ### Moderate
  } else if(multcorr==2){
    W <- apply(mc.pval,1,function(x){
      subp <- length(which(x<sig))
    })
    ### No correction
  } else if(multcorr==3){
    W <- apply(logratio.mat,1,function(x){
      subp <- length(which(x<sig))
    })
  }
  
  return(W)
  }



ANCOM.main = function(OTUdat,Vardat,
                      adjusted,repeated,
                      main.var,adj.formula,
                      repeat.var,longitudinal,
                      random.formula,
                      multcorr,sig,
                      prev.cut){
  
  p.zeroes=apply(OTUdat[,-1],2,function(x){
    s=length(which(x==0))/length(x)
  })
  
  zeroes.dist=data.frame(colnames(OTUdat)[-1],p.zeroes,row.names=NULL)
  colnames(zeroes.dist)=c("Taxon","Proportion_zero")
  
  zero.plot = ggplot(zeroes.dist, aes(x=Proportion_zero)) + 
    geom_histogram(binwidth=0.1,colour="black",fill="white") + 
    xlab("Proportion of zeroes") + ylab("Number of taxa") +
    theme_bw()
  
  #print(zero.plot)
  
  OTUdat.thinned=OTUdat
  OTUdat.thinned=OTUdat.thinned[,c(1,1+which(p.zeroes<prev.cut))]
  
  otu.names=colnames(OTUdat.thinned)[-1]
  
  W.detected   <- ancom.W(OTUdat.thinned,Vardat,
                          adjusted,repeated,
                          main.var,adj.formula,
                          repeat.var,longitudinal,random.formula,
                          multcorr,sig)
  
  W_stat       <- W.detected
  
  
  ### Bubble plot
  
  W_frame = data.frame(otu.names,W_stat,row.names=NULL)
  W_frame = W_frame[order(-W_frame$W_stat),]
  
  W_frame$detected_0.9=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.8=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.7=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.6=rep(FALSE,dim(W_frame)[1])
  
  W_frame$detected_0.9[which(W_frame$W_stat>0.9*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.8[which(W_frame$W_stat>0.8*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.7[which(W_frame$W_stat>0.7*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.6[which(W_frame$W_stat>0.6*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  
  final_results=list(W_frame,zero.plot)
  names(final_results)=c("W.taxa","PLot.zeroes")
  return(final_results)
}

# subsetting phyloseq objects in pairwise day comparisons
day_pairwise <- function(physeq) {
  physeq_1 <- subset_samples(physeq, Day == '1')
  physeq_5 <- subset_samples(physeq, Day == '5')
  physeq_9 <- subset_samples(physeq, Day == '9')
  physeq_13 <- subset_samples(physeq, Day == '13')
  physeq_17 <- subset_samples(physeq, Day == '17')
  physeq_22 <- subset_samples(physeq, Day == '22')
  return(list(physeq_15 <- merge_phyloseq(physeq_1, physeq_5),
              physeq_19 <- merge_phyloseq(physeq_1, physeq_9),
              physeq_113 <- merge_phyloseq(physeq_1, physeq_13),
              physeq_117 <- merge_phyloseq(physeq_1, physeq_17),
              physeq_131 <- merge_phyloseq(physeq_1, physeq_22),
              physeq_59 <- merge_phyloseq(physeq_9, physeq_5),
              physeq_513 <- merge_phyloseq(physeq_13, physeq_5),
              physeq_517 <- merge_phyloseq(physeq_17, physeq_5),
              physeq_522 <- merge_phyloseq(physeq_22, physeq_5),
              physeq_913 <- merge_phyloseq(physeq_9, physeq_13),
              physeq_917 <- merge_phyloseq(physeq_9, physeq_17),
              physeq_922 <- merge_phyloseq(physeq_9, physeq_22),
              physeq_1317 <- merge_phyloseq(physeq_17, physeq_13),
              physeq_1322 <- merge_phyloseq(physeq_17, physeq_22),
              physeq_1722 <- merge_phyloseq(physeq_17, physeq_22)))
}

#subset phyloseq object by family 
physeq_fam <- tax_glom(physeq, taxrank = 'Family')

#create list of all day pairwise comparisons
list_for_ancom_fam <- day_pairwise(physeq_fam)

#create function to pull OTU table for ANCOM analysis, and put it into usable format
# for function
otu_ancom_make <- function(physeq) {
  otu_ancom <- data.frame(otu_table(physeq))
  otu_ancom <- data.frame(t(otu_ancom))
  Sample.ID <- rownames(otu_ancom)
  rownames(otu_ancom) <- NULL
  otu_ancom <- cbind(Sample.ID, otu_ancom)
  return(otu_ancom)
}

# add metadata into correct format for ANCOM function
metadata_ancom <- metadata
colnames(metadata_ancom)[1] <- 'Sample.ID'

#make blank list to collect results
ancom_results_fam <- list()

for(i in 1:length(list_for_ancom_fam)) {
  otu_ancom <- otu_ancom_make(list_for_ancom_fam[[i]])
  comparison_test <- ANCOM.main(OTUdat = otu_ancom,
                                Vardat = metadata_ancom,
                                adjusted = FALSE,
                                repeated = F,
                                main.var = "Day",
                                adj.formula = NULL,
                                repeat.var=NULL,
                                longitudinal=FALSE,
                                random.formula=NULL,
                                multcorr=2,
                                sig=0.05,
                                prev.cut=0.90)
  w_values <- data.frame(comparison_test$W.taxa)
  tax <- data.frame(tax_table(physeq))
  tax <- tax %>% select(Kingdom, Phylum, Order, Class, Family)
  tax$otu.names <- rownames(tax)
  ancom_sign_taxa <- merge(w_values, tax, by='otu.names')
  ancom_sign_taxa <- ancom_sign_taxa[,-1]
  ancom_results_fam[[i]] <- ancom_sign_taxa
}

# Table was ultimately organized and created in excel

## Table 4
# Random Forest regression among days. 
# Random forest regressions were made with 500 trees.
# Models were created with out-of-bag (OOB) error rates as well as test set error rates,
# including all sampling locations, mouth and ventral sites, and dorsal sites.

# get OTU table into correct format for RF
otu <- as.data.frame(t(otu_table(physeq)))
otu$SampleID <- rownames(otu)
#metadata reduced to only necessary data
meta_sa <- metadata %>% select(SampleID, Day, Location)
otu <- merge(meta_sa, otu, by = 'SampleID')
otu <- otu[,-1]
names(otu) <- make.names(names(otu))

# OOB error
# multiple trees were tested, error leveled off after 500 trees
m1 <- randomForest(
  formula = Day ~ .,
  data    = otu,
  ntree= 500
)

#model results
m1
#mse
sqrt(m1$mse[which.min(m1$mse)])

#comparing OOB model to test/ training set
valid_split <- initial_split(otu, .7)
otu_train <- analysis(valid_split)
otu_valid <- assessment(valid_split)
x_test <- otu_valid[setdiff(names(otu_valid), "Day")]
y_test <- otu_valid$Day

rf_oob_comp <- randomForest(
  formula = Day~ .,
  data    = otu_train,
  xtest   = x_test,
  ytest   = y_test,
  ntree = 500
)

# compare error rates of OOB vs. test/ training set
oob <- sqrt(m1$mse)
validation <- sqrt(rf_oob_comp$mse)

tibble::tibble(
  `Out of Bag Error` = oob,
  `Test error` = validation,
  ntrees = 1:rf_oob_comp$ntree
) %>%
  gather(Metric, Error, -ntrees) %>%
  ggplot(aes(ntrees, Error, color = Metric)) +
  geom_line() +
  xlab("Number of trees")
#look at model
rf_oob_comp
#MSE 
sqrt(m1$mse[which.min(rf_oob_comp$mse)])

# Same process as before
# But, using mouth and ventral sites only
otu_l1 <- otu %>%  filter(Location == 'L1')
otu_l2 <- otu %>%  filter(Location == 'L2')
otu_l3 <- otu %>%  filter(Location == 'L3')

otu_l1.2.3 <- rbind(otu_l1, otu_l2, otu_l3)

m1 <- randomForest(
  formula = Day ~ .,
  data    = otu_l1.2.3,
  ntree= 500
)

m1
sqrt(m1$mse[which.min(m1$mse)])

valid_split <- initial_split(otu_l1.2.3, .8)
otu_train <- analysis(valid_split)
otu_valid <- assessment(valid_split)
x_test <- otu_valid[setdiff(names(otu_valid), "Day")]
y_test <- otu_valid$Day

rf_oob_comp <- randomForest(
  formula = Day~ .,
  data    = otu_train,
  xtest   = x_test,
  ytest   = y_test,
  ntree =500
)

oob <- sqrt(m1$mse)
validation <- sqrt(rf_oob_comp$mse)

# compare error rates
tibble::tibble(
  `Out of Bag Error` = oob,
  `Test error` = validation,
  ntrees = 1:rf_oob_comp$ntree
) %>%
  gather(Metric, Error, -ntrees) %>%
  ggplot(aes(ntrees, Error, color = Metric)) +
  geom_line() +
  xlab("Number of trees")

rf_oob_comp
sqrt(m1$mse[which.min(rf_oob_comp$mse)])

# Same process as before
# But, using dorsal sites only
otu_l4 <- otu %>%  filter(Location == 'L4')
otu_l5 <- otu %>%  filter(Location == 'L5')

otu_l4.5 <- rbind(otu_l4, otu_l5)

m1 <- randomForest(
  formula = Day ~ .,
  data    = otu_l4.5,
  ntree= 500
)

m1
sqrt(m1$mse[which.min(m1$mse)])

valid_split <- initial_split(otu_l4.5, .8)
otu_train <- analysis(valid_split)
otu_valid <- assessment(valid_split)
x_test <- otu_valid[setdiff(names(otu_valid), "Day")]
y_test <- otu_valid$Day

rf_oob_comp <- randomForest(
  formula = Day~ .,
  data    = otu_train,
  xtest   = x_test,
  ytest   = y_test,
  ntree =500
)

oob <- sqrt(m1$mse)
validation <- sqrt(rf_oob_comp$mse)

# compare error rates
tibble::tibble(
  `Out of Bag Error` = oob,
  `Test error` = validation,
  ntrees = 1:rf_oob_comp$ntree
) %>%
  gather(Metric, Error, -ntrees) %>%
  ggplot(aes(ntrees, Error, color = Metric)) +
  geom_line() +
  xlab("Number of trees")

rf_oob_comp
sqrt(m1$mse[which.min(rf_oob_comp$mse)])

# Table was organized and created in excel

## Table 5
# Quadratic regression summary of for relative abundances of three indicator taxa. 

# Model comparison code is available in supplmental tables

#looked at phylum level, changes over time in relative abundance
#picked out significant phylum found in random forest model
#built logisitcal model using phyla

otu <- as.data.frame(t(otu_table(physeq)))
otu$SampleID <- rownames(otu)
meta_sa <- metadata %>% select(SampleID, Day, Location)
otu <- merge(meta_sa, otu, by = 'SampleID')
otu <- otu[,-1]
names(otu) <- make.names(names(otu))

m1 <- randomForest(
  formula = Day ~ .,
  data    = otu,
  ntree= 500, importance = TRUE
)

m1

# extracting important taxa identified by RF
imp <- importance(m1, type=2)
imp <- data.frame(predictors = rownames(imp), imp)
imp.sort <- arrange(imp, desc(IncNodePurity))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
#found top 20 taxa
imp.20 <- imp.sort[1:20, ]
imp.20$predictors <- sub('X', '', imp.20$predictors)
colnames(imp.20) <- c('OTUID', 'IncNodePurity')
imp.20 <- merge(imp.20, tax, by = 'OTUID')
#most represented taxa were from Firmicutes, Bacteriodetes, and Proteobacteria

#phylum level taxa
physeq_phyla <- tax_glom(physeq, taxrank = 'Phylum')
#removing singletons and phyla not present in 10% of samples
physeq_phyla <- filter_taxa(physeq_phyla, function(x) sum(x > 1) > (0.10*length(x)), TRUE)
#calculating relative abundances 
physeq_phyla_rel <- transform_sample_counts(physeq_phyla, function(OTU) OTU/sum(OTU) )

#pull out important taxa
physeq_lr1 <- subset_taxa(physeq_phyla_rel, Phylum=='D_1__Firmicutes')
physeq_lr2 <- subset_taxa(physeq_phyla_rel, Phylum=='D_1__Bacteroidetes')
physeq_lr3 <- subset_taxa(physeq_phyla_rel, Phylum=='D_1__Proteobacteria')
#build phyloseq object with important taxa only
physeq_lr <- merge_phyloseq(physeq_lr1, physeq_lr2, physeq_lr3)
physeq_lr

#get data frame from phyloseq object
linreg <- as.data.frame(t(otu_table(physeq_lr)))
linreg$SampleID <- rownames(linreg)
colnames(linreg) <- c('Firmicutes', 'Bacteroidetes', 'Proteobacteria', 'SampleID')
#merge metadata with rel. abundance information 
linreg_f <- merge(linreg, metadata, by = 'SampleID')

#potential models
m0 <- lm(Day ~ 1, data = linreg_f)
m1 <- lm(Day ~ poly(Firmicutes + Bacteroidetes + Proteobacteria), data = linreg_f)

#see model comparison in supplemental code
summary(m1)

#identify and remove outliers 
cooks <- cooks.distance(m1)
plot(cooks, pch="*", cex=2) 
abline(h = 4*mean(cooks, na.rm=T), col="red")
text(x=1:length(cooks)+1, y=cooks, labels=ifelse(cooks>4*mean(cooks, na.rm=T),names(cooks),""), col="red")
influential <- as.numeric(names(cooks)[(cooks > 4*mean(cooks, na.rm=T))])
head(linreg_f[influential, ])  

linreg_test <- linreg_f%>%
  filter(SampleID != c('FP111'))

m0 <- lm(Day ~ 1, data = linreg_test)
m1 <- lm(Day ~ poly(Firmicutes + Bacteroidetes + Proteobacteria), data = linreg_test)

#model comparison using Chisq and AIC
anova(m0,m1)
AIC(m0,m1)
