## Cold Case Kaszubinski et al. 
# Supplemental Tables

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
library(outliers)
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

## S Table 1
air_temp <- read.csv("temp_data_coldcase_MSU.csv")
colnames(air_temp) <- c('Date', 'Temp', 'Location')
air_temp <- air_temp[complete.cases(air_temp),]
df_air_temp <- ddply(air_temp, c('Location'), summarise,
                     mean_temp = mean(Temp),
                     sd_temp   = sd(Temp),
                     min_temp = min(Temp),
                     max_temp = max(Temp)
)

air_temp_gr <- air_temp %>% filter(Location == 'Grand Rapids')
air_temp_l <- air_temp %>% filter(Location != 'Grand Rapids')

grubbs.test(air_temp_gr$Temp, type = 11)
grubbs.test(air_temp_l$Temp, type = 11)

## S Table 2
# See above

## S Table 3
# Table was built in excel 

## S Table 4
# Table was built in excel

## S Table 5
# A) Kruskal-Wallis test among alpha diversity metrics across days. 
# B) Mean and standard deviation of alpha diversity metrics across days. 
# C) Pair-wise posthoc Nemenyi test

# A
erich <- estimate_richness(physeq, measures = c("Observed", 'Chao1', "Shannon", "InvSimpson"))
erich <- add_rownames(erich, "SampleID")
erich_sums <- merge(erich, metadata)
erich_sums %>% group_by(Day) %>% summarise_at(c('Observed', 'Chao1', "Shannon", "InvSimpson"), funs(mean, sd))

# B and C

# get data organized and formatted for stats testing
erich <- erich %>%
  gather(Index, Observation, c("Observed", 'Chao1', "Shannon", "InvSimpson"), na.rm = TRUE)
rich = merge(erich, metadata)

prepare_samples_kw <- function(rich) {
  return(list(rich_obs <- rich %>% filter(Index == 'Observed'),
              rich_cha <- rich %>% filter(Index == 'Chao1'),
              rich_sha <- rich %>% filter(Index == 'Shannon'),
              rich_inv <- rich %>% filter(Index == 'InvSimpson')))
}

kw_values <- prepare_samples_kw(rich)

# KW and post-hoc Nemenyi 
for(i in 1:length(kw_values)) {
  print(kruskal.test(Observation ~ Day, data = kw_values[[i]]))
  out <- posthoc.kruskal.nemenyi.test(x=kw_values[[i]]$Observation, g=kw_values[[i]]$Day, dist='Tukey', p.adjust.method = 'bonf' )
  print(out$p.value)
}

## S Table 6
# A) PERMANOVA results for comparing beta-diversity and beta-dispersion among days for 999 permuations. 
# B) Pair-wise PERMANOVA results for comparing beta-diversity and beta-dispersion among days for 999 permuations. 

# A
beta_diversity_calc_day <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "unifrac")
  sampledf <- data.frame(sample_data(physeq))
  print(return(adonis(GPdist ~ Day, data = sampledf)))
}

beta_dispersion_calc_day <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "unifrac")
  sampledf <- data.frame(sample_data(physeq))
  beta <- betadisper(GPdist, sampledf$Day)
  print(return(permutest(beta)))
}

beta_diversity_calc_day(physeq_beta)
beta_dispersion_calc_day(physeq_beta)

# B
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

day_list <- day_pairwise(physeq_beta)

for(i in 1:length(day_list)) {
  print(beta_diversity_calc_day(day_list[[i]]))
  print(beta_dispersion_calc_day(day_list[[i]]))
}



## S Table 7
# Analysis of the composition of microbiomes (ANCOM) pairwise comparisons for differentially abundant phyla among days. 
# Only significant results at a 0.7 detection level or higher were included.  

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

# examining phylum level, reduce phyloseq object to phylum
physeq_phy <- tax_glom(physeq, taxrank = 'Phylum')

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

#create list of all day pairwise comparisons
list_for_ancom_phy <- day_pairwise(physeq_phy)

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

#blank list to collect results from ANCOM function in
ancom_results_phy <- list()

# run ANCOM for all day pairwise comparisons
for(i in 1:length(list_for_ancom_phy)) {
  otu_ancom <- otu_ancom_make(list_for_ancom_phy[[i]])
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
  tax <- tax %>% select(Kingdom, Phylum)
  tax$otu.names <- rownames(tax)
  ancom_sign_taxa <- merge(w_values, tax, by='otu.names')
  ancom_sign_taxa <- ancom_sign_taxa[,-1]
  ancom_results_phy[[i]] <- ancom_sign_taxa
}

## S Table 8
# nalysis of the composition of microbiomes (ANCOM) pairwise comparisons for differentially abundant families among days. 
# Only significant results at a 0.7 detection level or higher were included.  

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

## S Table 9
# Alpha diversity metrics among pig replicates. Averge and standard deviation for each 
# metric are included, as well as the chi-squared value, degrees of freedom, 
# and p-value for Kruskal-Wallis test comparing among pig replicates. 

erich <- estimate_richness(physeq, measures = c("Observed", 'Chao1', "Shannon", "InvSimpson"))
erich <- add_rownames(erich, "SampleID")

erich <- erich %>%
  gather(Index, Observation, c("Observed", 'Chao1', "Shannon", "InvSimpson"), na.rm = TRUE)
erich %>% group_by(Index) %>% summarise_all(funs(mean, sd))

rich = merge(erich, metadata)


prepare_samples_kw <- function(rich) {
  return(list(rich_obs <- rich %>% filter(Index == 'Observed'),
              rich_cha <- rich %>% filter(Index == 'Chao1'),
              rich_sha <- rich %>% filter(Index == 'Shannon'),
              rich_inv <- rich %>% filter(Index == 'InvSimpson')))
}

kw_values <- prepare_samples_kw(rich)

for(i in 1:length(kw_values)) {
  print(kruskal.test(Observation ~ Pig, data = kw_values[[i]]))
  out <- posthoc.kruskal.nemenyi.test(x=kw_values[[i]]$Observation, g=kw_values[[i]]$Pig, dist='Tukey', p.adjust.method = 'bonf' )
  print(out$p.value)
}

## S Table 10
# A) Mean and standard deviation of Unifrac beta-diversity values. 
# B) PERMANOVA results for comparing beta-diversity among pig replicates for 999 permuations. 
# C) Mean and standard deviation of Unifrac beta-dispersion values. 
# D) PERMANOVA results for comparing beta-dispersion among pig replicates for 999 permuations. 

GPdist=phyloseq::distance(physeq_beta, "unifrac")
sampledf <- data.frame(sample_data(physeq_beta))
a <- adonis(GPdist ~ Pig, data = sampledf)
beta_dist <- a$coef.sites
beta <- betadisper(GPdist, sampledf$Pig)
permutest(beta)

## S Table 11
# Random forest classification of pig (replicates 1-5) using out-of-bag (OOB) classification and 1,000 decision trees. 
# Confusion matrix shows classification of pig replicate and overall error rate.  

random_foresting_pig <- function(data_phy) {
  fd = data_phy
  predictors=t(otu_table(fd))
  dim(predictors)
  resp <- as.factor(sample_data(fd)$Pig)
  rf.data <- data.frame(resp, predictors)
  Forest <- randomForest(resp~., data=rf.data, ntree=1000)
  return(Forest)
}

m1 <- random_foresting_pig(physeq)
plot(m1)

## S Table 12
# Random forest classification of sampling location using out-of-bag (OOB) classification and 1,000 decision trees. 
# Confusion matrix shows classification of sampling location and overall error rate.  

otu <- as.data.frame(t(otu_table(physeq)))
otu$SampleID <- rownames(otu)
meta_sa <- metadata %>% select(SampleID, Location)
otu <- merge(meta_sa, otu, by = 'SampleID')
otu <- otu[,-1]
names(otu) <- make.names(names(otu))

m1 <- randomForest(
  formula = Location ~ .,
  data    = otu,
  ntree= 1000
)

m1
plot(m1)

## S Table 13
# A) Alpha diversity metrics among sampling locations. Averge and standard deviation for each metric are included, 
# as well as the chi-squared value, degrees of freedom, and p-value for Kruskal-Wallis test comparing among pig replicates. 
# B) Post-hoc pairwise Nemenyi test among sampling locations and alpha diversity metrics. 

erich <- estimate_richness(physeq, measures = c("Observed", 'Chao1', "Shannon", "InvSimpson"))
erich <- add_rownames(erich, "SampleID")

erich <- erich %>%
  gather(Index, Observation, c("Observed", 'Chao1', "Shannon", "InvSimpson"), na.rm = TRUE)
erich %>% group_by(Index) %>% summarise_all(funs(mean, sd))

rich = merge(erich, metadata)


prepare_samples_kw <- function(rich) {
  return(list(rich_obs <- rich %>% filter(Index == 'Observed'),
              rich_cha <- rich %>% filter(Index == 'Chao1'),
              rich_sha <- rich %>% filter(Index == 'Shannon'),
              rich_inv <- rich %>% filter(Index == 'InvSimpson')))
}

kw_values <- prepare_samples_kw(rich)
for(i in 1:length(kw_values)) {
  kw_values[[i]]$Location_simple <- 
    ifelse(kw_values[[i]]$Swab_Areas == 'ventral_skin', 'ventral',
           ifelse(kw_values[[i]]$Swab_Areas == 'ventral_cloth', 'ventral',
                  ifelse(kw_values[[i]]$Swab_Areas == 'dorsal_skin', 'dorsal',
                         ifelse(kw_values[[i]]$Swab_Areas == 'dorsal_cloth', 'dorsal',
                                ifelse(kw_values[[i]]$Swab_Areas == 'mouth', 'mouth', 999)))))
  kw_values[[i]]$Location_simple <- as.factor(kw_values[[i]]$Location_simple) 
}


for(i in 1:length(kw_values)) {
  print(kruskal.test(Observation ~ Location, data = kw_values[[i]]))
  out <- posthoc.kruskal.nemenyi.test(x=kw_values[[i]]$Observation, g=kw_values[[i]]$Location, dist='Tukey', 
                                      p.adjust.method = 'bonf' )
  print(out$p.value)
}


for(i in 1:length(kw_values)) {
  print(kruskal.test(Observation ~ Location_simple, data = kw_values[[i]]))
  out <- posthoc.kruskal.nemenyi.test(x=kw_values[[i]]$Observation, g=kw_values[[i]]$Location_simple, dist='Tukey',
                                      p.adjust.method = 'bonf' )
  print(out$p.value)
}

## S Table 14
# A) PERMANOVA results for comparing beta-diversity and beta-dispersion among sample locations for 999 permuations. 
# B) Pair-wise PERMANOVA results for comparing beta-diversity among sampling locations for 999 permuations. 

beta_diversity_calc_loc <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "unifrac")
  sampledf <- data.frame(sample_data(physeq))
  print(return(adonis(GPdist ~ Location, data = sampledf)))
}

beta_dispersion_calc_loc <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "unifrac")
  sampledf <- data.frame(sample_data(physeq))
  beta <- betadisper(GPdist, sampledf$Location)
  print(return(permutest(beta)))
}

beta_diversity_calc_loc(physeq_beta)
beta_dispersion_calc_loc(physeq_beta)

location_pairwise <- function(physeq) {
  physeq_L1 <- subset_samples(physeq, Location == 'L1')
  physeq_L2 <- subset_samples(physeq, Location == 'L2')
  physeq_L3 <- subset_samples(physeq, Location == 'L3')
  physeq_L4 <- subset_samples(physeq, Location == 'L4')
  physeq_L5 <- subset_samples(physeq, Location == 'L5')
  return(list(physeq_L12 <- merge_phyloseq(physeq_L1, physeq_L2),
              physeq_L13 <- merge_phyloseq(physeq_L1, physeq_L3),
              physeq_L14 <- merge_phyloseq(physeq_L1, physeq_L4),
              physeq_L15 <- merge_phyloseq(physeq_L1, physeq_L5),
              physeq_L23 <- merge_phyloseq(physeq_L3, physeq_L2),
              physeq_L24 <- merge_phyloseq(physeq_L4, physeq_L2),
              physeq_L25 <- merge_phyloseq(physeq_L5, physeq_L2),
              physeq_L34 <- merge_phyloseq(physeq_L4, physeq_L3),
              physeq_L35 <- merge_phyloseq(physeq_L3, physeq_L5),
              physeq_L45 <- merge_phyloseq(physeq_L4, physeq_L5)))
}

location_list <- location_pairwise(physeq_beta)

for(i in 1:length(location_list)) {
  print(beta_diversity_calc_loc(location_list[[i]]))
}

## S Table 15
# Random forest regression of ADH using out-of-bag (OOB) classification and 500 decision trees.
# Variation explained, mean square residuals, and mean squre error are reported. 

otu <- as.data.frame(t(otu_table(physeq)))
otu$SampleID <- rownames(otu)
meta_sa <- metadata %>% select(SampleID, ADH)
otu <- merge(meta_sa, otu, by = 'SampleID')
otu <- otu[,-1]
names(otu) <- make.names(names(otu))

m1 <- randomForest(
  formula = ADH ~ .,
  data    = otu,
  ntree= 500
)

m1
sqrt(m1$mse[which.min(m1$mse)])

valid_split <- initial_split(otu, .8)
otu_train <- analysis(valid_split)
otu_valid <- assessment(valid_split)
x_test <- otu_valid[setdiff(names(otu_valid), "ADH")]
y_test <- otu_valid$ADH

rf_oob_comp <- randomForest(
  formula = ADH~ .,
  data    = otu_train,
  xtest   = x_test,
  ytest   = y_test,
  ntree = 500
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


## Model Comparison for linear regression 
physeq_phyla <- tax_glom(physeq, taxrank = 'Phylum')

#removing singletons and phyla not present in 10% of samples
physeq_phyla <- filter_taxa(physeq_phyla, function(x) sum(x > 1) > (0.10*length(x)), TRUE)

physeq_phyla_rel <- transform_sample_counts(physeq_phyla, function(OTU) OTU/sum(OTU) )

physeq_lr1 <- subset_taxa(physeq_phyla_rel, Phylum=='D_1__Firmicutes')
physeq_lr2 <- subset_taxa(physeq_phyla_rel, Phylum=='D_1__Bacteroidetes')
physeq_lr3 <- subset_taxa(physeq_phyla_rel, Phylum=='D_1__Proteobacteria')

physeq_lr <- merge_phyloseq(physeq_lr1, physeq_lr2, physeq_lr3)
physeq_lr

linreg <- as.data.frame(t(otu_table(physeq_lr)))
linreg$SampleID <- rownames(linreg)
colnames(linreg) <- c('Firmicutes', 'Bacteroidetes', 'Proteobacteria', 'SampleID')
linreg_f <- merge(linreg, metadata, by = 'SampleID')

m0 <- lm(Day ~ 1, data = linreg_f)
m1 <- lm(Day ~ Firmicutes + Bacteroidetes + Proteobacteria, data = linreg_f)
m1.1 <- lm(Day~ Firmicutes, data = linreg_f)
m1.2 <- lm(Day~ Proteobacteria, data = linreg_f)
m1.3 <- lm(Day~ Bacteroidetes, data = linreg_f)
#Firmicutes not significant. Try without it
m2 <- lm(Day ~ Bacteroidetes + Proteobacteria, data = linreg_f)
# not doing great, add temp
m3 <- lm(Day ~ Proteobacteria + Bacteroidetes + Temp, data = linreg_f)
m3.3 <- lm(Day ~ Firmicutes + Proteobacteria + Bacteroidetes + Temp, data = linreg_f)
#not much better
#try quadratic 
m4 <- lm(Day ~ poly(Firmicutes + Bacteroidetes + Proteobacteria), data = linreg_f)
#Better!
m5 <- lm(Day ~ poly(Proteobacteria + Bacteroidetes + Temp), data = linreg_f)
m6 <- lm(Day ~ poly(Proteobacteria + Bacteroidetes), data = linreg_f)
m7 <- lm(Day ~ poly(Firmicutes + Bacteroidetes + Proteobacteria + Temp), data = linreg_f)

#add models as needed to look at summary
summary(m7)

AIC(m0,m1,m1.1, m1.2, m1.3, m2,m3,m3.3, m4,m5, m6, m7)

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

anova(m0,m1)
AIC(m0,m1)


