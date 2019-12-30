## Cold Case Kaszubinski et al. 
# Figures

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

## Figure 1
# Figure 1 was made in powerpoint

## Figure 2
# Microbial community metrics among days. 
#A) Alpha-diversity metrics among days including: observed richness (Observed), Chao1, 
# Shannon diversity, and Inverse Simpson diversity (InvSimpson). 
#B) Principal Coordinate Analysis (PCoA) of Unifrac distances among days. 

# A
#calculating alpha div metrics
erich <- estimate_richness(physeq, measures = c("Observed", 'Chao1', "Shannon", "InvSimpson"))
erich$SampleID <- rownames(erich)
# tidy data and merge additional metadata from samples
erich <- erich %>%
  gather(Index, Observation, c("Observed", 'Chao1', "Shannon", "InvSimpson"), na.rm = TRUE)
rich = merge(erich, metadata)
#make sure day is a factor rather than a number for plotting purposes
rich$Day <- factor(rich$Day, levels = c('1', '5', '9', '13', '17', '22'))

#make alpha div component of figure 3
p <- ggplot(rich, aes(x=Day, y=Observation, fill=Day)) +
  geom_boxplot() + 
  ylab('Alpha-Diversity Metric') +
  scale_fill_manual(values = c('#2F6E2B', '#E37D3E', '#886F98', 
                               '#E4BE18', '#A41400', '#0C677C')) + 
  facet_wrap(~Index, scales="free") 
p

# B
#calculate beta div with unifrac dissimilarity matrix
ord = ordinate(physeq_beta, method="PCoA", distance="unifrac")
#create plot
ordplot=plot_ordination(physeq_beta, ord, color="Day", shape = 'Swab_Areas') +
  scale_color_manual(values = c('#2F6E2B', '#E37D3E', '#886F98', '#E4BE18', '#A41400', '#0C677C')) +
  scale_shape_manual(values=c(17,8,16,15,18)) + geom_point(size = 5) + 
  labs(shape = "Sampling Location")

ordplot

#make final figure and export as a tiff 
theme_set(theme_classic(base_size = 18))
tiff("FIG2.TIF", width = 4500, height = 4000, res=300)
ggarrange(p,ordplot, 
          labels = c("A", "B"),
          nrow = 2, ncol = 2)
dev.off()

## Figure 3
# Relative abundances of differentially abundant taxa among day indicated by ANCOM. 
# A) Differentially abundant taxa among day 1, and days 5-22. 
# B) Differentially abundant taxa among day 5, and days 1 and 9-22. 
# C) Differentially abundant taxon among day 9, and days 1-5, and 13-22. 

# Differentially abundant taxa identified by ANCOM were determined. See supplemental
# tables for code

#calculate relative abundance
#removing singletons and phyla not present in 10% of samples
physeq_phy <- tax_glom(physeq, taxrank = 'Phylum')
physeq_phyla <- filter_taxa(physeq_phy, function(x) sum(x > 1) > (0.10*length(x)), TRUE)
physeq_phyla_rel <- transform_sample_counts(physeq_phyla, function(OTU) OTU/sum(OTU) )


# Relative abundances of phyla by day 
df <- psmelt(physeq_phyla_rel) 
Trtdata <- ddply(df, c("Phylum", 'Day'), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
Trtdata$Phylum <- as.character(Trtdata$Phylum)
#targets were identified by ANCOM 
targets <- c('D_1__Bacteroidetes', 'D_1__Chlamydiae', ' D_1__Cyanobacteria',
             'D_1__Dependentiae', 'D_1__Elusimicrobia', 'D_1__Fibrobacteres', 
             'D_1__Latescibacteria', 'D_1__Margulisbacteria', ' D_1__Tenericutes', 
             ' D_1__Thaumarchaeota')
#filter out non significant phyla
Trtdata_one <- filter(Trtdata, Phylum %in% targets)  

# A
adataplot=ggplot(Trtdata_one, aes(x=Day,y=mean))+
  geom_bar(aes(fill = Phylum),colour="black", stat="identity") + 
  scale_fill_manual(values = c('#955DB7', '#D266C8','#ED402F' , '#2FED6D',
                               '#ED9C2F', '#2FA4ED', '#EDC72F'
  )) + ylab('Relative Abundance of Phyla')
adataplot

# B
targets <- c('D_1__Latescibacteria',
             'D_1__Margulisbacteria', 
             'D_1__Bacteroidetes')

Trtdata_five <- filter(Trtdata, Phylum %in% targets)

bdataplot=ggplot(Trtdata_five, aes(x=Day,y=mean))+
  geom_bar(aes(fill = Phylum),colour="black", stat="identity") +
  scale_fill_manual(values=c('#955DB7', '#2FA4ED', '#EDC72F')) + 
  ylab('Relative Abundance of Phyla')
bdataplot

# C
targets <- c('D_1__Margulisbacteria')

Trtdata_nine <- filter(Trtdata, Phylum %in% targets)

cdataplot=ggplot(Trtdata_nine, aes(x=Day,y=mean))+
  geom_bar(aes(fill = Phylum),colour="black", stat="identity") +scale_fill_manual(values = c('#EDC72F')) + 
  ylab('Relative Abundance of Phyla')
cdataplot

# Make final figure and export as tiff
theme_set(theme_classic(base_size = 18))
tiff("FIG3.TIF", width = 2000, height = 4500, res=300)
ggarrange(adataplot, bdataplot, cdataplot, 
          labels = c("A", "B", 'C'),
          nrow = 3, ncol = 1)
dev.off()


## Figure 4
# Microbial community metrics among sampling locations (DC: dorsal cloth, DS: dorsal skin, M: mouth; VC: ventral cloth, VS: ventral skin).
# A) Alpha-diversity metrics among sampling locations including: observed richness (Observed), Chao1, Shannon diversity, and Inverse Simpson diversity (InvSimpson). 
# B) Principal Coordinate Analysis (PCoA) of Unifrac distances among sampling locations. 

erich <- estimate_richness(physeq, measures = c("Observed", 'Chao1', "Shannon", "InvSimpson"))
erich <- add_rownames(erich, "SampleID")

# get data organized
erich <- erich %>%
  gather(Index, Observation, c("Observed", 'Chao1', "Shannon", "InvSimpson"), na.rm = TRUE)
rich = merge(erich, metadata)

theme_set(theme_bw(base_size = 18))

p <- ggplot(rich, aes(x=Swab_Areas, y=Observation, fill=Swab_Areas)) +
  geom_boxplot() + xlab('Sampling Location') + ylab('Alpha-diversity Metrics') + 
  facet_wrap(~Index, scales="free") + labs(fill = "Sampling Location") +
  scale_x_discrete(labels=c('mouth' = 'M', 'ventral_cloth' = 'VC',
                            'ventral_skin' = 'VS', 'dorsal_cloth' = 'DC',
                            'dorsal_skin' = 'DS')) +
  scale_fill_manual(values = c('#AA8846', '#ECBE64', '#50CF9A', '#317585', '#50B7CF'))
p

ord = ordinate(physeq_beta, method="PCoA", distance="unifrac")
ordplot=plot_ordination(physeq_beta, ord, color="Swab_Areas")
ordplot

p2 <- ordplot +
  geom_point(size = 4) + stat_ellipse(alpha = 0.03, geom = "polygon", aes(fill = Swab_Areas)) +
  theme(legend.position = 'none')  +
  scale_color_manual(values = c('#AA8846', '#ECBE64', '#50CF9A', '#317585', '#50B7CF'))
p2

theme_set(theme_classic(base_size = 18))
tiff("FIG4.TIF", width = 4500, height = 3000, res=300)
ggarrange(p,p2, 
          labels = c("A", "B"),
          nrow = 2, ncol = 2)
dev.off()

## Figure 5
# Error of the best performing random forest regression model for PMSI. 
# True values were regressed against predicted values. 
# Error bars were made based on the mean square error and the linear regression model equation was reported.  

otu <- as.data.frame(t(otu_table(physeq)))
otu$SampleID <- rownames(otu)
meta_sa <- metadata %>% select(SampleID, Day, Location)
otu <- merge(meta_sa, otu, by = 'SampleID')
otu <- otu[,-1]
names(otu) <- make.names(names(otu))

m1 <- randomForest(
  formula = Day ~ .,
  data    = otu,
  ntree= 500
)

m1
pred <- m1$predicted
pred <- as.data.frame(pred)
colnames(pred) <- 'Predicted'

#using the mean square error to determine the true value
pred[pred$Predicted < 24.858, "True"] <- '22'
pred[pred$Predicted < 19.858, "True"] <- '17'
pred[pred$Predicted < 15.858, "True"] <- '13'
pred[pred$Predicted < 11.858, "True"] <- '9'
pred[pred$Predicted < 7.858, "True"] <- '5'
pred[pred$Predicted < 3.858, "True"] <- '1'

pred$True <- as.numeric(pred$True)

ggplotRegression <- function(fit){
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_abline(slope =0, intercept = 0, color='black') +
    geom_vline(xintercept = 0, color='black') +
    stat_smooth(method = "lm", col = "red") 
}

d <- ggplotRegression(lm(Predicted ~ True, data = pred)) +
  geom_pointrange(aes(ymin = Predicted-2.858, ymax = Predicted+2.858), 
                  position=position_jitter(width=0.5), alpha = 0.5) 

lm_eqn <- function(df){
  m <- lm(Predicted ~ True, df);
  eq <- substitute(italic(Predicted) == a + b %.% italic(True)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

d1 <- d + geom_text(x = 5, y = 23, label = lm_eqn(pred), parse = TRUE)
d1

theme_set(theme_classic(base_size = 18))
tiff("FIG5.TIF", width = 3000, height = 3000, res=300)
d1
dev.off()


## Figure 6
# Quadratic regression of PMSI indicator taxa. 
# A) Relative abundance of Firmicutes regressed across days. 
# B) Relative abundance of Bacteroidetes regressed across days. 
# C) Relative abundance of Proteobacteria regressed across days.    

# See Tables code for model comparisons identifying Firmicutes, Bacteroidetes, and
# Proteobacteria. 

#subset phyloseq object by important taxa
physeq_lr1 <- subset_taxa(physeq_phyla_rel, Phylum=='D_1__Firmicutes')
physeq_lr2 <- subset_taxa(physeq_phyla_rel, Phylum=='D_1__Bacteroidetes')
physeq_lr3 <- subset_taxa(physeq_phyla_rel, Phylum=='D_1__Proteobacteria')

#merge them together
physeq_lr <- merge_phyloseq(physeq_lr1, physeq_lr2, physeq_lr3)

#take out OTU table from object for plotting purposes and adding SampleIDS
linreg <- as.data.frame(t(otu_table(physeq_lr)))
linreg$SampleID <- rownames(linreg)
colnames(linreg) <- c('Firmicutes', 'Bacteroidetes', 'Proteobacteria', 'SampleID')
# add metadata to relative abundance information 
linreg_f <- merge(linreg, metadata, by = 'SampleID')

#separate plots for each taxa
a <- ggplot(linreg_f, aes(y=Firmicutes, x=Day)) +
  geom_point(color = '#3E8E18', size=4) + stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  ylab('Relative Abundance of Firmicutes')
a

b <- ggplot(linreg_f, aes(y=Bacteroidetes, x=Day)) +
  geom_point(color= '#1C5AA6', size = 4) + stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  ylab('Relative Abundance of Bacteroidetes')
b

c <- ggplot(linreg_f, aes(y=Proteobacteria, x=Day)) +
  geom_point(color = '#7732BF', size = 4) + stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  ylab('Relative Abundance of Proteobacteria')
c

#make final figure
theme_set(theme_classic(base_size = 10))
tiff("FIG6.TIF", width = 3000, height = 1000, res=300)
ggarrange(a,b,c, 
          labels = c("A", "B", "C"),
          nrow = 1, ncol = 3)
dev.off()