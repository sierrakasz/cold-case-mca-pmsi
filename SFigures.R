## Cold Case Kaszubinski et al. 
# Supplemental Figures

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


## S Figure 1
# Shannon diversity across sequencing depth among microbial samples. 
# Library size of 5,000 was chosen based on the plateau of 
# Shannon diversity metrics at larger sequencing depths. 


df <- read.csv('alpharare-pigs.csv')

#take first column and make it the rownames
rownames(df) <- df[,1]
#remove first column
df_new <- df[,-1]
#take out extra metadata
df_new <- df_new[,-c(101:112)]
#remove delware samples
df_new <- df_new[-c(85:118),]

#transpose dataframe
df_t <- t(df_new)
df_t <- as.data.frame(df_t)

#adding the sequencing depth as a column
df_t$depth <- rownames(df_t)

#make data tidy
boxplot <- df_t %>% 
  gather(key = "SampleID", value = "shannon_div")
#remove depth as samples
boxplot <- boxplot %>%  filter(SampleID != 'depth')

#re-add depth column from colnames of old dataframe
depth <- colnames(df_new)
depth <- c(rep(depth, 84))
boxplot$depth <- depth
#remove iteration from name
boxplot$depth = substr(boxplot$depth,1,nchar(boxplot$depth)-7)
boxplot$depth <- gsub( "_", "", as.character(boxplot$depth))

#not going to use all levels, just showing plateo 
boxplot$depth <- factor(boxplot$depth, levels= c('depth.1', 'depth.5556', 'depth.11111', 'depth.16667', 'depth.22222'))
boxplot <- boxplot[complete.cases(boxplot),]

#make sure Shannon div is numeric not a factor
boxplot$shannon_div <- as.numeric(boxplot$shannon_div)

#make figure
theme_set(theme_bw(base_size = 18))
tiff("SFIG1.TIF", width = 4000, height = 3000, res=300)
ggplot(boxplot, aes(x = depth, y = shannon_div)) + geom_boxplot(fill = '#F1AC37') + xlab('Sequencing Depth') +
  ylab('Shannon Diversity') + 
  stat_summary(fun.y=mean, geom="line", aes(group=1))
dev.off()

## S Figure 2

# Hourly air temperature for the Grand Rapid and Lansing areas during the 2005 and 2018 
# time periods of interest, respectively.  A t-test was used to determine significant 
# differences in mean air temperatures between the time periods of interest. 
# Each point represents the average within an hour. 

air_temp <- read.csv("temp_data_coldcase_MSU.csv")
colnames(air_temp) <- c('Date', 'Temp', 'Location')
air_temp <- air_temp[complete.cases(air_temp),]

df_air_temp <- ddply(air_temp, c('Date', 'Location'), summarise,
                     mean_temp = mean(Temp),
                     sd_temp   = sd(Temp)
)

shapiro.test(df_air_temp$mean_temp)
a <- aov(mean_temp ~ Location, df_air_temp)
summary(a)


theme_set(theme_bw(base_size = 18))
tiff("SFIG2.TIF", width = 4000, height = 3000, res=300)
ggplot(df_air_temp, aes(x = Date, y = mean_temp, col = Location, group=22)) + 
  geom_point(size = 4) +
  geom_errorbar(width=.5, aes(ymin=mean_temp-sd_temp, ymax=mean_temp+sd_temp), data=df_air_temp) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_color_manual(values = c('#6127E5', '#097629')) + geom_line() + 
  ylab('Mean Air Temperature in Celsius (SD)') +
  facet_wrap(~Location, scales = "free_x")
dev.off()


