# change working directory
setwd("C:/Yuchen/WCM/Courses/2023Spring/CMPB5005/homework/FIHTD/1 - Metabolomics/")
load("Simulated_metabolomics_data.RData")
library(Rcpm)
library(ggplot2)

#log data, transpose matrix and covert to dataframe for boxplot function
logged_dat <-data.frame(t(log(dat)))
#rename colnames as 1-906
colnames(logged_dat)<-seq(1:906)
#create boxplot, dont plot outliers
samples_boxplot <-boxplot(logged_dat, xlab = "Sample", ylab = "log(concentration)" , outline = FALSE)


#quotient normalization using Dieterle, F., Ross, A., Schlotterbeck, G. & Senn, H. Probabilistic Quotient Normalization as Robust Method to Account for Dilution of Complex Biological Mixtures. Application in H1 NMR Metabonomics. Anal. Chem. 78, 4281-4290 (2006).
normalized_dat <-pqn(dat,QC=NULL)
#log data, transpose matrix and covert to dataframe for boxplot function
norm_logged_dat <-data.frame(t(log(normalized_dat)))
#rename colnames as 1-906
colnames(norm_logged_dat)<-seq(1:906)
#create boxplot, dont plot outliers
norm_samples_boxplot <-boxplot(norm_logged_dat, xlab = "Sample", ylab = "log(concentration)" , outline = FALSE)


#get indices of males and females
males <- which(gender == 1)
females <-which(gender ==2)
#subset normalized matrix by male and female
male_metabolites_norm_dat <- normalized_dat[males,]
female_metabolites_norm_dat <-normalized_dat[females,]
#vector to hold p values for each metabolite
pvals <- c()
log2fold <- c()
#for each metabolite, do t test comparing concentation between males and females
for (i in 1:length(colnames(normalized_dat))){
  male_sample <- male_metabolites_norm_dat[,i]
  female_sample <- female_metabolites_norm_dat[,i]
  temp <- t.test(male_sample,female_sample)
  pvals[i] <- temp$p.value
  log2fold[i] <- log2(mean(male_sample))-log2(mean(female_sample))
}

# create dataframe for volcano plot, contains log2fold change, p-value, metabolites name, and labels for differentially expressed metabolites
de <- data.frame(log2fold,pvals,annotations$name)
names(de) <- c('log2FoldChange','pvalue','name')
# categorize whether metabolites are expressed differently in M or F
de$diffexp <- 'No Diff'
de$diffexp[de$pvalue<0.05 & de$log2FoldChange>0] <- 'High in Male'
de$diffexp[de$pvalue<0.05 & de$log2FoldChange<0] <- 'High in Female'
# create labels for highly differentially expressed metabolites
de$delabel <- NA
de$delabel[de$diffexp!='No Diff' & abs(de$log2FoldChange)>1] <- de$name[de$diffexp!='No Diff' & abs(de$log2FoldChange)>1]

ggplot(de, aes(x=pvalue)) + geom_histogram(bins=50)

ggplot(data=de, aes(x=log2FoldChange,y=-log10(pvalue),col=diffexp,label=delabel)) + geom_point() + theme_minimal() + geom_text(check_overlap=T,size=3,vjust=0,nudge_y=0.25)


#copy DE dataframe and keep only metabolite name, logFC, p value
de_adj <-de
de_adj <- de_adj[,c(1,2,3)]
#add adjusted p value using FDR
de_adj$p.adj <-p.adjust(de_adj$pvalue, method = "fdr")

#Relabel metabolites after correction
de_adj$diffexp <- 'No Diff'
de_adj$diffexp[de_adj$p.adj<0.05 & de_adj$log2FoldChange>0] <- 'High in Male'
de_adj$diffexp[de_adj$p.adj<0.05 & de_adj$log2FoldChange<0] <- 'High in Female'
# create labels for highly differentially expressed metabolites
de_adj$delabel <- NA
de_adj$delabel[de_adj$diffexp!='No Diff' & abs(de_adj$log2FoldChange)>1] <- de_adj$name[de_adj$diffexp!='No Diff' & abs(de_adj$log2FoldChange)>1]

#count significant metabolites for each group
sig_female_before = length(which(de$diffexp == "High in Female"))
sig_male_before = length(which(de$diffexp == "High in Male"))

sig_female_after = length(which(de_adj$diffexp == "High in Female"))
sig_male_after = length(which(de_adj$diffexp == "High in Male"))

#create dataframe for barplot
counts <- c(sig_female_before,sig_male_before, sig_female_after, sig_male_after)
group <- c("Sig_Female_Before", "Sig_Male_Before", "Sig_Female_After", "Sig_Male_After")
barplot_df <- data.frame(group,counts )

#barplot
ggplot(data=barplot_df, aes(x=group, y=counts)) +
  geom_bar(stat="identity") + scale_x_discrete(limits=c("Sig_Female_Before","Sig_Female_After","Sig_Male_Before", "Sig_Male_After"))