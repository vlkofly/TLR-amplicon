# script that takes files with properties of segregating sites of TLR and output heatmap and PCA
#install.packages(c("grid","gridExtra"))
library(devtools)
library(ggplot2)
#install_github("slowkow/ggrepel")
library(ggrepel)
library(grid)
library(gridExtra)
library(gtable)
library(RColorBrewer)
library(reshape2)

#?install_github
###set variables####
tlr<-"TLR1B" # supply locus name
tlrdirectory<-"~/TLR_mockers/illumina_data/phase_TLR1B/translation/" # select directory according the locus
#galahap<-c(3,4,10,12,13,16,17,18) # haplotype numbers that occur in galapagos TLR4
#galahap<-c(1,5,6,8,10,11,13,15,18,20,22) # haplotype numbers that occur in galapagos TLR15
#hapfreq<-c(1,1,2,43,1,1,1,1,1,360,2,1,9,1,24,1,2,6) # haplotype freq TLR4 from frequency_dna_proteins_tlr4.ods
hapfreq<-c(391,1,3,1,30) # haplotype freq TLR4 from frequency_dna_proteins_tlr4.ods
#hapfreq<-c(2,2,5,2,1,3,2,59,1,170,84,1,6,1,1,4,2,2,1,2,1,50) # haplotype freq TLR4 from frequency_dna_proteins_tlr4.ods
# vector with population origin follows again from the table parsed by paste -d, ---| sed 's/,/","/g'
#pop.tlr<-c("gala","cont","cont","cont","gala","gala","cont","gala","cont","gala","gala","cont","gala","cont","gala","cont","cont","gala","cont","cont","cont","gala") #TLR15
#pop.tlr<-c("cont","cont","gala","gala","cont","cont","cont","cont","cont","gala","cont","gala","gala","cont","cont","gala","gala","gala") # TLR4
pop.tlr<-c("both","cont","gala","gala","cont")

####varibles for MHC##### probably added later so now it wont work for TLR?
if (tlr == "MHC") {
       print("analysing MHC")	
mhc_prot_table<-read.delim("prot_table.tsv") # key for DNA-protein code
mhc_prot_table<-mhc_prot_table[,-3]
mhc_prot_table[2]
dnc.tmp<-apply(mhc_prot_table[2],2,as.character)
mhc_prot_table[2]<-gsub("\\)","",gsub("c\\(","",dnc.tmp))
mhc_prot_table[2]<-apply(mhc_prot_table[2],2,as.numeric)
frequency_dna_proteins_MHC<-merge(mhc_prot_table,MHCfreqpop,by.x ="dna_code",by.y = "MHCallelename", all.x = T )
dim(frequency_dna_proteins_MHC)
#continue here
head (frequency_dna_proteins_MHC)
frequency_dna_proteins_MHC[frequency_dna_proteins_MHC$Champion == "NA",]
MHCfreqpop[,1]
# fill lines where more alleles represented just one protein haplotype
line40.tmp<-MHCfreqpop[MHCfreqpop$MHCallelename == 3214,-1] + MHCfreqpop[MHCfreqpop$MHCallelename == 1870,-1]
frequency_dna_proteins_MHC[72,-c(1,2)] <-line40.tmp
line53.tmp<-MHCfreqpop[MHCfreqpop$MHCallelename == 3008,-1] + MHCfreqpop[MHCfreqpop$MHCallelename == 5,-1]
frequency_dna_proteins_MHC[53,-c(1,2)] <-line53.tmp
frequency_dna_proteins_MHC<-frequency_dna_proteins_MHC[order(frequency_dna_proteins_MHC$prot_code),] # sort the table so protein names are ascending
#get frequency of each protein haplotype
rownames(frequency_dna_proteins_MHC) <- frequency_dna_proteins_MHC$prot_code
hapfreq<-apply(frequency_dna_proteins_MHC[,-c(1,2)],1,sum)
pop.mhc<-c()
#function that checks if haplotype is present only in continent, or only in galapagos or in both
pop.tlr<-apply(frequency_dna_proteins_MHC[,-c(1,2)],1, function(x) if (x[2] == 0) append(pop.mhc,"gala") else if ((x[2] != 0)& ( 1 %in% x[-2])) append(pop.mhc,"both") else append(pop.mhc, "cont") )
}

property_files<-c("formal_charge.tsv","hydropathy.tsv","miyazawa-jernigan.tsv","prot_segregating.tsv","volume.tsv")

setwd(tlrdirectory)


property_df<-(paste(tlr,gsub(".tsv","",property_files),sep = "")) # generate data.frame names

for (i in 1:length(property_files)) { # read tables and assign them with names
  assign(property_df[i], read.csv(paste(tlrdirectory,property_files[i],sep ="")))
  }

property_df_traits<-property_df[-4] # use only proproperties, not amk

###heatmapfunction###

ht<-function (x) { # function that generates heatmap for each property
  d<-get(x)
  #print(x)
  m<-data.matrix(d)
  heatmap(m,Colv=NA, main = x)
}

###generata alignment view of segreg sites####
tlramk<-(get(property_df[4])[-1])
uniqAMK<-apply(tlramk,2,unique) # get unique amk
uniqAMK<-unique(unlist(uniqAMK)) # for TLR15
#uniqAMK<-unique(c(uniqAMK[1,],uniqAMK[2,]))
spect <- palette(rainbow(length(uniqAMK))) # vector of colors with the size of unique amk
tlramkf<-cbind(tlramk,hapfreq,pop.tlr)

amkcol<-tlramk # a new matrix with the same dimension as the alignment 
amkcol[]<-spect[match(apply(tlramk,2,as.character),uniqAMK)] # fill the matrix with color codes acording amk
amkcol<-cbind(amkcol,rep("black",length(hapfreq)),rep("black",length(hapfreq)))
tt <- ttheme_default(base_size = 16,core=list(fg_params = list(col = apply(amkcol,2,as.character)),
                                    bg_params = list(col=NA)),
                          rowhead=list(bg_params = list(col=NA)),
                          colhead=list(bg_params = list(col=NA))) # theme mapping for grid graphics

amkt <- tableGrob(tlramkf,theme = tt) # generate tablegrob

###add title to the grid graphics###
title <- textGrob(paste(tlr,"segregating aminoacid positions and heatmaps of their properties"),gp=gpar(fontsize=25))
padding <- unit(5,"mm")
amkt <- gtable_add_rows(amkt, heights = grobHeight(title) + padding,pos = 0)
amkt <- gtable_add_grob(amkt, title,1, 1, 1, ncol(amkt))

pdf(paste(tlr,"_alignment.pdf",sep=""), width=23, height = 50, pointsize = 20) # separated from properties in MHC due to size
grid.draw(amkt)
dev.off()
### print the grid with alignment and heatmaps ### 
pdf(paste(tlr,"_properties.pdf",sep=""), width=20, height = 20, pointsize = 20)
lapply(property_df_traits,ht)
dev.off()

###pca with all the properties###
all_property_matrix<- matrix(, nrow=dim(get(property_df[1]))[1], ncol=0,) # empty matrix for pca
for (i in property_df_traits) { # load all property matrices to a single one
  d<-get(i)
  #print(i)
  m<-data.matrix(d)
  colnames(m)<-paste(colnames(m),i,sep="")
  all_property_matrix<-cbind(all_property_matrix,m)
}

rownames(all_property_matrix)<-c(1:dim(all_property_matrix)[1]) # assign row names (haplotype name)

#all_property_matrix_gala<-all_property_matrix[galahap,] # define only galapagos matrix

tlr.pca<-prcomp(all_property_matrix) # compute pca, change input if you want only gala or all
#print(tlr.pca)
svg("variance_dist_pca.svg") # in MHC it would need three axes, the variablity is high
plot(tlr.pca,type="l")
dev.off()
sumpc<-summary(tlr.pca) # variance explained by each axis ####!!!!this needs to be fixed a bit so it works as a running script not stupind Rstudio clicking 
PC1_var<-round(sumpc$importance[2,1],3)
PC2_var<-round(sumpc$importance[2,2],3)
PC3_var<-round(sumpc$importance[2,3],3)

xpc1name<-paste("PC1_",PC1_var,sep="") # axis labels
xpc2name<-paste("PC2_",PC2_var,sep="")#shift+alt+arrow (duplicate line)
xpc3name<-paste("PC3_",PC3_var,sep="")#shift+alt+arrow (duplicate line)

scores<-as.data.frame(tlr.pca$x)
scores<-cbind(scores,pop.tlr)
svg(paste("PCA_PC12_physiochemical_all_",tlr,".svg",sep="")) #export to svg
ggplot(data=scores, aes(x = PC1, y = PC2, label = rownames(scores), colour = pop.tlr))+
  geom_point()+
  labs(x=xpc1name,y=xpc2name) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_text_repel(colour = "black", alpha = 0.9, size = 3 ) + # text is not overlapping
  ggtitle(paste(tlr," all haplotypes in physiochemical space",sep="")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()
