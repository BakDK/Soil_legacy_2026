# Determining the core taxa and plotting figure 1

library(tidyverse)
library(phyloseq)
library(ggVennDiagram)
library(ggplot2)
library(ampvis2)

#Import file
setwd("..");phyl_gen<-readRDS("Input_files/Phyl_rare_genus.rds");setwd("Figures") # This data only contains root samples.

unique(sample_data(phyl_gen)$Fertilizer )
unique(sample_data(phyl_gen)$Zone)
M1P1_phyl<-subset_samples(phyl_gen, Fertilizer %in% "M1P1")

core_M1P1 <-filter_taxa(M1P1_phyl, function(x) sum(x>1) > (0.95*length(x)), TRUE)
M1P1_df<-data.frame(tax_table(core_M1P1)) # 16 genera

N1K1_phyl<-subset_samples(phyl_gen, Fertilizer %in% "N1K1")
core_N1K1<-filter_taxa(N1K1_phyl, function(x) sum(x>1) > (0.95*length(x)), TRUE)
N1K1_df<-data.frame(tax_table(core_N1K1)) #26

N1P2K2_phyl<-subset_samples(phyl_gen, Fertilizer %in% "N1P2K2")
core_N1P2K2<-filter_taxa(N1P2K2_phyl, function(x) sum(x>1) > (0.95*length(x)), TRUE)
N1P2K2_df<-data.frame(tax_table(core_N1P2K2)) #26

M1P1_df$List<-paste(M1P1_df$Family,M1P1_df$Genus, sep = "_")
N1K1_df$List<-paste(N1K1_df$Family,N1K1_df$Genus, sep = "_")
N1P2K2_df$List<-paste(N1P2K2_df$Family,N1P2K2_df$Genus, sep = "_")


genus_list<-list(M1P1 = M1P1_df$List,
                 N1K1 =N1K1_df$List,
                 N1P2K2 =N1P2K2_df$List)
venn_dia_core<-ggVennDiagram(genus_list, category.names = c("M1P1","N1K1","N1P2K2"), label = "count", label_alpha = 0)+
  scale_fill_distiller(palette = "PuBu")+ theme(legend.position = "none")

#x11(width = 5, height =5)
#venn_dia_core
#ggsave("Plots/venndiagram.png", dpi = 300)


# relative abundance in each sample

core_sample_rel_N1P2K2<-sample_sums((core_N1P2K2))/sample_sums(N1P2K2_phyl)*100
mean(core_sample_rel_N1P2K2) # 45%
sd(core_sample_rel_N1P2K2) # 12%

core_sample_rel_N1K1<-sample_sums((core_N1K1))/sample_sums(N1K1_phyl)*100
mean(core_sample_rel_N1K1) # 47%
sd(core_sample_rel_N1K1) # 13%

core_sample_rel_M1P1<-sample_sums((core_M1P1))/sample_sums(M1P1_phyl)*100
mean(core_sample_rel_M1P1) # 39%
sd(core_sample_rel_M1P1) # 10%


core_abundance_1<-as.data.frame(core_sample_rel_M1P1) 
kolonnenavne<-c("Rel_ab")
colnames(core_abundance_1)<-kolonnenavne
core_abu_2<-as.data.frame(core_sample_rel_N1K1)
colnames(core_abu_2)<-kolonnenavne
core_abu3<-as.data.frame(core_sample_rel_N1P2K2)
colnames(core_abu3)<-kolonnenavne

All_rel<-rbind(core_abundance_1,core_abu_2,core_abu3)
All_rel$SampleNames<-rownames(All_rel)



abun_meta<-merge(data.frame(sample_data(phyl_gen)),All_rel, by.x = "SampleID",by.y = "SampleNames")

head(abun_meta)

abun_meta$Time<-gsub("DAS","", abun_meta$Time)
abun_meta$Time<-factor(abun_meta$Time, 
                                   levels = c("T0","T1","3DBS", "11",
                                              "14","17","31",
                                              "158","187","215",
                                              "227","256","318"))

core_rel_over_time<-abun_meta %>% 
    ggplot(aes(x = Time, y = Rel_ab, color = Fertilizer))+
  geom_point(position =position_dodge(width = 0.7))+theme_bw()+
  labs(color = "Sample")+theme(legend.position = "bottom")+
scale_color_manual(values=c("#E69F00","#999999",  "#56B4E9"))+
  ylab("Rel. abundance (%)")+
  guides(color=guide_legend(title="Fertilizer"),shape =guide_legend(title="Fertility level"))+
  theme(legend.margin = margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,0,-10))
  

#Subset the genera from the full data
list_of_taxa<-unique(c(rownames(N1P2K2_df),rownames(M1P1_df),rownames(N1K1_df)))

tax_table(phyl_gen)[,8]<-rownames(tax_table(phyl_gen))

phyl_gen_rel<-transform_sample_counts(phyl_gen, function(x) x/sum(x)*100)
total_core<-subset_taxa(phyl_gen_rel, ASV %in% list_of_taxa )

sample_sums(total_core)

sample_data(total_core)<-sample_data(total_core)[,c(8,1,2,3,4,5,6,7)]


total_core_amp<-amp_load(total_core)

amp_boxplot(total_core_amp, normalise = FALSE, tax_show =31)


# Identify the members that belong to the different groups
grp_mem<-process_data(Venn(genus_list))

full<-grp_mem$regionLabel
full$item

Shared_by_all<-full$item[[7]]
shared_inorg<-full$item[[6]]
Shared_M1N1<-full$item[[4]]
uniqueMan<-full$item[[1]]
uniqueP2<-full$item[[3]]
uniqueN1K1<-full$item[[2]]

List_of_all<-c(Shared_by_all,shared_inorg,Shared_M1N1,uniqueMan,uniqueP2,uniqueN1K1)

List_of_all2<-ifelse(grepl("_NA", List_of_all), List_of_all, sub("^[^_]+_", "",List_of_all))

List_of_all2[grep("_NA",List_of_all2)]<-paste("Unclass.",List_of_all2[grep("_NA",List_of_all2)], sep = " ")
List_of_all2<-sub("_NA", "",List_of_all2)

#Modify missing genus names for alignment of the vector and matrix

total_core_amp$tax$Genus[total_core_amp$tax$Genus == ""]<-paste("Genus_NA_Rep_",total_core_amp$tax$OTU[total_core_amp$tax$Genus %in% ""],sep ="")

total_core_amp$tax$Genus[grep("Genus_",total_core_amp$tax$Genus)]<-paste("Unclass.",total_core_amp$tax$Family[grep("Genus_",total_core_amp$tax$Genus)], sep = " ")

total_core_amp$tax$Genus<-factor(total_core_amp$tax$Genus, levels = List_of_all2)

total_core_amp2<-total_core_amp
total_core_amp2$tax<-total_core_amp$tax[order(total_core_amp$tax$Genus),]
#total_core_amp2$tax$Genus<-total_core_amp$tax$Comb

amp_boxplot(total_core_amp2, normalise = FALSE, tax_show =31,tax_aggregate = "Genus")
amp_heatmap(total_core_amp2, normalise = FALSE, tax_show =31,tax_aggregate = "Genus")

first_ver_box<-amp_boxplot(total_core_amp2, normalise = FALSE, tax_show =31, order_y = rev(total_core_amp2$tax$Genus),
            tax_aggregate = "Genus")+
  ylab("Rel. abundance (%)")+geom_vline(xintercept=c(18.5,8.5,7.5,5.5,2.5))
  annotate("text",label = c("Shared", "N1K1,N1P2K2","M1P1,N1K1","Unique M1P1","Unique N1P2K2","Unique N1K1"),
                            x = c(20,9,8,6,3,1), yend= c(55,55,56,57,58,51), size = 8)


sec_ver_box<-  first_ver_box+ annotate("text",label = c("Shared", "N1K1,N1P2K2","M1P1,N1K1","Unique M1P1","Unique N1P2K2","Unique N1K1"),
                          x = c(20,10,8,6.5,4,1.5), y= c(72,72,72,72,72,72), size = 2.7, hjust = 1)


#Combine all three plots
library(ggpubr)
all_together_plot<-ggarrange(ggarrange(venn_dia_core,core_rel_over_time, ncol= 2,
                                       labels = c ("A","B")),sec_ver_box, heights = c(0.4,0.6), ncol = 1, nrow = 2, labels = c("A", "C"))

x11(width = 10, height =9)
all_together_plot
ggsave("Plots/Figure1.png",dpi = 300)
dev.off()



