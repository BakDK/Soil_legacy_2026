# Plot genus data from RF modelling 
# This are based on both rarefied and non-rarefied data
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(ampvis2)
library(here) # Crucial for robust file paths

# 1. Define Paths relative to the Git Root
# 'here()' automatically finds the 'Future_cropping' folder
input_dir  <- here("outputs", "legacy")

output_dir <- here("Figures")

file_path_N13<-file.path(input_dir,"legacy_non_rarefied_block_train13_test24/03_biomarker_rankings/final_biomarker_ranking_main.csv")
N13<-read.csv(file_path_N13)
file_path_R13<-file.path(input_dir,"legacy_rarefied_block_train13_test24/03_biomarker_rankings/final_biomarker_ranking_main.csv")
R13<-read.csv(file_path_R13)

file_path_N12<-file.path(input_dir,"legacy_non_rarefied_block_train12_test34/03_biomarker_rankings/final_biomarker_ranking_main.csv")
N12<-read.csv(file_path_N12)
file_path_R12<-file.path(input_dir,"legacy_rarefied_block_train12_test34/03_biomarker_rankings/final_biomarker_ranking_main.csv")
R12<-read.csv(file_path_R12)

# add a variable
N12$Rare<-"No"
R12$Rare<-"Yes"
N13$Rare<-"No"
R13$Rare<-"Yes"

N12_mer<-rbind(N12,R12)
ggplot(N12_mer, aes(x= robust_potency, y = Feature, color = Rare, fill = Rare))+
  geom_bar(stat="identity",position="dodge")+facet_grid(.~Class,scales = "free_y")
N12_mer_genus_both<-N12_mer %>%
  group_by(Feature,Class) %>%
  add_count() %>%
  filter(n==2) %>% ungroup()


ggplot(N12_mer_genus_both, aes(x= robust_potency, y = Feature, color = Rare, fill = Rare))+
  geom_bar(stat="identity",position="dodge")+facet_grid(.~Class,scales = "free_y")

# repeat for the 1,3 train set

N13_mer<-rbind(N13,R13)
ggplot(N13_mer, aes(x= robust_potency, y = Feature, color = Rare, fill = Rare))+
  geom_bar(stat="identity",position="dodge")+facet_grid(.~Class,scales = "free_y")
N13_mer_genus_both<-N13_mer %>%
  group_by(Feature,Class) %>%
  add_count() %>%
  filter(n==2) %>% ungroup()


ggplot(N13_mer_genus_both, aes(x= robust_potency, y = Feature, color = Rare, fill = Rare))+
  geom_bar(stat="identity",position="dodge")+facet_grid(.~Class,scales = "free_y")

## Merge the two data sets
N13_mer_genus_both$train<-"t13"
N12_mer_genus_both$train<-"t12"

Genus_both_mer<-rbind(N13_mer_genus_both,N12_mer_genus_both)

# identify the taxa that are found in the rare data for both training sets. 
#Since filtering of taxa that do not appear in non-rarefied AND rarefied have been removed, filtering is only done on the "rare" variable. 
genus_all_rare<-Genus_both_mer %>%
  subset(Rare %in% "Yes") %>%
  group_by(Feature,Class) %>%
  add_count(name = "All") %>%
  filter(All >=2 ) %>% ungroup()

genus_all_rare$Class<-factor(genus_all_rare$Class, levels = c("M1P1","N1K1","N1P2K2"))

class.labs<-c("M1P1","N1K1","N1P2K2")
names(class.labs)<-c("M1P1","N1K1","N1P2K2")

# 
Shared_ind_gen_rare<-genus_all_rare %>% 
  ggplot( aes(x= robust_potency, y = Feature, color = train, fill = train))+
  geom_bar(stat="identity",position="dodge")+facet_grid(.~Class,labeller=labeller(Class =class.labs))

# Find those taxa that appear in all 4 data sets.
genus_all_4<-Genus_both_mer %>%
  select(Feature, Class, robust_potency, Rare, train) %>%
  group_by(Feature,Class) %>%
  add_count(name = "All") %>%
  filter(All >=4 ) %>% ungroup()

genus_all_4$Class<-factor(genus_all_4$Class, levels = c("M1P1","N1K1","N1P2K2"))


#Calculate mean robust potency
genus_all_4_mean<-genus_all_4 %>% 
  group_by(Feature,Class) %>% 
  summarise(Mean = mean(robust_potency),
            max = max(robust_potency),
            min = min(robust_potency)) %>% ungroup()

genus_all_4_mean # 
# order based on their values first in M1P1, then N1K1  and finally N1P2K2.
genus_ranked_mean<-genus_all_4_mean %>%
  select(-max,-min)%>%
  pivot_wider(names_from = "Class", values_from ="Mean") %>%
  arrange(desc(M1P1),desc(N1K1),desc(N1P2K2))

genus_ranked_mean$Feature<-factor(genus_ranked_mean$Feature, levels = genus_ranked_mean$Feature)

genus_ranked_long<-genus_ranked_mean %>% pivot_longer(!Feature,names_to = "Class", values_to = "Mean")
genus_ranked_long$Class<-factor(genus_ranked_long$Class, levels = c("M1P1","N1K1","N1P2K2"))

# Plot draft
robust_potency_plot_v1<-genus_ranked_long %>% 
  ggplot( aes(x= Mean, y = Feature, fill = Class))+
  geom_bar(stat="identity",position="dodge")+scale_y_discrete(limits = rev)+
  facet_grid(.~Class,labeller=labeller(Class =class.labs))


# Modify slightly 

genus_ranked_long_red<-genus_ranked_long %>% filter(Mean > 0.0025)

robust_potency_plot_v2<-genus_ranked_long_red %>% 
  ggplot( aes(x= Mean, y = Feature, fill = Class))+
  geom_bar(stat="identity",position="dodge")+scale_y_discrete(limits = rev)+
  facet_grid(.~Class,labeller=labeller(Class =class.labs))

# improve aesthetics
plot_1st_ed<-robust_potency_plot_v2 + theme(axis.text = element_text(size = 10), axis.title = element_text(size = 22), 
                                            legend.title = element_text(size = 20),legend.text = element_text(size = 15),
                                            strip.text = element_text(size = 20),
                                            legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                            panel.background =  element_blank(),
                                            panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
                                            strip.background = element_rect(fill="white", color = "white"))+
  scale_fill_manual(values=c("#E69F00","Grey",  "#56B4E9") )+
  scale_color_manual(values=c("#E69F00","Grey",  "#56B4E9"))+
   ylab("")+xlab("Robust Potency")+
  theme(panel.spacing = unit(2, "lines"))+
  theme(plot.margin = unit(c(0, 0.25, 0.25, 0.25),
                           "inches"))+
  scale_x_continuous(limits=c(0,0.04),breaks = c(0,0.015, 0.03),expand=c(0,0))

plot_1st_ed

x11(width = 10, height = 8)
plot_1st_ed
ggsave("Figures/Plots/Suppl_fig_S12.png", dpi = 300)
dev.off()


# make a subset of the ampvis2 object
genus_all_4_list<-unique(genus_ranked_long_red$Feature) # using only those above 0.0025
genus_ranked_mean$Feature<-sub("\\."," ",genus_ranked_mean$Feature)

#Import the phyloseq object
setwd("..");phyl_ob<-readRDS("Input_files/Phyloseq_genus.rds")
colnames(sample_data(phyl_ob))[1]<-"Sample_ID"
sample_data(phyl_ob)<-sample_data(phyl_ob)[,-8]
tax_table(phyl_ob)[,8]<-rownames(tax_table(phyl_ob))

amp_obj2<-amp_load(phyl_ob)
amp_obj2$metadata$Time<-factor(amp_obj2$metadata$Time, 
                               levels = c("T0","T1","3DBS", "11DAS",
                                          "14DAS","17DAS","31DAS",
                                          "158DAS","187DAS","215DAS",
                                          "227DAS","256DAS","318DAS"))


#Modify missing genus names for alignment of the vector and matrix
amp_obj2$tax$Genus[amp_obj2$tax$Genus %in% ""]<- paste("Genus_NA_Rep_",amp_obj2$tax$OTU[amp_obj2$tax$Genus %in% ""],sep ="")
# Check if the names are identical
setdiff(genus_all_4_list,amp_obj2$tax$Genus) 
genus_all_4_list<-gsub("Clostridium.sensu.stricto.1","Clostridium sensu stricto 1", genus_all_4_list)
genus_all_4_list<-gsub("Cylindrospermum.PCC.7417","Cylindrospermum PCC-7417", genus_all_4_list)
genus_all_4_list<-gsub("SH.PL14","SH-PL14", genus_all_4_list)
genus_all_list_2<-gsub("CL500.29.marine.group" ,"CL500-29 marine group", genus_all_4_list)

amp_obj2$tax$Genus[grep("CL500",amp_obj2$tax$Genus)]

setdiff(genus_all_list_2,amp_obj2$tax$Genus) 


#subset the data only using the genera identified with treeSHAP

amp2_all4<-amp_obj2 %>% amp_subset_taxa(tax_vector = genus_all_list_2 , normalise = TRUE) # 36 genera

amp2_all4$tax$Genus<-factor(amp2_all4$tax$Genus, levels = genus_all_list_2)
setdiff(genus_all_list_2,amp2_all4$tax$Genus) # 0 differences


str(genus_all_list_2)
# Plot
heatmap_v1<-amp2_all4 %>% amp_subset_samples(!Time %in% c("T0","T1")) %>%
  amp_heatmap(facet_by = "Fertilizer", group_by = "Time", tax_show = 45, 
              tax_aggregate = "Genus",
              normalise = FALSE,
              plot_values = FALSE,
              #tax_add = "Family",
              order_y_by = rev(genus_all_list_2),
              color_vector = c("White","Blue"))+ 
  theme(strip.background = element_rect(fill="white", color = "white"))+
  theme(axis.text.x = element_text(size = 9, angle = 90), axis.title = element_text(size = 22), axis.text.y = element_text(size = 10))+
  labs(fill = "Rel. abundance (%) ")


ggsave("Figures/Plots/Figure6.svg",heatmap_v1, units = "in" ,width = 9, height = 7)
ggsave("Figures/Plots/Figure6.png",heatmap_v1, units = "in" ,width = 9, height = 7, dpi = 300)



