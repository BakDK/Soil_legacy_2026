# Plot Figure 3

# This are based on both rarefied and non-rarefied data
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(ampvis2)

NR12<-read.csv("BiomarkersForGrowthStages/A_rf_top_15_final_biomarkers.csv" )
R12<-  read.csv("BiomarkersForGrowthStages/B_rf_top_15_final_biomarkers.csv" )
NR13<-read.csv("BiomarkersForGrowthStages/C_rf_top_15_final_biomarkers.csv" )
R13<-read.csv("BiomarkersForGrowthStages/D_rf_top_15_final_biomarkers.csv" )

NR12$Rare<-"No"
R12$Rare<-"Yes"
NR13$Rare<-"No"
R13$Rare<-"Yes"

N12_mer<-rbind(NR12,R12)
ggplot(N12_mer, aes(x= robust_potency, y = Feature, color = Rare, fill = Rare))+
  geom_bar(stat="identity",position="dodge")+facet_grid(.~Class,scales = "free_y")
N12_mer_genus_both<-N12_mer %>%
  group_by(Feature,Class) %>%
  add_count() %>%
  filter(n==2) %>% ungroup()


ggplot(N12_mer_genus_both, aes(x= robust_potency, y = Feature, color = Rare, fill = Rare))+
  geom_bar(stat="identity",position="dodge")+facet_grid(.~Class,scales = "free_y")

# repeat for the 1,3 train set

N13_mer<-rbind(NR13,R13)
ggplot(N13_mer, aes(x= robust_potency, y = Feature, color = Rare, fill = Rare))+
  geom_bar(stat="identity",position="dodge")+facet_grid(.~Class,scales = "free_y")
N13_mer_genus_both<-N13_mer %>%
  group_by(Feature,Class) %>%
  add_count() %>%
  filter(n==2) %>% ungroup()


ggplot(N13_mer_genus_both, aes(x= robust_potency, y = Feature, color = Rare, fill = Rare))+
  geom_bar(stat="identity",position="dodge")+facet_grid(.~Class,scales = "free_y")

#Quite some overlap

## Merge the two data sets
N13_mer_genus_both$train<-"t13"
N12_mer_genus_both$train<-"t12"


Genus_both_mer<-rbind(N13_mer_genus_both,N12_mer_genus_both)

genus_all_rare<-Genus_both_mer %>%
  subset(Rare %in% "Yes") %>%
  group_by(Feature,Class) %>%
  add_count(name = "All") %>%
  filter(All >=2 ) %>% ungroup()

genus_all_rare$Class<-factor(genus_all_rare$Class, levels = c("early","middle","late"))

class.labs<-c("Early","Middle","Late")
names(class.labs)<-c("early","middle","late")

# de er i forvejen valgt på baggrund af at de er fundet i begge tests
Shared_ind_gen_rare<-genus_all_rare %>% 
  ggplot( aes(x= robust_potency, y = Feature, color = train, fill = train))+
  geom_bar(stat="identity",position="dodge")+facet_grid(.~Class,labeller=labeller(Class =class.labs))

# prøv med alle 4
genus_all_4<-Genus_both_mer %>%
  select(Feature, Class, robust_potency, Rare, train) %>%
  #subset(Rare %in% "Yes") %>%
  group_by(Feature,Class) %>%
  add_count(name = "All") %>%
  filter(All >=4 ) %>% ungroup()

genus_all_4$Class<-factor(genus_all_4$Class, levels = c("early","middle","late"))


#Order the plots

genus_all_4_mean<-genus_all_4 %>% 
  group_by(Feature,Class) %>% 
  summarise(Mean = mean(robust_potency),
            max = max(robust_potency),
            min = min(robust_potency)) %>% ungroup()

genus_all_4_mean # max or min are not very different for the 27 genera


# make a subset of the ampvis2 object
genus_all_4_mean$Feature<-sub("\\."," ",genus_all_4_mean$Feature)
genus_all_4_list<-unique(genus_all_4_mean$Feature)

setwd("..");saveRDS(genus_all_4_list,"Outputs/biomarker_list.rds");setwd("Figures")

#Import the phyloseq object
phyl_ob<-readRDS("Input_files/Phyloseq_genus.rds")
colnames(sample_data(phyl_ob))[1]<-"Sample_ID"
sample_data(phyl_ob)<-sample_data(phyl_ob)[,-8]
tax_table(phyl_ob)[,8]<-rownames(tax_table(phyl_ob))

amp_obj2<-amp_load(phyl_ob)

#Modify missing genus names for alignment of the vector and matrix
amp_obj2$tax$Genus[amp_obj2$tax$Genus %in% ""]<- paste("Genus_NA_Rep_",amp_obj2$tax$OTU[amp_obj2$tax$Genus %in% ""],sep ="")
# Check if the names are identical
setdiff(genus_all_4_list,amp_obj2$tax$Genus) 
#subset the data only using the genera identified with treeSHAP

amp2_all4<-amp_obj2 %>% amp_subset_taxa(tax_vector = genus_all_4_list , normalise = TRUE) # 27 genera

#Insert a new column in the tax table
amp2_all4$tax$Feature<-amp2_all4$tax$Genus

amp2_all4$tax$Genus[grep("Genus_",amp2_all4$tax$Genus)]<-paste("Unclass.",amp2_all4$tax$Family[grep("Genus_",amp2_all4$tax$Genus)], sep = " ")

Genus_mod_all<-merge(genus_all_4_mean,amp2_all4$tax, by = "Feature")

# order based on their values first in early, then middle and finally late.
genus_ranked_mean<-Genus_mod_all %>%
  select(-max,-min)%>%
  pivot_wider(id_cols = "Genus",names_from = "Class.x", values_from ="Mean") %>%
  arrange(desc(early),desc(middle),desc(late))

genus_ranked_mean$Genus<-factor(genus_ranked_mean$Genus, levels = genus_ranked_mean$Genus)


genus_ranked_long<-genus_ranked_mean %>% pivot_longer(!Genus,names_to = "Class.x", values_to = "Mean")

genus_ranked_long$Class.x<-factor(genus_ranked_long$Class.x, levels = c("early","middle","late"))


class.labs<-c("Early","Middle","Late")
names(class.labs)<-c("early","middle","late")

#Plot potency
robust_potency_plot_v1<-genus_ranked_long %>% 
  ggplot( aes(x= Mean, y = Genus, fill = Class.x))+
  geom_bar(stat="identity",position="dodge")+scale_y_discrete(limits = rev)+
  facet_grid(.~Class.x,labeller=labeller(Class.x =class.labs))
robust_potency_plot_v1


# improve aesthetics
plot_1st_ed<-robust_potency_plot_v1 + theme(axis.text = element_text(size = 12), axis.title = element_text(size = 22), 
                                            legend.title = element_text(size = 20),legend.text = element_text(size = 15),
                                            strip.text = element_text(size = 20),
                                            legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                            panel.background =  element_blank(),
                                            panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
                                            strip.background = element_rect(fill="white", color = "white"))+
  scale_fill_manual(values=c("black", "lightgrey", "darkgrey"), labels = c("Early","Middle","Late"))+
  scale_color_manual(values=c("black", "lightgrey", "darkgrey"), labels = c("Early","Middle","Late"))+
  labs(fill = "Growth Period",color = "Growth Period")+ ylab("")+xlab("Robust Potency")+
  theme(panel.spacing = unit(2, "lines"))+
  theme(plot.margin = unit(c(0, 0.25, 0.25, 0.25),
                           "inches"))+
  scale_x_continuous(limits=c(0,0.04),breaks = c(0,0.015, 0.03),expand=c(0,0))

plot_1st_ed

ggsave("Plots/Figure3.png",units = "in",width = 8, height = 7.8, dpi = 300)
plot_1st_ed
dev.off()
