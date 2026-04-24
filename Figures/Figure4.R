#Figure 4

# This are based on both rarefied and non-rarefied data
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(ampvis2)

setwd("..");genus_all_4_list<-readRDS("Outputs/biomarker_list.rds");setwd("Figures")

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

# relative abundance of these biomarkers
amp2_all4$metadata$Time<-gsub("DAS","",amp2_all4$metadata$Time)
amp2_all4$metadata$Time <-factor(amp2_all4$metadata$Time, 
                                 levels = c("T0","T1","3DBS", "11",
                                            "14","17","31",
                                            "158","187","215",
                                            "227","256","318"))

genus_ranked_mean<-readRDS("Outputs/genus_ranked.rds")

amp2_all4$tax$Genus<-factor(amp2_all4$tax$Genus, levels = genus_ranked_mean$Genus)
amp2_all4$metadata$Fert[amp2_all4$metadata$Fertilizer %in% "M1P1"]<-"M"
amp2_all4$metadata$Fert[amp2_all4$metadata$Fertilizer %in% "N1K1"]<-"NK"
amp2_all4$metadata$Fert[amp2_all4$metadata$Fertilizer %in% "N1P2K2"]<-"NPK"

amp2_all4$metadata$Period[amp2_all4$metadata$Time %in% c("11","14","17","31")]<-"Early"
amp2_all4$metadata$Period[amp2_all4$metadata$Time %in% c("158","187","215","227")]<-"Middle"
amp2_all4$metadata$Period[amp2_all4$metadata$Time %in% c("256","318")]<-"Late"
amp2_all4$metadata$Period[amp2_all4$metadata$Time %in% c("3DBS")]<-"Bulk Soil"

# Plot
heatmap_v1<-amp2_all4 %>% amp_subset_samples(!Time %in% c("T0","T1")) %>%
  amp_heatmap(facet_by = c("Time"), group_by = "Fert", tax_show = 45, 
              tax_aggregate = "Genus",
              normalise = FALSE,
              plot_values = FALSE,
              order_y_by = rev(genus_ranked_mean$Genus),
              color_vector = c("White","Blue"))+ 
  theme(strip.background = element_rect(fill="white", color = "white"))+
  theme(axis.text.x = element_text(size = 9, angle = 90), axis.title = element_text(size = 22), axis.text.y = element_text(size = 10))+
  labs(colour = "Rel. abundance (%)")
  
x11(width = 9, height =7)
heatmap_v1
ggsave("Plots/Figure4.png", dpi = 300)
ggsave("Plots/Figure4.svg")
dev.off()