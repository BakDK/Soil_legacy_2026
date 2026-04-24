library(Rtsne)
library(vegan)
library(ampvis2)
library(phyloseq)
library(ggforce)
library(ggpubr)

##======= FIGURE 2 ======== ##
# tSNE of the different time points

setwd("..");ordination_data<-readRDS("Diversity_Analysis_Results/Non_Rarefied_Genus/full_analysis_results.rds");setwd("Figures")
tsne_data_sam<-ordination_data$tsne_data

tsne_data_sam$response_var<-as.character(tsne_data_sam$response_var)
tsne_data_sam$response_var[tsne_data_sam$response_var %in% "late"]<- "Late"
tsne_data_sam$response_var[tsne_data_sam$response_var %in% "middle"]<- "Middle"
tsne_data_sam$response_var[tsne_data_sam$response_var %in% "early"]<- "Early"

tsne_data_sam$Time<-gsub("DAS","", tsne_data_sam$Time)

tsne_data_sam$Time<-factor(tsne_data_sam$Time, 
                       levels = c( "11",
                                  "14","17","31",
                                  "158","187","215",
                                  "227","256","318"))


tsne_root_samples<-ggplot(tsne_data_sam,aes(x = Dim1, y = Dim2 )) +geom_point(aes(shape = Fertilizer , color = Time))+
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 18), 
        legend.title = element_text(size = 18),legend.text = element_text(size = 13),
        legend.position = "right", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background =  element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        strip.background = element_rect(fill="white", color = "white"))+
  scale_colour_brewer(palette = "Spectral")+
  ggforce::geom_mark_ellipse(alpha = 0.1,aes(fill = as.factor(response_var)), expand = (unit(0.5,"mm")))+
  guides(shape = guide_legend(title = "Fertility Level"),fill = guide_legend(title ="Growth period"), color = guide_legend(title = "Time (DAS)"))+
  ylab("t-SNE 2")+xlab("t-SNE 1")

#combined_plot<-ggarrange(fertilizer_soil_plot, tsne_root_samples, ncol =1, labels = c("a","b"))
x11(width = 7, height = 6)
tsne_root_samples
ggsave("Figures/figure_2.png",  dpi = 300)


# Supplementary figure S10
umap_data<-ordination_data$umap_data

umap_data$response_var<-as.character(umap_data$response_var)
umap_data$response_var[umap_data$response_var %in% "late"]<- "Late"
umap_data$response_var[umap_data$response_var %in% "middle"]<- "Middle"
umap_data$response_var[umap_data$response_var %in% "early"]<- "Early"

umap_data$Time<-gsub("DAS","", umap_data$Time)

umap_data$Time<-factor(umap_data$Time, 
                           levels = c( "11",
                                       "14","17","31",
                                       "158","187","215",
                                       "227","256","318"))


umap_root_samples<-ggplot(umap_data,aes(x = Dim1, y = Dim2 )) +geom_point(aes(shape = Fertilizer , color = Time))+
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 18), 
        legend.title = element_text(size = 18),legend.text = element_text(size = 13),
        legend.position = "right", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background =  element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        strip.background = element_rect(fill="white", color = "white"))+
  scale_colour_brewer(palette = "Spectral")+
  ggforce::geom_mark_ellipse(alpha = 0.05,aes(fill = as.factor(response_var)), expand = (unit(0.5,"mm")))+
  guides(shape = guide_legend(title = "Fertility Level"),fill = guide_legend(title ="Growth period"),color = guide_legend(title = "Time (DAS)"))+
  ylab("UMAP 2")+xlab("UMAP 1")

x11(width = 7, height = 6)
umap_root_samples
ggsave("Figures/Supp_fig_s10.png",  dpi = 300)
dev.off()
#Supplementary figure

## import data
phyl_gen<-readRDS("Input_files/Phyloseq_genus.rds")

## Soil 
Soil_only<-subset_samples(phyl_gen, Zone %in% "BulkSoil" )
abund_soil<-otu_table(Soil_only)
dm_soil<-vegdist(abund_soil, "bray")

# Run T-SNE
set.seed(5)
tsne_res_soil<-Rtsne(dm_soil, dims = 2, peplexity = 50, verbose = TRUE, max_iter = 5000)


#Plot 
tsne_df_soil <- data.frame(
  X = tsne_res_soil$Y[, 1],
  Y = tsne_res_soil$Y[, 2],
  Time = sample_data(Soil_only)$Time,
  Fert = sample_data(Soil_only)$Fertilizer,
  Ino = sample_data(Soil_only)$Inoculum.
)
tsne_df_soil$Time[tsne_df_soil$Time %in% "T0"]<-"22 DBS"
tsne_df_soil$Time[tsne_df_soil$Time %in% "T1"]<-"19 DBS"
tsne_df_soil$Time[tsne_df_soil$Time %in% "3DBS"]<-"3 DBS"

#
fertilizer_soil_plot<-ggplot(tsne_df_soil, aes(x = X, y = Y,color = factor(Fert), shape = factor(Time)))+
  geom_point(size = 3)+xlab("t-SNE 1")+ylab("t-SNE 2")+theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  guides(color=guide_legend(title="Fertility Level"),shape =guide_legend(title="Sampling Day"))+
  xlim(-15,15)+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 20))

x11(width =7, height = 6)
fertilizer_soil_plot
ggsave("Plots/tsne_soil_Fertility_supp_plot.png",dpi = 300)


