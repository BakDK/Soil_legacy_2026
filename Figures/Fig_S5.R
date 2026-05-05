# Determining the Shannon diversity and the Chao1 richness
library(phyloseq)
library(tidyverse)

phyl_ob<-readRDS("Input_files/Filtered_phyloseq25.rds")
min(sample_sums(phyl_ob))

sample_data(phyl_ob)$Time<-factor(sample_data(phyl_ob)$Time, 
                                     levels = c("T0","T1","3DBS", "11DAS",
                                                "14DAS","17DAS","31DAS",
                                                "158DAS","187DAS","215DAS",
                                                "227DAS","256DAS","318DAS"))


#First rarefy 100 times and take the average
## Create a matrix with all the samples as rows and 100 empty columns
shan_ra<-matrix(NA,nsamples(phyl_ob),100)

# Then rarefy once and then calculate  Shannon diversity, repeat 100 times. Do not set rngseed as the values will the be the same
#Do not set seed, otherwise it will be the same every time
for(i in 1:100){
  {Full16_rare<-rarefy_even_depth(phyl_ob, sample.size = min(sample_sums(phyl_ob)))
  Shannon<-plot_richness(Full16_rare ,x = "Sample_ID", measures = "Shannon")
    shan_ra[,i]<-Shannon$data$value}
}

print(shan_ra)
saveRDS(shan_ra,"Outputs/shannon.rds")


#Calculate mean of the 100 rarefied Shannon indices
shan_ra_me<-apply(shan_ra,MARGIN =  1, mean)

shan_mean_met<-cbind(sample_data(phyl_ob),Shannon =shan_ra_me)
unique(shan_mean_met$Time)
#  Plot the data
#shan_mean_met$Plot_var<-paste(shan_mean_met$Variety,shan_mean_met$Treatment, sep = "_")

Shann_plot_all_points<-shan_mean_met %>% 
  #filter(!Treatment %in% "Mock" & !Compartment %in% "Mock") %>% 
  ggplot(aes(x = Time, y = Shannon, color = Fertilizer))+
  geom_point(position =position_dodge(width = 0.4))+theme_bw()+
  labs(color = "Sample")+theme(legend.position = "bottom")+
  scale_color_manual(values=c("#E69F00","#999999",  "#56B4E9"))+
  labs(colour = "Fertily Level")

Shann_plot_all_points
ggsave("Figures/Plots/FigS5.png",Shann_plot_all_points,width = 10, height = 6, dpi = 300)

#Statistical testing

shan_mean_rhiz <- shan_mean_met %>% filter(Zone %in% "Root")

shan_mean_rhiz$Block<-sapply(strsplit(shan_mean_rhiz$Sample,"\\."),'[',2)
shan_mean_rhiz$Block2<-(gsub("[A-Z]","",shan_mean_rhiz$Block))
shan_mean_rhiz$Days<-as.integer(gsub("[A-Z]","",shan_mean_rhiz$Time))

str(shan_mean_rhiz)
shan_mean_rhiz$Fertilizer<-as.factor(shan_mean_rhiz$Fertilizer)
shan_mean_rhiz$Inoculum.<-as.factor(shan_mean_rhiz$Inoculum.)


m1_null<-lme4::lmer(Shannon ~ Days  + Inoculum. + (1| Block), data = shan_mean_rhiz, REML = FALSE)
m1_test_shan<-lme4::lmer(Shannon ~Days +Fertilizer + Inoculum. + (1|Block), data = shan_mean_rhiz, REML = FALSE) 
summary(m1_null)
anova(m1_null,m1_test_shan)
# This means that there is a small affect (chi square = 9.0125, p = 0.044)

# For time
shan_null<-lme4::lmer(Shannon ~Fertilizer +Inoculum. + (1|Block), data = shan_mean_rhiz, REML = FALSE)

anova(shan_null,m1_test_shan)
# Also significant ( chi square = 0.498, p = 0.46)


## Richness

#First I will rarefy 100 times and take the average
# Create a matrix with all the samples as rows and 100 empty columns
rich_ra<-matrix(NA,nsamples(phyl_ob),100)

# Then rarefy once and then calculate  Shannon diversity, repeat 100 times. Do not set rngseed as the values will the be the same
#Do not set seed, otherwise it will be the same every time
for(i in 1:100){
  {Full16_rare<-rarefy_even_depth(phyl_ob, sample.size = min(sample_sums(phyl_ob)))
  Richness<-plot_richness(Full16_rare ,x = "Sample_ID", measures = "Chao1")
  rich_ra[,i]<-Richness$data$value}
}


#Calculate mean of the 100 rarefied Shannon indices
rich_ra_me<-apply(rich_ra,MARGIN =  1, mean)

rich_mean_met<-cbind(sample_data(phyl_ob),Richness =rich_ra_me)
unique(rich_mean_met$Time)
#  Plot the data

richness_plot<-rich_mean_met %>% 
    ggplot(aes(x = Time, y = Richness, color = Fertilizer))+
  geom_point(position =position_dodge(width = 0.4))+theme_bw()+
  labs(color = "Sample")+theme(legend.position = "bottom")

ggsave("Test_plots/richness_all_data.png", richness_plot, width = 10,height = 6, dpi = 300)


## permanova
library(vegan)

rhizo_phyl<-subset_samples(phyl_ob, Zone %in% "Root" )
ASV_table<-data.frame((otu_table(rhizo_phyl))) # samples as rows
min(sample_sums(rhizo_phyl))
rhizo_dist_av<-avgdist(ASV_table, meanfun = "mean", dmethod = "bray", sample = min(sample_sums(rhizo_phyl)))
saveRDS(rhizo_dist_av, "rhizo_dist_av.rds")
head(rhizo_dist_av)

sample_data(rhizo_phyl)$Block<-sapply(strsplit(sample_data(rhizo_phyl)$Sample,"\\."),'[',2)
h1 <- with(data.frame(sample_data(rhizo_phyl)), how(nperm = 999, blocks = Block)) # Accounting for block effect
set.seed(5)
model<-adonis2(rhizo_dist_av~Time*Fertilizer*Inoculum., data = data.frame(sample_data(rhizo_phyl)), by = "terms", permutations = h1)

saveRDS(model,"Permanova_res.rds")
model
#saved as alpha_beta_data