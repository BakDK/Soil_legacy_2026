# Raw reads were processed in all_dada_r431.R

#Initial data processing of the 16S rRNA data from the future cropping project.
# The phyloseq object contains seed samples as well, which will be excluded from this analysis.
library(phyloseq)
library(tidyverse)
fut_crop<-readRDS("Fut_crop_phyloseq_25.rds")
fut_crop
min(sample_sums(fut_crop))

hist(sample_sums(fut_crop))
#Import metadata
library(xlsx)
setwd("..");metadata<-read.xlsx("20230119_FutureCropping_RawSequences_Metadata.xlsx",sheetIndex = 1);setwd("Data_processing")
head(metadata)
rownames(metadata)<-lapply(lapply(strsplit(metadata$Fastq_filename, "_"), `[`, 1:2), function(x) paste(x,collapse = "_"))
head(metadata)
metadata<-metadata[,-7] #Remove column with no relevance

metadata_phy<-sample_data(metadata)
sample_names(metadata_phy)
#Merging the metadata with the phyloseq object
fut_phy<-merge_phyloseq(fut_crop,metadata_phy)
fut_phy
#Number of reads
sample_sums(fut_phy)

# start by removing all Seed samples
Fut_root_soil<-subset_samples(fut_phy,!Zone %in% "Seed" )
Fut_root_soil<-prune_taxa(taxa_sums(Fut_root_soil)>0, Fut_root_soil) #16,114 ASVs (removed, i.e. only present in seed samples)


#Some samples have 0 reads. They are removed:
Fut_root_soil<-prune_samples(sample_sums(Fut_root_soil)>0,Fut_root_soil)
Fut_root_soil # down to 688 samples (from 691)

min(sample_sums(Fut_root_soil)) #1
max(sample_sums(Fut_root_soil)) #118,484

# Chloroplast, mitochondria and non-bacterial ASVs are removed
Fut_root_soil<-subset_taxa(Fut_root_soil, Family != "Mitochondria") #Removes 18,024 ASV
Fut_root_soil<-subset_taxa(Fut_root_soil, Class != "Chloroplast") #Removes 0 ASV
Fut_root_soil<-subset_taxa(Fut_root_soil, Kingdom = "Bacteria") #Removes 0 ASV

library(ampvis2)
Fut_amp<-amp_load(Fut_root_soil)

Fut_amp %>% amp_subset_samples(Zone %in% "NegControl") %>%
  amp_heatmap(group_by = "SampleID",
              tax_aggregate = "OTU",
              tax_show = 15,
              tax_add = "Genus") 

Fut_amp %>% amp_subset_samples(Zone %in% "BulkSoil" & Time %in% "T0") %>%
  amp_heatmap(group_by = "SampleID",
              tax_aggregate = "OTU",
              tax_show = 15,
              tax_add = "Genus") 

#Check negative controls and compare with other samples within each sequencing run:
E6vW_sample_names<-Fut_amp$metadata$SampleID[grep("Ev6W",Fut_amp$metadata$Fastq_filename)]
E6vW_samples<-amp_subset_samples(Fut_amp,SampleID %in% E6vW_sample_names)
E6vW_samples %>% amp_heatmap(group_by = "SampleID",
                             tax_aggregate = "OTU",
                             tax_show = 20,
                             tax_add = "Genus") #29 samples
#The three most abundant ASVs in the negative controls were not found in the real samples. 
# The fourth most abundant ASV was found, but at low rel. abundance (less than 0.6%) and not even in all negative controls. 
# These negative samples are just removed

#Check negative controls and compare with other samples within each sequencing run:
Ex1F_sample_names<-Fut_amp$metadata$SampleID[grep("Ex1F",Fut_amp$metadata$Fastq_filename)]
Ex1F_samples<-amp_subset_samples(Fut_amp,SampleID %in% Ex1F_sample_names)
Ex1F_samples %>% amp_heatmap(group_by = "Time",
                             tax_aggregate = "OTU",
                             tax_show = 20,
                             tax_add = "Genus") #100 samples
#Here, ASV3 (~45%) and ASV19 (~10%) were the most abundant in the negative controls. They were 0 in the regular samples. 
# There were some ASVs found in both samples and controls, but not in all controls, so nothing seems problematic.

Ex1B_sample_names<-Fut_amp$metadata$SampleID[grep("Ex1B",Fut_amp$metadata$Fastq_filename)]
Ex1B_samples<-amp_subset_samples(Fut_amp,SampleID %in% Ex1B_sample_names)
Ex1B_samples %>% amp_heatmap(group_by = "SampleID",
                             tax_aggregate = "OTU",
                             tax_show = 20) #100 samples
#Here, ASV3 (~55%) and ASV19 (~10%) were the most abundant in the negative controls. They were 0 in the regular samples. 
# There were some ASVs found in both samples and controls, but not in all controls, so nothing seems problematic.
#Same as for the previous sequencing run.

Ex1C_sample_names<-Fut_amp$metadata$SampleID[grep("Ex1C",Fut_amp$metadata$Fastq_filename)]
Ex1C_samples<-amp_subset_samples(Fut_amp,SampleID %in% Ex1C_sample_names)
Ex1C_samples %>% amp_heatmap(group_by = "Time",
                             tax_aggregate = "OTU",
                             tax_show = 20)
#Same pattern, again ASV3 and 19 are in high abundance in the controls

Ex6_sample_names<-Fut_amp$metadata$SampleID[grep("Ex6",Fut_amp$metadata$Fastq_filename)]
Ex6_samples<-amp_subset_samples(Fut_amp,SampleID %in% Ex6_sample_names)
Ex6_samples %>% amp_heatmap(group_by = "SampleID",
                             tax_aggregate = "OTU",
                            tax_add = "Genus",
                             tax_show = 25) #96 samples
  #ASV 88 only in the negative controls. 
# But ASV46 (Afipia), ASV44 (Udaeobacter) ASV30 (Udaeobacter) and ASV20 (Bradyrhizobium) are in comparable concentrations across most samples
# The Afipia could be an issue (but see below). The other three are common soil bacteria, and the samples from this sequencing run are soil samples

Ex8_sample_names<-Fut_amp$metadata$SampleID[grep("Ex8",Fut_amp$metadata$Fastq_filename)]
Ex8_samples<-amp_subset_samples(Fut_amp,SampleID %in% Ex8_sample_names)
Ex8_samples %>% amp_heatmap(group_by = "SampleID",
                            tax_aggregate = "OTU",
                            tax_add = "Genus",
                            tax_show = 25) # 75 samples
#Here the ASV46 (Afipia) is found in comparable proportions across the samples, but not negative controls.
# Only Halomonas (ASV43) should be removed 

Ex9_sample_names<-Fut_amp$metadata$SampleID[grep("Ex9",Fut_amp$metadata$Fastq_filename)]
Ex9_samples<-amp_subset_samples(Fut_amp,SampleID %in% Ex9_sample_names)
Ex9_samples %>% amp_heatmap(group_by = "SampleID",
                            tax_aggregate = "OTU",
                            tax_add = "Genus",
                            tax_show = 25) # 101 samples
# Here I will remove ASV 428 as a contaminant, ASV78 and maybe ASV5 

ExF_sample_names<-Fut_amp$metadata$SampleID[grep("ExF",Fut_amp$metadata$Fastq_filename)]
ExF_samples<-amp_subset_samples(Fut_amp,SampleID %in% ExF_sample_names)
ExF_samples %>% amp_heatmap(group_by = "SampleID",
                            tax_aggregate = "OTU",
                            tax_add = "Genus",
                            tax_show = 25)
#Here, I will also remove ASV 78

Fut_amp %>%
  amp_subset_taxa(tax_vector = "ASV78", normalise = TRUE) %>%
  amp_heatmap(group_by = "Time", normalise = FALSE, round = 4)
# 31 DAS is from the ExF sequencing run
# 11, 17 and 14 DAS are from the Ex9 run where ASV 78 was also a problem. The ASV will be removed

Fut_amp %>%
  amp_subset_taxa(tax_vector = "ASV428", normalise = TRUE) %>%
  amp_heatmap(group_by = "Time", normalise = FALSE, round = 4, tax_add = "Genus")
# remove as well
Fut_amp %>%
  amp_subset_taxa(tax_vector = "ASV5", normalise = TRUE) %>%
  amp_heatmap(group_by = "Time", normalise = FALSE, round = 4, tax_add = "Genus")

Fut_amp %>%
  amp_subset_taxa(tax_vector = "Pantoea", normalise = TRUE) %>%
  amp_heatmap(group_by = "Time", normalise = FALSE, round = 4, tax_add = "Genus", tax_aggregate = "OTU")
# Pantoea is a normal seed-borne bacteria, so it will be kept as it is not certain it is a contamination

Fut_amp %>%
  amp_subset_taxa(tax_vector = "ASV43", normalise = TRUE) %>%
  amp_heatmap(group_by = "Time", normalise = FALSE, round = 4, tax_add = "Genus")
# This ASV has the highest relative abundance in the samples where I suspect contamination so it will be removed.


#Remove ASVs 43, 428 and 78 and all the negative controls.

Futu_minus_ctrl<-Fut_amp %>% amp_subset_samples(!Zone %in% "NegControl")
contaminants<-c("ASV43","ASV428","ASV78")
Futu_almfinal<-Futu_minus_ctrl %>% amp_subset_taxa(tax_vector = contaminants, remove = TRUE)
#Filter samples with less than 1,000 reads.
Futu_1000<-amp_filter_samples(Futu_almfinal, minreads = 1000) # removes four samples

Futu_1000 %>% amp_ordinate(type = "PCA", distmeasure = "bray", transform = "none", filter_species = 0,
                          sample_color_by = "Time" , sample_label_by = "SampleID")

#Here, a few samples are clearly outliers. Check the rarefaction curves per day

Futu_1000 %>% amp_subset_samples(Time %in% "T0") %>%
  amp_rarefaction_curve(facet_by = "Fertilizer", color_by = "SampleID")
Futu_1000 %>% amp_subset_samples(Time %in% "T1") %>%
  amp_rarefaction_curve(facet_by = "Fertilizer", color_by = "SampleID")
Futu_1000 %>% amp_subset_samples(Time %in% "3DBS") %>%
  amp_rarefaction_curve(facet_by = "Fertilizer", color_by = "SampleID")
Futu_1000 %>% amp_subset_samples(Time %in% "11DAS") %>%
  amp_rarefaction_curve(facet_by = "Fertilizer", color_by = "SampleID")
Futu_1000 %>% amp_subset_samples(Time %in% "14DAS") %>%
  amp_rarefaction_curve(facet_by = "Fertilizer", color_by = "SampleID")
Futu_1000 %>% amp_subset_samples(Time %in% "17DAS") %>%
  amp_rarefaction_curve(facet_by = "Fertilizer", color_by = "SampleID")
Futu_1000 %>% amp_subset_samples(Time %in% "31DAS") %>%
  amp_rarefaction_curve(facet_by = "Fertilizer", color_by = "SampleID")
Futu_1000 %>% amp_subset_samples(Time %in% "158DAS") %>%
  amp_rarefaction_curve(facet_by = "Fertilizer", color_by = "SampleID")
Futu_1000 %>% amp_subset_samples(Time %in% "187DAS") %>%
  amp_rarefaction_curve(facet_by = "Fertilizer", color_by = "SampleID")
Futu_1000 %>% amp_subset_samples(Time %in% "215DAS") %>%
  amp_rarefaction_curve(facet_by = "Fertilizer", color_by = "SampleID")
Futu_1000 %>% amp_subset_samples(Time %in% "227DAS") %>%
  amp_rarefaction_curve(facet_by = "Fertilizer", color_by = "SampleID")
Futu_1000 %>% amp_subset_samples(Time %in% "256DAS") %>%
  amp_rarefaction_curve(facet_by = "Fertilizer", color_by = "SampleID")
Futu_1000 %>% amp_subset_samples(Time %in% "318DAS") %>%
  amp_rarefaction_curve(facet_by = "Fertilizer", color_by = "SampleID")
# Still  some samples with few reads. I will filter at 1000 reads

#Subset a few samples based on the PCoA:
Futu_1000 %>% amp_subset_samples(SampleID %in% "S1620-215-B-Ex9_S116") %>% .$metadata # 14 DAS, N1P2K2
# ~14,000 reads, but very low amount of ASVs compared to the rest of the samples
Futu_1000 %>% amp_subset_samples(Time %in% "14DAS" & Fertilizer %in% "N1P2K2") %>%
  amp_heatmap(group_by="SampleID", tax_aggregate = "Genus") # 73% Bacillus in that sample. Sample is removed

#Subset a few samples based on the PCoA:
Futu_1000 %>% amp_subset_samples(SampleID %in% "S1620-155-A-Ex9_S11") %>% .$metadata # 11 DAS, N1K1
# ~7,000 reads, but very low amount of ASVs compared to the rest of the samples
Futu_1000 %>% amp_subset_samples(Time %in% "11DAS" & Fertilizer %in% "N1K1") %>%
  amp_heatmap(group_by="SampleID", tax_aggregate = "Genus") # 60% Bacillus in that sample. Sample is removed

Futu_1000 %>% amp_subset_samples(SampleID %in% "S1620-311-B-ExF_S168") %>% .$metadata # 31 DAS, M1P1
# ~42,600 reads, but  low amount of ASVs compared to the rest of the samples
Futu_1000 %>% amp_subset_samples(Time %in% "31DAS" & Fertilizer %in% "M1P1") %>%
  amp_heatmap(group_by="SampleID", tax_aggregate = "Genus") # 85% Serratia in that sample. Sample is removed

Futu_1000_red<- amp_subset_samples(Futu_1000, !SampleID %in% c("S1620-215-B-Ex9_S116","S1620-155-A-Ex9_S11","S1620-311-B-ExF_S168"), removeAbsentOTUs = TRUE)

Futu_1000_red %>% amp_ordinate(type = "PCOA", distmeasure = "bray", transform = "none", filter_species = 0.1,
                           sample_color_by = "Time" , sample_label_by = "Time")
Futu_1000_red %>% amp_ordinate(type = "PCOA", distmeasure = "bray", transform = "none", filter_species = 0.1,
                               sample_color_by = "Fertilizer" )

saveRDS(Futu_1000_red, "Input_files/Future_crop_ampvis_2025.rds")

##=====##
# filter the phyloseq object as well
Fut_root_raw<-subset_samples(Fut_root_soil, !Zone %in% "NegControl")
sample_data(Fut_root_raw)$SampleID<-sample_names(Fut_root_raw)
Fut_root_raw2<-subset_samples(Fut_root_raw,!SampleID %in% c("S1620-215-B-Ex9_S116","S1620-155-A-Ex9_S11","S1620-311-B-ExF_S168"))
taxo_table<-as.data.frame(tax_table(Fut_root_raw2))
taxo_table$ASV<-rownames(tax_table(Fut_root_raw2))
tax_ok<-tax_table(as.matrix(taxo_table))
tax_table(Fut_root_raw2)<-tax_ok

Fut_root_raw3<-subset_taxa(Fut_root_raw2, !ASV %in% c("ASV43","ASV428","ASV78"))
Fut_phyl_fin<-prune_samples(sample_sums(Fut_root_raw3)>= 1000,Fut_root_raw3)
Fut_phyl_fin<-prune_taxa(taxa_sums(Fut_phyl_fin)>0,Fut_phyl_fin)

saveRDS(Fut_phyl_fin, "../Input_files/Filtered_phyloseq25.rds")

#Save a objects where ASVs are agglomerated at genus.
Fut_crop_genus<-tax_glom(Fut_phyl_fin, taxrank = "Genus", NArm = FALSE) #Reduces the data from +40k ASVs to 985 genera

saveRDS(Fut_crop_genus,"../Input_files/Phyloseq_genus.rds")

##==== below may not needed



#Export otu table
otutabel<-otu_table(fut_1000)
rownames(otutabel)<-sample_data(fut_1000)$Sample
saveRDS(otutabel,"otu_table.rds")
head(t(otu_table(fut_1000)))

library(tidyverse)
sample_data(fut_1000) %>%
  data.frame()%>%
  count(Time,Zone)
# Looks okay in terms of number of samples.

# After a successful first round of LASSO regression (September 1st 2023, data presented by Thomas Martini)
# 02-09-2023
#As opposed to last time, I will keep more samples

sample_data(fut_1000)$Combined<-paste(sample_data(fut_1000)$Time,sample_data(fut_1000)$Zone,sep = "_")
unique(sample_data(fut_1000)$Combined)
fut_1000_red<-subset_samples(fut_1000,!Combined %in% c("T0.Seed_Seed","NegControl_NegControl",
                                                      "318DAS_Seed","158DAS_Seed","187DAS_Seed", "215DAS_Seed","227DAS_Seed",
                                                      "256DAS_Seed"))

unique(sample_data(fut_1000_red)$Combined)

fut_1000_red<-prune_taxa(taxa_sums(fut_1000_red)>0, fut_1000_red)
#38,436 ASVs. I will agglomerate to Genus level




# Agglomerate at Genus level:
gen_fut_1000_red<-tax_glom(fut_1000_red, taxrank = "Genus" ) #647 genera in 806 samples


# Get the ASV count table

gen_fut_1000_red_rel<-transform_sample_counts(gen_fut_1000_red, function(x) x / sum(x) *100)

abundance_large<-as.data.frame(t(otu_table(gen_fut_1000_red_rel)))
#Change sample names to something meaningful:
colnames(abundance_large)<-sample_data(gen_fut_1000_red_rel)$Sample
tax_table(gen_fut_1000_red_rel)
navne_og_antal<-merge(as.data.frame(tax_table(gen_fut_1000_red_rel)),abundance_large, by ="row.names")
head(navne_og_antal)


write.csv(navne_og_antal,"Large_data_setabundance.csv")
#Import yield data
udbytte<-read.xlsx("Winter_wheat_yield.xlsx", sheetIndex = 1)
udbytte_red<- udbytte %>% select(Yield_kg.ha,SampleID)

#edit the metadata
infodata_big<-as.matrix(sample_data(gen_fut_1000_red_rel)) %>%
  as.data.frame()%>%
  select(Sample,Time,Fertilizer,Inoculum.,Zone)

rownames(infodata_big)<-sample_data(gen_fut_1000_red_rel)$Sample

#Make a new column with SampleID so the data frame can be merged with the yield
infodata_big$SampleID<-gsub("T[0-9]+.","",infodata_big$Sample)
infodata_big$SampleID<-gsub("\\.B","",infodata_big$SampleID)
infodata_big$SampleID<-gsub("\\.A","",infodata_big$SampleID)

#
infodata_big %>%
  subset(Zone %in% "BulkSoil")

udbytte_og_meta<-merge(infodata_big,udbytte_red, by.x = "SampleID", all.x = TRUE)

gsub("\\..*","","Ab.X.y")
udbytte_og_meta$Block<-gsub("\\..*","",udbytte_og_meta$SampleID)

write.csv(udbytte_og_meta,"udbytte_almost_full_data.csv")

# For simplicity, I will try to reduce the data set to contain only one type of fertilizer regime 
# and then agglomerate to Genus level.

sample_data(fut_1000)
fut_1000_N1P2K2<-subset_samples(fut_1000, Fertilizer %in% "N1P2K2")
sample_data(fut_1000_N1P2K2)

#Keep only the root samples
fut_1000_N1P2K2_root<-subset_samples(fut_1000_N1P2K2, Zone %in% "Root")
fut_1000_N1P2K2_root<-prune_taxa(taxa_sums(fut_1000_N1P2K2_root)>0, fut_1000_N1P2K2_root)
# 11955 ASVs left, which might still be too much
# Agglomerate at Genus level:
gen_N1P2K2_root<-tax_glom(fut_1000_N1P2K2_root, taxrank = "Genus" )
# This leaves 409 genera which should be okay
head(tax_table(gen_N1P2K2_root))
Oversigt_gen_N1P2K2_root<-as.data.frame(tax_table(gen_N1P2K2_root))
Oversigt_gen_N1P2K2_root$Kingdom<-paste("k__",Oversigt_gen_N1P2K2_root$Kingdom, sep ="")
Oversigt_gen_N1P2K2_root$Phylum<-paste("p__",Oversigt_gen_N1P2K2_root$Phylum, sep ="")
Oversigt_gen_N1P2K2_root$Class<-paste("c__",Oversigt_gen_N1P2K2_root$Class, sep ="")
Oversigt_gen_N1P2K2_root$Order<-paste("o__",Oversigt_gen_N1P2K2_root$Order, sep ="")
Oversigt_gen_N1P2K2_root$Family<-paste("f__",Oversigt_gen_N1P2K2_root$Family, sep ="")
Oversigt_gen_N1P2K2_root$Species<-paste(paste("s__",Oversigt_gen_N1P2K2_root$Genus, sep =""),"unclassified",sep="_")
Oversigt_gen_N1P2K2_root$Genus<-paste("g__",Oversigt_gen_N1P2K2_root$Genus, sep ="")

# Get the ASV count table

gen_N1P2K2_root_rel<-transform_sample_counts(gen_N1P2K2_root, function(x) x / sum(x) *100)

abundance<-as.data.frame(t(otu_table(gen_N1P2K2_root_rel)))
#Change sample names to something meaningful:
colnames(abundance)<-sample_data(gen_N1P2K2_root_rel)$Sample
navn_and_abundance<-merge(Oversigt_gen_N1P2K2_root,abundance, by = 'row.names')
navn_abundance<-navn_and_abundance[,-1]
#Then the different taxonomic levels need to be merged
abundance_fil<-cbind(paste(paste(paste(paste(paste(paste(navn_abundance$Kingdom,navn_abundance$Phylum,sep ="|"),navn_abundance$Class,sep = "|"),navn_abundance$Order,sep ="|"),navn_abundance$Family,sep ="|"),
navn_abundance$Genus,sep = "|"),navn_abundance$Species,sep="|"),navn_abundance)
colnames(abundance_fil)[1]<-"Taxonomy"
abundance_fil<-abundance_fil[,-2:-8]
write.csv(abundance_fil,"abundance.csv")
#Import yield data
udbytte<-read.xlsx("Winter_wheat_yield.xlsx", sheetIndex = 1)
udbytte_red<- udbytte %>% select(Yield_kg.ha,SampleID)

#edit the metadata
infodata<-as.matrix(sample_data(gen_N1P2K2_root_rel)) %>%
  as.data.frame()%>%
  select(Sample,Time,Inoculum.)
rownames(infodata)<-sample_data(gen_N1P2K2_root_rel)$Sample
infodata$SampleID<-sample_data(gen_N1P2K2_root_rel)$Sample
infodata$SampleID<-gsub("\\.B","",infodata$SampleID)
infodata$SampleID<-gsub("T[0-9]+.","",infodata$SampleID)

udbytte_samlet<-merge(infodata,udbytte_red, by.x = "SampleID")
write.csv(udbytte_samlet,"udbytte.csv")

