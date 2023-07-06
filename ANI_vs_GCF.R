### Cytoscape network analyses on BigSCAPE results

library(tidyverse)
library(dplyr)
library(stringr)
library(ggpubr)
library(RCy3)


# 
# #Import network file
# df_network_0.30 <- read.csv(file = "../for_cytoscape_antismash_6.1.1/2023-05-12 18_58_45_df_network_0.30.csv")
# df_GCF_0.30 <- read.csv(file = "../for_cytoscape_antismash_6.1.1/2023-05-12 18_58_45_df_clusters_0.30.csv")
# df_fam_name_0.30 <- read.csv(file = "../for_cytoscape_antismash_6.1.1/2023-05-12 18_58_45_df_families_0.30.csv")
# df_knowns_0.30 <- read.csv(file = "../for_cytoscape_antismash_6.1.1/2023-05-12 18_58_45_df_known_0.30.csv")
# 
# #df_abs_pres_0.30 <- read.csv(file = "../for_cytoscape_antismash_6.1.1/2023-05-12 18_58_45_df_family_presence_0.30.csv")
# df_network_0.30$bgc_id <- df_network_0.30$Clustername.1
# #df_network_0.30$bgc_id_clustername.2 <- df_network_0.30$Clustername.2
# 
# 
# df_network_0.30_meta <- df_GCF_0.30 %>% select(bgc_id, product, genome_id, gcf_0.30) %>% left_join(df_network_0.30, ., by="bgc_id")
# 
# #remove _contigs-fixed from genome_id 
# df_network_0.30_meta$genome_id <- gsub("\\_contigs-fixed", "", df_network_0.30_meta$genome_id)
# 
# 
# df_network_0.30_meta <- df_knowns_0.30 %>% 
#   select(bgc_id, biosyn_class, compounds, gcf_0.30) %>% left_join(df_network_0.30_meta, ., by="bgc_id")
# 
# #Add Cluster no to the table
# df_network_0.30_meta <- nMDS_1_2_meta %>% select(ClusterNumber, genome_id, Spot_S.aur., S3320) %>% left_join(df_network_0.30_meta, ., by="genome_id")
# 
# 
# #df_network_0.30_meta$genome_id 
# #df_network_0.30_meta %>% filter(Clustername.2 == "c_000000000013.region002")
# #df_network_0.30_meta$Clustername.2 %>% grepl("c_000000000013.region002",.)
# 
# str(df_knowns_0.30)
# 
# nodes=df_network_0.30[,2:3]
# 
# # edges <- data.frame(source = df_network_0.30$Clustername.1, target= df_network_0.30$Clustername.2, 
# #                     score = df_network_0.30$Jaccard.index)
# 
# # edges <- data.frame(source=df_network_0.30[,"Clustername.1"],target=df_network_0.30[,"Clustername.2"],
# #                      Jaccard.index=df_network_0.30[,"Jaccard.index"],
# #                     stringsAsFactors=FALSE)
# 
# edges <- data.frame(source=df_network_0.30[,"Clustername.1"],target=df_network_0.30[,"Clustername.2"],
#                     Squared.similarity=df_network_0.30[,"Squared.similarity"],
#                     stringsAsFactors=FALSE)
# 
# 
# createNetworkFromDataFrames(edges=edges, title="GCF network", collection = "anvio vibrio collection")
# 
# #data.key.column is the node 
# loadTableData(df_network_0.30_meta, data.key.column="Clustername.1")
# 
# 
# #Andrimid bgc fam cluster 2514
# df_network_0.30_meta %>% filter(gcf_0.30.y == "2514")
# 
# 
# getTableColumns('node')
# 


#Import network file
df_network_0.30 <- read.csv(file = "../for_cytoscape_antismash_6.1.1/2023-05-12 18_58_45_df_network_0.30.csv")
df_GCF_0.30 <- read.csv(file = "../for_cytoscape_antismash_6.1.1/2023-05-12 18_58_45_df_clusters_0.30.csv")
df_fam_name_0.30 <- read.csv(file = "../for_cytoscape_antismash_6.1.1/2023-05-12 18_58_45_df_families_0.30.csv")
df_knowns_0.30 <- read.csv(file = "../for_cytoscape_antismash_6.1.1/2023-05-12 18_58_45_df_known_0.30.csv")


edges <- data.frame(source=df_network_0.30[,"Clustername.1"],target=df_network_0.30[,"Clustername.2"],
                    Squared.similarity=df_network_0.30[,"Squared.similarity"],
                    stringsAsFactors=FALSE)


createNetworkFromDataFrames(edges=edges, title="GCF network", collection = "anvio vibrio collection")

nodes <- getTableColumns('node')
nodes_2 <- nodes

nodes_2$bgc_id <- nodes$name

nodes_2 <- df_GCF_0.30 %>% select(bgc_id, product, genome_id, gcf_0.30) %>% left_join(nodes_2, ., by="bgc_id")

#remove _contigs-fixed from genome_id 
nodes_2$genome_id <- gsub("\\_contigs-fixed", "", nodes_2$genome_id)


nodes_2 <- df_knowns_0.30 %>% 
  select(bgc_id, biosyn_class, compounds, gcf_0.30) %>% left_join(nodes_2, ., by="bgc_id")



#Add Cluster no to the table
nodes_2 <- df_Vibrio %>% select(genome_id, Spot_S.aur., S3320, ribS2_taxon, Latitude) %>% left_join(nodes_2, ., by="genome_id")
nodes_2$Latitude <- abs(nodes_2$Latitude)

nodes_2 %>% filter(is.na(Spot_S.aur.))

createNetworkFromDataFrames(nodes=nodes_2, edges=edges, title="GCF network", collection = "anvio vibrio collection")

setBackgroundColorDefault('white', style.name = "Sample2")



###

df_abs_pres_0.30 <- read.csv(file = "../for_cytoscape_antismash_6.1.1/2023-05-12 18_58_45_df_family_presence_0.30.csv", header = T)

unique(colnames(df_abs_pres_0.30))

###

df_knowns_0.30 %>% filter(compounds == "andrimid") #BGF 2514

nodes_2 %>% filter(gcf_0.30.x == "2514")

df_network_0.30 %>% filter(Clustername.2 == "c_000000000013.region002")

nodes_2 %>% filter(!is.na(genome_id)) %>% summarise(n_distinct(name))


##### ANI vs BGC similarity #####
df_network_0.30 <- read.csv(file = "../for_cytoscape_antismash_6.1.1/2023-05-12 18_58_45_df_network_0.30.csv")


#removing diagonal combinations
df_network_0.30_filt <- df_network_0.30[df_network_0.30$Clustername.1!=df_network_0.30$Clustername.2,]

#mutate genome_id info to the dataframe
df_GCF_0.30 <- read.csv(file = "../for_cytoscape_antismash_6.1.1/2023-05-12 18_58_45_df_clusters_0.30.csv")
df_GCF_0.30$genome_id <- gsub("\\_contigs-fixed", "", df_GCF_0.30$genome_id)

df_GCF_0.30 <- df_GCF_0.30 %>% mutate(Clustername.1 = bgc_id, Clustername.2 = bgc_id)
df_GCF_0.30 <- df_GCF_0.30 %>% select(genome_id, Clustername.1, Clustername.2)

#mutate genome_id to dataframe
df_network_0.30_filt <- df_network_0.30_filt %>% left_join(., df_GCF_0.30, by = "Clustername.1") 
#clean up
df_network_0.30_filt <- df_network_0.30_filt[,-ncol(df_network_0.30_filt)] 
colnames(df_network_0.30_filt)[3] <- "Clustername.2"
df_network_0.30_filt <- df_network_0.30_filt %>% left_join(., df_GCF_0.30, by = "Clustername.2") 
colnames(df_network_0.30_filt)[2] <- "Clustername.1"
df_network_0.30_filt <- df_network_0.30_filt[,-ncol(df_network_0.30_filt)] 

colnames(df_network_0.30_filt)[ncol(df_network_0.30_filt)] <- "quary"
colnames(df_network_0.30_filt)[ncol(df_network_0.30_filt)-1] <- "genome_id"



#Import ANI data
ANIb_dist <- read.table(file = "../ANIb_output/ANIb_percentage_identity.tab", header = TRUE, row.names = 1,  sep = "\t")
str(ANIb_dist)

#clean up table
colnames(ANIb_dist) <-  gsub("\\_contigs.fixed", "", colnames(ANIb_dist))
rownames(ANIb_dist) <-  gsub("\\_contigs.fixed", "", rownames(ANIb_dist))

#convert to long format
ANIb_dist_long <- ANIb_dist %>%
  mutate(genome_id = rownames(.)) %>% 
  pivot_longer(cols = -genome_id, names_to = "quary", values_to = "ANI_similarity")

#mutate ANI similarity to df_network_0.30_filt
df_network_0.30_filt_ani <- left_join(df_network_0.30_filt, ANIb_dist_long)

#add comparison info
df_network_0.30_filt_ani$genome_comparison <- paste0(df_network_0.30_filt_ani$genome_id," - ",df_network_0.30_filt_ani$quary)
df_network_0.30_filt_ani$GCF_comparison <- paste0(df_network_0.30_filt_ani$genome_id," - ",df_network_0.30_filt_ani$quary)


#Plot ANI similarity vs. BGC fam similarity 
df_network_0.30_filt_ani$genome_comparison

df_network_0.30_filt_ani %>% filter(Shared.group!="") %>% group_by(Shared.group) %>%
  filter(n() >= 5, !is.na(ANI_similarity) ) %>%
  ggplot() +
  #ggplot(aes(abs(Latitude), as.numeric(bgcs_count), col=as.numeric(cluster))) +
  #geom_point(aes(size=bgcs_count)) + scale_size(range = c(-5, 10)) 
  geom_point(aes(1-Jaccard.index, 1-ANI_similarity, col=Shared.group), alpha=0.3) + facet_wrap("Shared.group")#+ stat_smooth(aes(Jaccard.index, ANI_similarity, col=Shared.group), method = "lm")

df_network_0.30_filt_ani %>% filter(Shared.group!="") %>% group_by(Shared.group) %>%
  filter(n() >= 10, !is.na(ANI_similarity) ) %>%
ggplot(aes(1-Jaccard.index, 1-ANI_similarity, col=Shared.group)) +
  geom_point() +
  geom_smooth(se = FALSE, method = lm) + facet_wrap("Shared.group") +  stat_cor()

df_network_0.30_filt_ani %>% filter(Shared.group!="") %>% group_by(Shared.group) %>%
  filter(n() >= 10, !is.na(DSS.index) ) %>%
  ggplot(aes(1-DSS.index, 1-ANI_similarity, col=Shared.group)) +
  geom_point() +
  geom_smooth(se = FALSE, method = lm) + facet_wrap("Shared.group") +  stat_cor(method = "spearman")

#Vizualize group betalactams GCF 2687, 2676 (largest BGF)
nodes_2_2676 <- nodes_2 %>% filter(gcf_0.30.x %in% c("2676"))
nodes_2_2687 <- nodes_2 %>% filter(gcf_0.30.x %in% c("2676"))

library(ggrepel)
BGF_2676_2687_plot <- df_network_0.30_filt_ani %>% filter(Clustername.1 %in% nodes_2_2676$name, Clustername.2 %in% nodes_2_2676$name,
                                    Clustername.1 %in% nodes_2_2687$name, Clustername.2 %in% nodes_2_2687$name) %>% 
  #filter(n() >= 10, !is.na(ANI_similarity) ) %>%
  ggplot(aes(1-Squared.similarity, 1-ANI_similarity)) +
  geom_point( )+
  geom_smooth(se = FALSE, method = lm) + 
  #geom_text_repel() + labs(title = "geom_text_repel()") + 
  stat_cor(method = "spearman")

BGF_2676_2687_plot


#Vizualize group GCF 2558 (aryl polyene)
nodes_2_2558 <- nodes_2 %>% filter(gcf_0.30.x %in% c("2558"))

df_network_0.30_filt_ani %>% filter(Clustername.1 %in% nodes_2_2558$name, Clustername.2 %in% nodes_2_2558$name
                                    ) %>% 
  #filter(n() >= 10, !is.na(ANI_similarity) ) %>%
  ggplot(aes(1-Squared.similarity, 1-ANI_similarity, label = genome_comparison)) +
  geom_point() +
  geom_smooth(se = FALSE, method = lm) + stat_cor(method = "spearman")# + geom_text_repel() + labs(title = "geom_text_repel()")

#Vizualize group GCF 2495 RiPP-like
nodes_2_2495 <- nodes_2 %>% filter(gcf_0.30.x %in% c("2495"))

df_network_0.30_filt_ani %>% filter(Clustername.1 %in% nodes_2_2495$name, Clustername.2 %in% nodes_2_2495$name
) %>% na.omit() %>%
  #filter(n() >= 10, !is.na(ANI_similarity) ) %>%
  ggplot(aes(1-DSS.index, 1-ANI_similarity, label = genome_comparison)) +
  geom_point() +
  geom_smooth(se = FALSE, method = lm) + 
  stat_cor(method = "spearman")# geom_text_repel() + labs(title = "geom_text_repel()")

#Vizualize group GCF 886 Vibrioferrin siderophore
nodes_2_886 <- nodes_2 %>% filter(gcf_0.30.x %in% c("886"))

df_network_0.30_filt_ani %>% filter(Clustername.1 %in% nodes_2_886$name, Clustername.2 %in% nodes_2_886$name
) %>% na.omit() %>%
  #filter(n() >= 10, !is.na(ANI_similarity) ) %>%
  ggplot(aes(1-Squared.similarity, 1-ANI_similarity, label = genome_comparison)) +
  geom_point() +
  geom_smooth(se = FALSE, method = lm) +
  stat_cor(method = "spearman"


           
           
##### Heatmap presence absence###

df_abs_pres_0.30 <- read.csv(file = "../for_cytoscape_antismash_6.1.1/2023-05-12 18_58_45_df_family_presence_0.30.csv")

#convert to long format
df_abs_pres_0.30_long <- df_abs_pres_0.30 %>%
#  mutate(genome_id = rownames(.)) %>% 
  pivot_longer(cols = -genome_id, names_to = "fam_id_0.30", values_to = "Presence")

#clean up
df_abs_pres_0.30_long$genome_id<- gsub("\\_contigs-fixed", "", df_abs_pres_0.30_long$genome_id)
df_abs_pres_0.30_long$fam_id_0.30<- gsub("X", "", df_abs_pres_0.30_long$fam_id_0.30)


# add fam_name and gcf_0.30 to the dataframe
df_fam_name_0.30 <- read.csv(file = "../for_cytoscape_antismash_6.1.1/2023-05-12 18_58_45_df_families_0.30.csv")
df_fam_name_0.30$fam_id_0.30 <- as.character(df_fam_name_0.30$fam_id_0.30)
df_abs_pres_0.30_long<- df_fam_name_0.30 %>% select(fam_id_0.30, fam_name) %>% left_join(df_abs_pres_0.30_long, .)

str(df_abs_pres_0.30_long)
#and gcf_0.30 to the dataframe
df_GCF_0.30 <- read.csv(file = "../for_cytoscape_antismash_6.1.1/2023-05-12 18_58_45_df_clusters_0.30.csv")
str(df_GCF_0.30)
df_GCF_0.30$fam_id_0.30 <- as.character(df_GCF_0.30$fam_id_0.30 )
df_abs_pres_0.30_long<- df_GCF_0.30 %>% select(fam_id_0.30, bigscape_class, gcf_0.30) %>% left_join(df_abs_pres_0.30_long, ., by = "fam_id_0.30")



# Heatmap 
df_abs_pres_0.30_long %>%
ggplot(aes(genome_id, fam_name, fill= Presence)) + 
  geom_tile()



