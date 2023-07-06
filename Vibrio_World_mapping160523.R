library(tidyverse)
library(dplyr)
library(stringr)
library(ggpubr)
library(ggrepel)


#Import Phylogenomic distances
#Import autoMLST distance matrix and use S3320
anvio_dist <- read.table(file = "anvio_tree160523.csv", header = FALSE, sep = ";")
#anvio_dist <- read.table(file = "anvio_tree230523_noRef_noPhotobacterium.csv", header = FALSE, sep = ";")

rownames(anvio_dist) <- anvio_dist$V1
anvio_dist <- anvio_dist[,-1]
anvio_dist_t <- as.data.frame(t(anvio_dist))
rownames(anvio_dist_t) <- rownames(anvio_dist) 
anvio_dist_S3320 <- as.data.frame(anvio_dist_t$S3320)
rownames(anvio_dist_S3320) <- rownames(anvio_dist_t)


anvio_dist_S3320 <- anvio_dist_S3320 %>% mutate(samples = rownames(anvio_dist_S3320))
colnames(anvio_dist_S3320)[1] <- "S3320"
colnames(anvio_dist_S3320)[2] <- "genome_id"

#Import genome assembly stats
df_seqfu_stats <- read.table(file = "../bgcflow/df_seqfu_stats.csv", header = TRUE, sep = ",")


#Import Metadata
# Clean up metadata table
Meta_data_Galathea <- read.table(file = "../Metadata/Meta_data_latitude.txt", header = TRUE, sep = "\t")

#Filter Vibrio genome data from the data table
Meta_data_Galathea_Vibrios <- Meta_data_Galathea %>% filter(Strain %in% anvio_dist_S3320$genome_id)

#Meta_data_Galathea_Vibrios$X16SrRNA_gene <- gsub("\\s+$", "", Meta_data_Galathea_Vibrios$X16SrRNA_gene)

#Convert antibiotic activity to 1 to 3 scale
Meta_data_Galathea_Vibrios <- Meta_data_Galathea_Vibrios[,-c(1, ncol(Meta_data_Galathea_Vibrios))]
colnames(Meta_data_Galathea_Vibrios)[1] <- "genome_id"

Meta_data_Galathea_Vibrios[] <- lapply(Meta_data_Galathea_Vibrios, function(x) gsub("[()]", "", x))


replace_plus_values <- function(x) {
  ifelse(x == "-", 0,
         ifelse(x == "+", 1, 
                ifelse(x == "++", 2, 
                       ifelse(x == "++++", 4, 
                              ifelse(x == "+++", 3, x)))))
}

# Use mutate_all() to apply the function to all columns
Meta_data_Galathea_Vibrios <- Meta_data_Galathea_Vibrios %>% mutate_all(replace_plus_values)


# Add taxonomic informaton
df_tax <- read.table(file = "../bgcflow_antismash_gbk_regions/collection-ref_filtered.txt", header = TRUE, sep = "\t")
colnames(df_tax)[1] <- "genome_id"

df_tax_ribS2  <- read.table(file = "../external-genomes_refs_sgc_species-LONG-FORMAT.txt", header = TRUE, sep = "\t")
colnames(df_tax_ribS2 )[2] <- "genome_id"

df_tax_ribS2 <- filter(df_tax_ribS2 , times_observed ==1) %>% select("genome_id", "taxon")
colnames(df_tax_ribS2 )[2] <- "ribS2_taxon"

df_tax <- left_join(df_tax_ribS2, df_tax, by="genome_id")

# antiSMASH data frame ########

df_antismash <- read.table(file = "../bgcflow/tables/df_antismash_6.1.1_summary.csv", header = TRUE, sep = ",")
df_antismash$bgcs_count


colnames(df_antismash)[1] <- "genome_id"


#Join all data frames
df_Vibrio <- left_join(Meta_data_Galathea_Vibrios, df_seqfu_stats)
df_Vibrio <- left_join(df_Vibrio, anvio_dist_S3320)
df_Vibrio <- left_join(df_Vibrio, df_antismash)
df_Vibrio <- left_join(df_Vibrio, df_tax)


colnames(df_Vibrio)[22] <- "S3320"

#Convert cat to num
str(df_Vibrio)
df_Vibrio$Spot_V.ang <- as.numeric(df_Vibrio$Spot_V.ang)
df_Vibrio$Spot_S.aur.  <- as.numeric(df_Vibrio$Spot_S.aur. )
df_Vibrio$mmZone_V.ang <- as.numeric(df_Vibrio$mmZone_V.ang)
df_Vibrio$mmZone_S.aur. <- as.numeric(df_Vibrio$mmZone_S.aur.)
df_Vibrio$Longitude <-  as.numeric(df_Vibrio$Longitude)
df_Vibrio$Latitude <-  as.numeric(df_Vibrio$Latitude)

df_Vibrio$Longitude <-  as.numeric(df_Vibrio$Longitude)

#Add labels for plotting
df_Vibrio$labels <- ifelse(df_Vibrio$genome_id=="S3320","S3320", NA)
str(df_Vibrio)

#Add groups described by phylogenetic distance:
df_Vibrio %>% filter(Taxonomy == "aquimaris")
anvio_dist_vs_latitude_plot <-  
  df_Vibrio %>% 
  ggplot(aes(S3320, Latitude, col=S3320)) +
  #ggplot(aes(abs(Latitude), as.numeric(bgcs_count), col=as.numeric(cluster))) +
  #geom_point(aes(size = bgcs_count)) + ylim(-80,80) +
  geom_point(aes(size = bgcs_count)) + ylim(-80,80) +
  #geom_point(size=2) + ylim(-80,80) +
  #geom_smooth() +
  #geom_smooth(method = lm, se=FALSE) +
  #stat_cor(method = "spearman", p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name="rho", size = 2)+
  #scale_size(range = c(0.1, 10)) +
  labs(x = "\nCore genes phylogenetic distance", y = "\nLatitude ") +
  theme_bw(base_size = 10) +
  scale_colour_gradient(low = "blue", high = "red", na.value = NA) +
  # scale_colour_steps2(low = "white",
  #                     mid = "blue",
  #                     high = "red",
  #                     midpoint = 3) + 
  #scale_size(range = c(-5, 10)) +
  #scale_fill_manual(values = c("#f25a13", "#580a82"))+
  #scale_color_manual(values = c("#f25a13", "#580a82"))+
  geom_label_repel(aes(label = labels),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme(#panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    #axis.title.x = element_text(margin = margin(t = -10)))
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "bottom"
  ) + scale_size(range = c(-5, 10)) +
  annotate('rect', xmin=0, xmax=0.015, ymin=-Inf, ymax=Inf, alpha=.5, fill='red') +
  annotate('rect', xmin=0.016, xmax=0.06, ymin=-Inf, ymax=Inf, alpha=.5, fill='blue') +
  annotate('rect', xmin=0.13, xmax=0.17, ymin=-Inf, ymax=Inf, alpha=.3, fill='green') +
  annotate('rect', xmin=0.18, xmax=0.20, ymin=-Inf, ymax=Inf, alpha=.3, fill='purple')
#+ facet_grid(.~cluster)
#scale_y_continuous(limits = c(0,1000))
anvio_dist_vs_latitude_plot



## Adding clade info to the table
df_Vibrio$clade <- ifelse(df_Vibrio$S3320<0.1, "Clade_1", "Clade_2")
df_Vibrio$clade <- ifelse(df_Vibrio$S3320>0.18, "Clade_3", df_Vibrio$clade)



# df_Vibrio %>% filter(clade == "Clade_3") %>% ggplot(aes(bgcs_count, genome_id, fill=Taxonomy)) +
#   #ggplot(aes(abs(Latitude), as.numeric(bgcs_count), col=as.numeric(cluster))) +
#   geom_bar(position="stack", stat="identity") 


df_Vibrio$Spot_S.aur.
df_Vibrio %>% 
  ggplot(aes(Latitude, Spot_S.aur., col=Taxonomy)) +
  #ggplot(aes(abs(Latitude), as.numeric(bgcs_count), col=as.numeric(cluster))) +
  #geom_point(aes(size = bgcs_count)) + ylim(-80,80) +
  geom_point(aes(size = bgcs_count)) + scale_size(range = c(-5, 10))


df_Vibrio %>% filter(clade %in% c("Clade_4", "Clade_3", "Clade_2", "Clade_1") ) %>% 
  ggplot(aes(S3320, bgcs_count, col=Taxonomy)) +
  #ggplot(aes(abs(Latitude), as.numeric(bgcs_count), col=as.numeric(cluster))) +
  #geom_point(aes(size = bgcs_count)) + ylim(-80,80) +
  geom_point(aes(size = bgcs_count)) + scale_size(range = c(-5, 10))



anvio_dist_vs_latitude_plot_SpotS <-  
  df_Vibrio %>% 
  ggplot(aes(S3320, Latitude, col=S3320)) +
  #ggplot(aes(abs(Latitude), as.numeric(bgcs_count), col=as.numeric(cluster))) +
  #geom_point(aes(size = bgcs_count)) + ylim(-80,80) +
  geom_point(aes(size = Spot_S.aur.)) + ylim(-80,80) +
  #geom_point(size=2) + ylim(-80,80) +
  #geom_smooth() +
  #geom_smooth(method = lm, se=FALSE) +
  #stat_cor(method = "spearman", p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name="rho", size = 2)+
  #scale_size(range = c(0.1, 10)) +
  labs(x = "\nCore genes phylogenetic distance", y = "\nLatitude ") +
  theme_bw(base_size = 10) +
  scale_colour_gradient(low = "blue", high = "red", na.value = NA) +
  # scale_colour_steps2(low = "white",
  #                     mid = "blue",
  #                     high = "red",
  #                     midpoint = 3) + 
  #scale_size(range = c(-5, 10)) +
  #scale_fill_manual(values = c("#f25a13", "#580a82"))+
  #scale_color_manual(values = c("#f25a13", "#580a82"))+
  geom_label_repel(aes(label = labels),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme(#panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    #axis.title.x = element_text(margin = margin(t = -10)))
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "bottom"
  ) + scale_size(range = c(-5, 10))
#+ facet_grid(.~cluster)
#scale_y_continuous(limits = c(0,1000))
anvio_dist_vs_latitude_plot_SpotS

df_Vibrio %>% filter(S3320 <= 0.18 & bgcs_count >=6 & N50 >=100000) %>% dplyr::summarise(n())




###Filtering ######
df_Vibrio %>% filter(S3320 <= 0.18 & N50 >=100000) %>% 
  ggplot(aes(S3320, Latitude, col=S3320)) +
  #ggplot(aes(abs(Latitude), as.numeric(bgcs_count), col=as.numeric(cluster))) +
  #geom_point(aes(size = bgcs_count)) + ylim(-80,80) +
  geom_point(aes(size = N50)) +
  annotate('rect', xmin=0, xmax=0.015, ymin=-Inf, ymax=Inf, alpha=.5, fill='red') +
  annotate('rect', xmin=0.016, xmax=0.06, ymin=-Inf, ymax=Inf, alpha=.5, fill='blue') +
  annotate('rect', xmin=0.12, xmax=0.2, ymin=-Inf, ymax=Inf, alpha=.3, fill='green')

#For Pangenomic analyses pick representative:
# - relative close related with phylogenetic distance from root to max=0.18
# - pick those whitin lat -60 to -40: Southern
# - pick those wihtin lat 60 to 40: Northern 
# - pick those wihtin lat -20 to 20: Tropical  
df_Vibrio_Southern <- df_Vibrio %>% filter(S3320 <= 0.18 & bgcs_count >=6 & N50 >=100000 ) %>%
  filter(Latitude <= -40) #%>% dplyr::summarise(n()) # N=23
df_Vibrio_Southern$Placement <- rep("Southern")
df_Vibrio_Southern %>% dplyr::summarise(n()) # N=23

df_Vibrio_northern <- df_Vibrio %>% filter(S3320 <= 0.18 & bgcs_count >=6 & N50 >=100000 ) %>%
  filter(Latitude >= 40) #%>% dplyr::summarise(n()) # N=16
df_Vibrio_northern$Placement <- rep("Northern")
df_Vibrio_northern %>% dplyr::summarise(n()) # N=14


df_Vibrio_tropical <- df_Vibrio %>% filter(S3320 <= 0.18 & bgcs_count >=6 & N50 >=100000 ) %>%
  filter(Latitude <= 20 & Latitude >= -20) #%>% dplyr::summarise(n()) # N=51
df_Vibrio_tropical$Placement <- rep("Tropical")
df_Vibrio_tropical %>% dplyr::summarise(n()) # N=41


df_Vibrio_filter_pan <- rbind(df_Vibrio_Southern, df_Vibrio_northern, df_Vibrio_tropical)
df_Vibrio_filter_pan %>% summarise(n()) #78

df_Vibrio_filter_pan$Clade <- ifelse(df_Vibrio_filter_pan$S3320 < 0.015, "Clade_1", "Clade_2")
df_Vibrio_filter_pan$Clade <- ifelse(df_Vibrio_filter_pan$S3320 > 0.12, "Clade_3", df_Vibrio_filter_pan$Clade)

#Create list of names to substract from the external-genomes_refs.txt
list_of_names <- as.data.frame(df_Vibrio_filter_pan$genome_id)
write_csv(list_of_names, file = "list_of_names.txt")

#Create layers data to import to pan genome analyse in Anvio
misc_data_layers <- df_Vibrio_filter_pan[,c("genome_id", "Clade", "Placement", "Taxonomy")]
colnames(misc_data_layers)[1] <- "Item_Name"
write.table(misc_data_layers, file = "../Vibrio_filter_pan/misc_data_layers.txt", sep = "\t", quote = F, row.names = F)



 # WORLD mapping######
world_map <- map_data("world")

world_map


world_plot <- ggplot() +
  geom_map(
    data = world_map, map = world_map,
    aes(long, lat, map_id = region),
    color = "white", fill = "lightgray", size = 0.1
  ) + labs(x = "\nLongitude", y = "\nLatitude ") +
  geom_point(data = df_Vibrio
             , aes(Longitude, Latitude,  
                   color = S3320, 
                   #color = as.factor(Sub.sample), 
                   #size = Spot_S.aur.
             ), 
             #alpha = 0.5, 
             size=4
  ) +
  scale_colour_gradient(low = "blue", high = "red", na.value = NA) +
  # scale_colour_steps2(low = "white",
  #                     mid = "blue",
  #                     high = "red",
  #                     midpoint = 3) 
  # scale_colour_gradient2(low = "white",
  #                     mid = "blue",
  #                     high = "red",
  #                     midpoint = 4)+
  theme_bw(base_size = 10) +
  #scale_size(range = c(-5, 10)) + 
  theme(#panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    #axis.title.x = element_text(margin = margin(t = -10)),
    legend.position = "bottom"
  ) 
world_plot


library(cowplot)
world_plot_phylo <- plot_grid(world_plot, anvio_dist_vs_latitude_plot, ncol = 2, align = "h")
world_plot_phylo_SpotS <- plot_grid(world_plot, anvio_dist_vs_latitude_plot_SpotS, ncol = 2, align = "h")

anvio_dist_vs_latitude_plot_SpotS

ggsave(world_plot_phylo, file = "../Figures/world_plot_anvio_phylo.png", width = 9, height =5, limitsize = FALSE)
ggsave(world_plot_phylo_SpotS, file = "../Figures/world_plot_phylo_SpotS.png", width = 9, height =5, limitsize = FALSE)


#### Import anvi-compute-functional-enrichment-in-pan results
anvio_fun_enrich_placement_antiSMASH <- read.table(file = "../enriched-functions_placement_antiSMASH.txt", header = TRUE, sep = "\t")

anvio_fun_enrich_placement_antiSMASH %>% filter(enrichment_score >= 20 & associated_groups == "Tropical")


## Import summary from pan genome and looks for SGCs BGCs


library("data.table") 
#Vibrio_filter_pan_gene_clusters_summary <- read.table(file = "../../../Downloads/Vibrio_filter_pan_gene_clusters_summary.txt", header = FALSE, sep = "\t", row.names = F)
Vibrio_filter_pan_gene_clusters_summary <- fread(file = "../../../Downloads/Vibrio_filter_pan_gene_clusters_summary.txt", sep = "\t", select = c(1:22) ) 
head(Vibrio_filter_pan_gene_clusters_summary)
Vibrio_filter_pan_gene_clusters_summary <- as.data.frame(Vibrio_filter_pan_gene_clusters_summary)

SCG_BGC_genes <- Vibrio_filter_pan_gene_clusters_summary %>% filter(SCG==1, antiSMASH!="")
colnames(SCG_BGC_genes)[5] <- "gene_callers_id_1"

### Import antiSMASH function table and combine them to a whole:
library(dplyr)
library(readr)

# Get the list of files matching the pattern "*mod*"
file_names <- list.files(path = "../anvio/gbk_bgcflow_antismash/functions", pattern = "*mod*", full.names = T)
#file_names <- head(file_names)

# Create an empty list to store dataframes
data_list <- list()

# Loop through each file
for (file_name in file_names) {
  # Read the file into a dataframe
  df <- read_tsv(file_name, col_names = F, skip = 1)
  
  # Change the header of the 6th column to "contig_id"
  colnames(df) <- c("gene_callers_id","annotation","antiSMASH_ACC","function","e-val","contig_id")
  
  # Add the dataframe to the list
  data_list[[file_name]] <- df
}


# Merge all dataframes into a single dataframe
merged_df <- as.data.frame(bind_rows(data_list, .id = "file_name"))

# split last column in genome_id and contig_id
merged_df <- merged_df %>%
  separate(contig_id, into = c("genome_id_1", "contig_id"), "_c_") %>%
  mutate(contig_id = paste0("c_", contig_id))

merged_df <- merged_df[,-1]


## import bgc overview made from the bgcflow pipeline
df_regions_antismash_6.1.1 <- read.csv(file = "../bgcflow_antismash_gbk_regions/df_regions_antismash_6.1.1.csv")
colnames(df_regions_antismash_6.1.1)[4] <- "contig_id"

head(df_regions_antismash_6.1.1)

head(SCG_BGC_genes)

head(merged_df)

SCG_BGC_genes_antiSMASH_hits <- left_join(SCG_BGC_genes, merged_df) %>% filter(!is.na(contig_id)) 

SCG_BGC_genes_antiSMASH_hits <- SCG_BGC_genes_antiSMASH_hits %>% select(gene_cluster_id, genome_name, num_genomes_gene_cluster_has_hits, functional_homogeneity_index,
                                        antiSMASH_ACC, antiSMASH, contig_id)

df_regions_antismash_6.1.1_selected <- df_regions_antismash_6.1.1 %>% select(contig_id, product, most_similar_known_cluster_id, most_similar_known_cluster_description, similarity)

SCG_BGC_genes_antiSMASH_hits_descp <- left_join(SCG_BGC_genes_antiSMASH_hits,df_regions_antismash_6.1.1_selected, by="contig_id") #%>% 

SCG_BGC_genes_antiSMASH_hits_descp <- SCG_BGC_genes_antiSMASH_hits_descp[!duplicated(SCG_BGC_genes_antiSMASH_hits_descp$gene_cluster_id),] 

SGCs_antiSMASH_genes <- SCG_BGC_genes_antiSMASH_hits_descp %>% group_by(product, most_similar_known_cluster_description, most_similar_known_cluster_id) %>% summarise(n())
write.table(SGCs_antiSMASH_genes, "Figures/SGCs_antiSMASH_genes.txt", sep = "\t")
summarise(SCG_BGC_genes_antiSMASH_hits_descp, )

dim(SCG_BGC_genes_antiSMASH_hits_descp)


