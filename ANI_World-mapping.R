library(tidyverse)
library(dplyr)
library(stringr)
library(ggpubr)

ANI <- read.table(file = "../ANI/fastANI-out-raw.tsv", header = FALSE, sep = "\t")

#Clean up ANI table

colnames(ANI) <- c("query","reference","ANI","num_frags_mapped","total_query_frags")

ANI$query <-  gsub("\\.fna", "", ANI$query)
ANI$query <-  gsub("\\_genomic.fna.gz", "", ANI$query)
ANI$query <-  gsub("\\_genomic.gz", "", ANI$query)
ANI$reference <-  gsub("\\.fna", "", ANI$reference)
ANI$reference <-  gsub("\\_genomic.fna.gz", "", ANI$reference)
ANI$reference <-  gsub("\\_genomic.gz", "", ANI$reference)

ANI_new_ref <- separate(ANI, reference, into = c("Assembly.Accession", "Assembly.Name"), sep = "_ASM")
ANI_new_ref <- ANI_new_ref[1:3]
ANI_new_ref$Assembly.Name <- ifelse(!is.na(ANI_new_ref$Assembly.Name), paste0("ASM", ANI_new_ref$Assembly.Name), ANI_new_ref$Assembly.Name)

#Add taxonomic informaton:
ref_tax <- read.table(file = "../Metadata/ncbi_dataset.tsv", header = TRUE, sep = "\t")
ref_tax$Organism.Name
ref_tax$Assembly.Accession
ref_tax$Assembly.Name

ANI_new_ref <- left_join(ANI_new_ref, ref_tax)

# Clean up metadata table

Meta_data_Galathea <- read.table(file = "../Metadata/Meta_data_latitude.txt", header = TRUE, sep = "\t")


str(df)
Meta_data_Galathea_Vibrios <- Meta_data_Galathea %>% filter(Strain %in% ANI$query)

Meta_data_Galathea_Vibrios$X16SrRNA_gene <- gsub("\\s+$", "", Meta_data_Galathea_Vibrios$X16SrRNA_gene)


#Convert antibiotic activity to 1 to 3 scale
Meta_data_Galathea_Vibrios <- Meta_data_Galathea_Vibrios[,-c(1, ncol(Meta_data_Galathea_Vibrios))]
colnames(Meta_data_Galathea_Vibrios)[1] <- "samples"

Meta_data_Galathea_Vibrios[] <- lapply(Meta_data_Galathea_Vibrios, function(x) gsub("[()]", "", x))

unique(Meta_data_Galathea_Vibrios$Spot_S.aur.)


lreplace_plus_values <- function(x) {
  ifelse(x == "-", 0,
         ifelse(x == "+", 1, 
                ifelse(x == "++", 2, 
                       ifelse(x == "++++", 4, 
                       ifelse(x == "+++", 3, x)))))
}

# Use mutate_all() to apply the function to all columns
Meta_data_Galathea_Vibrios <- Meta_data_Galathea_Vibrios %>% mutate_all(replace_plus_values)

write.table(Meta_data_Galathea_Vibrios, "../Metadata/Meta_data_Galathea_Vibrios_metadata.txt", col.names = TRUE, row.names = F, quote = FALSE, sep = "\t")


unique(Meta_data_Galathea_Vibrios$Spot_S.aur.)

# Investigate ANI clusters with Pheatmap ####
library(pheatmap)

#Convert fastANI format to matrix

# pivot wider the ANI dataframe
ANI_matrix <- pivot_wider(ANI[1:3], names_from = reference, values_from = ANI)
ANI_matrix <- as.data.frame(ANI_matrix)

rownames(ANI_matrix) <- ANI_matrix$query
ANI_matrix <- ANI_matrix[,-1]

ANI_matrix <- ANI_matrix[,row.names(ANI_matrix)]
max(ANI_matrix)
ANI_matrix[is.na(ANI_matrix)] <- 75

#pheatmap_ANI <- pheatmap(ANI_matrix)

#Cluster according to hclust
Vibrio_hclust <- hclust(dist(ANI_matrix), method = "ward.D2")
dist(ANI_matrix)
plot(as.dendrogram(Vibrio_hclust))
#rect.hclust(Vibrio_hclust, k = 6)
rect.hclust(Vibrio_hclust, k = 6,
            which = c(1,2,3,4,5))

Vibrio_hclust$height
Vibrio_hclust$merge

# load package
library(dendextend)

as.dendrogram(Vibrio_hclust) %>%
  plot(horiz = TRUE)


#Vibrio_clusters_10 <- data.frame(cluster = cutree(pheatmap_ANI$tree_row, k = 6))
Vibrio_clusters_10 <- data.frame(cutree(tree = as.dendrogram(Vibrio_hclust), k = 6))
colnames(Vibrio_clusters_10)[1] <- "cluster"

#Add sidebar to pheatmap

Vibrio_clusters_10$cluster <- as.character(Vibrio_clusters_10$cluster)
ref_tax$reference <- paste0(ref_tax$Assembly.Accession, sep = "_", ref_tax$Assembly.Name)

Vibrio_clusters_10 <- Vibrio_clusters_10 %>% mutate(reference = rownames(Vibrio_clusters_10)) 

Vibrio_clusters_10 <- left_join(Vibrio_clusters_10, ref_tax)
Meta_data_Galathea_Vibrios <- mutate(Meta_data_Galathea_Vibrios, reference = Meta_data_Galathea_Vibrios$samples)

Vibrio_clusters_10 <-left_join(Vibrio_clusters_10, Meta_data_Galathea_Vibrios)


rownames(Vibrio_clusters_10) <- Vibrio_clusters_10$reference
colnames(Vibrio_clusters_10)
Vibrio_clusters_10 <- Vibrio_clusters_10[,c(1, 5, 19)]

Vibrio_clusters_10$Ref_strain <- ifelse(is.na(Vibrio_clusters_10$Organism.Name),
                                        "0", 1)

#Colors for metadata layer
colnames(Vibrio_clusters_10)
meta_data_layer_pheatmap <- Vibrio_clusters_10[,c(1,3,4)]
meta_data_layer_pheatmap$cluster <- factor(meta_data_layer_pheatmap$cluster, levels = unique(meta_data_layer_pheatmap$cluster))
meta_data_layer_pheatmap$Type <-  factor(meta_data_layer_pheatmap$Type, levels = unique(meta_data_layer_pheatmap$Type))

# Purplerain <- c("#e0b62b","#f25a13", "#2b55e0", "#240785", "#9111ab", "#e33963", "#e35b39", "#e3b039",
#                 "#580a82", "#060896", "#0d4582", "#c439e3", "#4d44ab")

Purplerain <- c("black","#f25a13", "#2b55e0", "#240785", "#e33963", "#e35b39", "#e3b039",
                 "#c439e3", "#4d44ab")
scales::show_col(Purplerain)

extra_col_4 <- c("#8e0152","#c51b7d","#de77ae","#f1b6da")

library(RColorBrewer)
Paired <- brewer.pal(11, "Paired")

Colors_22 <- c(Purplerain, Paired,extra_col_4)
Colors_2 <- c("white", "#240785")
Colors_6 <- c("black","#f25a13", "#2b55e0", "#240785", "#e33963", "#e35b39","#e3b039")
scales::show_col(Colors_22)
names(Colors_22) <- unique(meta_data_layer_pheatmap$cluster)
names(Colors_2) <-  unique(meta_data_layer_pheatmap$Ref_strain)
names(Colors_6) <-  unique(meta_data_layer_pheatmap$Type)
Colors_pheatmap <-list(cluster=Colors_22, Type=Colors_6, Ref_strain=Colors_2)

meta_data_layer_pheatmap$Type[5] <- "Seawater"


#Our results reveal clear genetic discontinuity, with 99.8% 
#of the total 8 billion genome pairs analyzed conforming to 
#>95% intra-species and <83% inter-species ANI values

ANI_matrix_pheatmap <- pheatmap(ANI_matrix, 
         show_colnames     = FALSE,
         show_rownames     = FALSE,
         #show_rownames=T, cluster_cols=T, cluster_rows=T, scale="row",
         #cex=1, 
         #clustering_distance_rows="canberra", cex=1,
         #clustering_distance_cols="euclidean", 
         clustering_method="ward.D2", cutree_cols = 6,
         #border_color=FALSE,
         annotation_col=meta_data_layer_pheatmap, 
         #annotation_colors = Colors_pheatmap
         )
ANI_matrix_pheatmap

ggsave(ANI_matrix_pheatmap, file = "../../ANI_matrix_pheatmap.png", width = 10, height = 10)


# CAZyme data frame######
#Load CAZyme info
Vibrio_CAZymes_table <- readRDS(file = "Vibrio_CAZymes_table.rds")
rownames(Vibrio_CAZymes_table) <- Vibrio_CAZymes_table$CAZyme
Vibrio_CAZymes_table <- Vibrio_CAZymes_table[,-1]
Vibrio_CAZymes_Summery <- readRDS(file = "Vibrio_CAZymes_Summery.rds")
rownames(Vibrio_CAZymes_Summery) <- Vibrio_CAZymes_Summery$CAZyme_fam
Vibrio_CAZymes_Summery <- Vibrio_CAZymes_Summery[,-1]
#Vibrio_CAZymes_table_long <- pivot_longer(Vibrio_CAZymes_table, cols = -CAZyme, names_to = "samples", values_to = "Value")
#Vibrio_CAZymes_Summery_long <- pivot_longer(Vibrio_CAZymes_Summery, cols = -CAZyme_fam, names_to = "samples", values_to = "Value")

Vibrio_CAZymes_table <- as.data.frame(t(Vibrio_CAZymes_table))
str(Vibrio_CAZymes_table)
rownames(Vibrio_CAZymes_table)
rownames(Vibrio_CAZymes_table) <- gsub("\\_Strain", "", rownames(Vibrio_CAZymes_table))
Vibrio_CAZymes_table$samples <- rownames(Vibrio_CAZymes_table)
Vibrio_CAZymes_table[is.na(Vibrio_CAZymes_table)] <- 0



Vibrio_CAZymes_Summery <- as.data.frame(t(Vibrio_CAZymes_Summery))
str(Vibrio_CAZymes_Summery)
rownames(Vibrio_CAZymes_Summery)
rownames(Vibrio_CAZymes_Summery) <- gsub("\\_summary", "", rownames(Vibrio_CAZymes_Summery))
Vibrio_CAZymes_Summery$samples <- rownames(Vibrio_CAZymes_Summery)
Vibrio_CAZymes_Summery[is.na(Vibrio_CAZymes_Summery)] <- 0



Vibrio_clusters_10
rownames(Meta_data_Galathea_Vibrios) <- Meta_data_Galathea_Vibrios$samples
colnames(Meta_data_Galathea_Vibrios)
Vibrio_clusters_10 <- Vibrio_clusters_10 %>% mutate(samples = rownames(Vibrio_clusters_10))
ANI_mat_long_meta <- left_join(Meta_data_Galathea_Vibrios, Vibrio_clusters_10)

ANI_mat_long_meta_cazy <- left_join(ANI_mat_long_meta, Vibrio_CAZymes_Summery)

#Add Cazy sums to data frame
ANI_mat_long_meta_cazy$Cazy_sum <- rowSums(ANI_mat_long_meta_cazy[15:23])

colnames(ANI_mat_long_meta_cazy)


ANI_mat_long_meta_cazy <- left_join(ANI_mat_long_meta_cazy, Vibrio_CAZymes_table)
max(ANI_mat_long_meta_cazy$Latitude)

#ANI_mat_long_meta_cazy[4:5] <- ANI_mat_long_meta_cazy[,c(5,4)] 
ANI_mat_long_meta_cazy$Latitude <- as.numeric(ANI_mat_long_meta_cazy$Latitude)
ANI_mat_long_meta_cazy$Longitude <- as.numeric(ANI_mat_long_meta_cazy$Longitude)
ANI_mat_long_meta_cazy$Spot_V.ang <- as.numeric(ANI_mat_long_meta_cazy$Spot_V.ang)
ANI_mat_long_meta_cazy$Spot_S.aur. <- as.numeric(ANI_mat_long_meta_cazy$Spot_S.aur.)
#ANI_mat_long_meta_cazy$GH <- as.numeric(ANI_mat_long_meta_cazy$GH)
str(ANI_mat_long_meta_cazy)



# antiSMASH data frame ########

df_antismash <- read.table(file = "../bgcflow/antiSMASH/df_antismash_6.1.1_summary.csv", header = TRUE, sep = ",")
df_antismash$bgcs_count


colnames(df_antismash)[1] <- "samples"
ANI_mat_long_meta_cazy_antismash <- left_join(ANI_mat_long_meta_cazy, df_antismash)
# genome info data frame ########
df_seqfu_stats <- read.table(file = "../bgcflow/antiSMASH/df_seqfu_stats.csv", header = TRUE, sep = ",")

colnames(df_seqfu_stats)[1] <- "samples"
ANI_mat_long_meta_cazy_antismash_seqfu <- left_join(ANI_mat_long_meta_cazy_antismash, df_seqfu_stats)

colnames(df_seqfu_stats)
ANI_mat_long_meta_cazy_antismash_seqfu$Latitude
ANI_mat_long_meta_cazy_antismash_seqfu$mmZone_V.ang

# Correlation patterns #####


#Investigate if rich (unique clusters per latitue is correlated)

ANI_mat_long_meta_cazy_antismash_seqfu %>% select(Latitude, cluster) %>% group_by(Latitude) %>% 
  dplyr::summarise(richness = n_distinct(cluster)) %>%
  #ggplot(aes(abs(Latitude), as.numeric(cluster), col=as.numeric(cluster))) +
  ggplot(aes(abs(Latitude), richness)) + ylim(0,6) +
  geom_point(size = 2) +
  geom_smooth(method = lm)

#Higher richness around -10 latitude
ANI_mat_long_meta_cazy_antismash_seqfu %>% select(Latitude, cluster) %>% group_by(Latitude) %>% 
  dplyr::summarise(richness = n_distinct(cluster)) %>%
ggplot(aes(x=(Latitude))) + 
  #geom_density()
  geom_histogram(binwidth=10)


#Import autoMLST distance matrix and use S3704
autoMLST_dist <- read.table(file = "../bgcflow/automlst_dist.txt", header = TRUE,  sep="\t")
autoMLST_dist_S3704 <- as.data.frame(autoMLST_dist$S3704)
rownames(autoMLST_dist_S3704) <- colnames(autoMLST_dist)[-1]

autoMLST_dist_S3704 <- autoMLST_dist_S3704 %>% mutate(samples = rownames(autoMLST_dist_S3704))
colnames(autoMLST_dist_S3704)[1] <- "S3704"
ANI_mat_long_meta_cazy_antismash_seqfu <- left_join(ANI_mat_long_meta_cazy_antismash_seqfu, autoMLST_dist_S3704)

autoMLST_dist_S3320 <- as.data.frame(autoMLST_dist$S3320)
rownames(autoMLST_dist_S3320) <- colnames(autoMLST_dist)[-1]
autoMLST_dist_S3320 <- autoMLST_dist_S3320 %>% mutate(samples = rownames(autoMLST_dist_S3320))
colnames(autoMLST_dist_S3320)[1] <- "S3320"
ANI_mat_long_meta_cazy_antismash_seqfu <- left_join(ANI_mat_long_meta_cazy_antismash_seqfu, autoMLST_dist_S3320)


autoMLST_dist_S1106 <- as.data.frame(autoMLST_dist$S1106)
rownames(autoMLST_dist_S1106) <- colnames(autoMLST_dist)[-1]
autoMLST_dist_S1106 <- autoMLST_dist_S1106 %>% mutate(samples = rownames(autoMLST_dist_S1106))
colnames(autoMLST_dist_S1106)[1] <- "S1106"
ANI_mat_long_meta_cazy_antismash_seqfu <- left_join(ANI_mat_long_meta_cazy_antismash_seqfu, autoMLST_dist_S1106)

ANI_mat_long_meta_cazy_antismash_seqfu$labels <- ifelse(ANI_mat_long_meta_cazy_antismash_seqfu$samples=="S3320", 
                                                        "S3320", NA)

str(ANI_mat_long_meta_cazy_antismash_seqfu$Ref_strain)
library(ggrepel)

ANI_mat_long_meta_cazy_antismash_seqfu$Spot_V.ang
automlst_dist_vs_latitude_plot <-  
  ANI_mat_long_meta_cazy_antismash_seqfu %>% filter(Ref_strain == 0) %>% 
  ggplot(aes(S3320, Latitude, col=S3320)) +
  #ggplot(aes(abs(Latitude), as.numeric(bgcs_count), col=as.numeric(cluster))) +
  #geom_point(aes(size = bgcs_count)) + ylim(-80,80) +
  geom_point(aes(size = Spot_S.aur.)) + ylim(-80,80) +
  #geom_point(size=2) + ylim(-80,80) +
  #geom_smooth() +
  #geom_smooth(method = lm, se=FALSE) +
  #stat_cor(method = "spearman", p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name="rho", size = 2)+
  #scale_size(range = c(0.1, 10)) +
  labs(x = "\nautomlst phylogenetic distance", y = "\nLatitude ") +
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
automlst_dist_vs_latitude_plot


# automlst_dist_vs_latitude_plot_cazy <-  
#   ANI_mat_long_meta_cazy_antismash_seqfu %>% filter(Ref_strain == 0) %>% 
#   ggplot(aes(S3320, Latitude, col=S3320)) +
#   #ggplot(aes(abs(Latitude), as.numeric(bgcs_count), col=as.numeric(cluster))) +
#   geom_point(aes(size = PL)) + ylim(-80,80) +
#   #geom_smooth() +
#   #geom_smooth(method = lm, se=FALSE) +
#   #stat_cor(method = "spearman", p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name="rho", size = 2)+
#   #scale_size(range = c(0.1, 10)) +
#   labs(x = "\nautomlst phylogenetic distance", y = "\nLatitude ") +
#   theme_bw(base_size = 10) +
#   scale_colour_gradient(low = "blue", high = "red", na.value = NA) +
#   # scale_colour_steps2(low = "white",
#   #                     mid = "blue",
#   #                     high = "red",
#   #                     midpoint = 3) + 
#   #scale_size(range = c(-5, 10)) +
#   #scale_fill_manual(values = c("#f25a13", "#580a82"))+
#   #scale_color_manual(values = c("#f25a13", "#580a82"))+
#   geom_label_repel(aes(label = labels),
#                    box.padding   = 0.35, 
#                    point.padding = 0.5,
#                    segment.color = 'grey50') +
#   theme(#panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.title = element_text(face = "bold"),
#     axis.text = element_text(face = "bold"),
#     #axis.title.x = element_text(margin = margin(t = -10)))
#     axis.title.y = element_blank(),
#     axis.text.y = element_blank(),
#     legend.position = "bottom"
#   ) + scale_size(range = c(-10, 10))
# #+ facet_grid(.~cluster)
# #scale_y_continuous(limits = c(0,1000))
# automlst_dist_vs_latitude_plot_cazy

#ggsave(Cluster_vs_latitude_plot, file = "../../Cluster_vs_latitude_plot.png", width = 5, height =5, limitsize = FALSE)

unique(ANI_mat_long_meta_cazy_antismash_seqfu$Type)
ANI_mat_long_meta_cazy_antismash_seqfu %>% 
  filter(Type == "Seawater CTD") %>% group_by(Sub.sample) %>% dplyr::summarise(n()) %>% dplyr::summarise(n())

ANI_mat_long_meta_cazy_antismash_seqfu %>% 
  filter(Type != "Seawater CTD") %>% group_by(Sub.sample) %>% dplyr::summarise(n()) %>% dplyr::summarise(n())


# Other correlation patterns ######
#colnames(ANI_mat_long_meta_cazy_antismash_seqfu)

#Cazy_sum_vs_Genome_size_plot <- 
  #ANI_mat_long_meta_cazy_antismash_seqfu %>% 
  #ggplot(aes(abs(Latitude), bgcs_count, col=as.numeric(cluster), group = cluster)) +
  #ggplot(aes(bgcs_count, Total, col=as.numeric(cluster), group = cluster)) +
  #ggplot(aes(Cazy_sum, Total, col=as.numeric(cluster), group = cluster)) +
  #ggplot(aes(abs(Latitude), PL, col=as.numeric(cluster), group = cluster)) +
  #ggplot(aes(abs(Latitude), as.numeric(cluster), col=as.numeric(cluster), group = cluster)) +
  #ggplot(aes(bgcs_count, bgcs_on_contig_edge, col=as.numeric(cluster), group = cluster)) +
  #geom_point(size = 2) +
  #geom_smooth(method = lm, se=FALSE)+
  #stat_cor(method = "spearman", p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name="rho", size = 2)+
  #scale_size(range = c(0.1, 10)) +
  #labs(x = "\n\n16S richness ", y = "\nAD richness ", subtitle = "") +
  #theme_bw(base_size = 10) +
  #scale_color_gradient2(low = "yellow", mid="blue", high = "red", na.value = NA)
  # scale_colour_steps2(low = "white",
  #                     mid = "blue",
  #                     high = "red",
  #                     midpoint = 3) + 
  #scale_size(range = c(-5, 10)) +
  #scale_fill_manual(values = c("#f25a13", "#580a82"))+
  #scale_color_manual(values = c("#f25a13", "#580a82"))+
  # theme(panel.grid.major = element_blank(),
  #       panel.grid.minor = element_blank(),
  #       axis.title = element_text(face = "bold"),
  #       axis.text = element_text(face = "bold"),
  #       axis.title.x = element_text(margin = margin(t = -10)),
  #       legend.position = "none") 
#+ facet_grid(.~cluster)
#scale_y_continuous(limits = c(0,1000))

#ggsave(Cazy_sum_vs_Genome_size_plot, file = "../../Cazy_sum_vs_Genome_size_plot.png", width = 5, height =5, limitsize = FALSE)


df <- ANI_mat_long_meta_cazy_antismash_seqfu

library(Hmisc)
colnames(df %>% select(where(is.numeric)) %>% select(!N90))
rownames(df) <- df$samples
#only select cols with numeric values
df_num <- df %>% select(where(is.numeric)) 

# Remove columns with fewer than 3 non-missing observations
df_num <- df_num[, colSums(!is.na(df_num)) >= 4]

df_num[is.na(df_num)] <- 0


# Remove columns with sd less than 3
df_num <- df_num %>% 
  select_if(~ sd(.) >= 3)

df_num$Latitude <- abs(df_num$Latitude)
df_corr <- rcorr(as.matrix(df_num))

str(df_corr)



# Convert matrix to long format
library(reshape2)
df_corr_r <- melt(df_corr$r)
df_corr_n <- melt(df_corr$n)
df_corr_P <- melt(df_corr$P)

df_corr_long <- cbind(df_corr_r, n=df_corr_n[,3], P=df_corr_P[,3])
colnames(df_corr_long)[3] <- "r"

df_corr_long_filt <- df_corr_long %>% 
                    filter(r!=0, r!=-Inf, r!=Inf, 
                           n>3, P<=0.01) %>% distinct()

df_corr_long_filt %>% filter(Var1 == "Latitude")






# Richness temperate vs. tropical climate ####
# Devide genomes in 23.5° to 66.8° N/S of Equator
max(ANI_mat_long_meta_cazy_antismash_seqfu$Latitude)
max(ANI_mat_long_meta_cazy_antismash_seqfu$Latitude)

ANI_mat_long_meta_cazy_antismash_seqfu$Climate <- 
  ifelse(abs(ANI_mat_long_meta_cazy_antismash_seqfu$Latitude)>23.5, "Temperate", "Tropic")

ANI_mat_long_meta_cazy_antismash_seqfu %>% group_by(Climate) %>% 
  dplyr::summarise(n())

library(dplyr)
library(broom)

df <- ANI_mat_long_meta_cazy_antismash_seqfu[,-c(717:724, 726:727)]

colnames(df)
library(Hmisc)
rownames(df) <- df$samples
#only select cols with numeric values
df_num <- df %>% select(where(is.numeric)) 

# Remove columns with fewer than 4 non-missing observations
df_num <- df_num[, colSums(!is.na(df_num)) >= 4]

df_num[is.na(df_num)] <- 0

df_num <- apply(df_num, 2, as.numeric)
rownames(df_num) <- df$samples

df_num <- as.data.frame(df_num)

df_num$Climate <- df$Climate
unique(df_num$Climate)
str(df_num)


df_num_t_test <- df_num %>% 
  #select(tip, total_bill, sex) %>% 
  gather(key = variable, value = value, -Climate) %>% 
  group_by(Climate, variable) %>% 
  summarise(value = list(value)) %>% 
  spread(Climate, value) %>% 
  group_by(variable) %>% 
  mutate(p_value = t.test(unlist(Temperate), unlist(Tropic))$p.value,
         t_value = t.test(unlist(Temperate), unlist(Tropic))$statistic) %>%
  mutate(padj = p.adjust(p_value, "BH")) %>% filter(padj <= 0.01 )


df_num_t_test

a <- data.frame(variable=colnames(Vibrio_CAZymes_table), group=rep("Cazymes_subfam"))
b <- data.frame(variable=colnames(Vibrio_CAZymes_Summery), group=rep("Cazymes_fam"))
c <- data.frame(variable=colnames(df_antismash), group=rep("antiSMASH_prediction"))

d <- rbind(a, b, c)

colnames(d)

df_num_t_test <- left_join(df_num_t_test, d)

df_num_t_test_names <- df_num_t_test %>% filter(group %in% c("Cazymes_fam", "antiSMASH_prediction"))

df_num_t_test_names <- df_num_t_test_names$variable

# WORLD mapping######
world_map <- map_data("world")

world_map

# colnames(ANI_mat_long_meta_cazy_antismash_seqfu)
# 
# ANI_mat_long_meta_cazy_antismash_seqfu %>% filter(Spot_S.aur.>2) %>% select(Sub.sample)
# unique(ANI_mat_long_meta_cazy_antismash_seqfu$Spot_S.aur.)

world_plot <- ggplot() +
  geom_map(
    data = world_map, map = world_map,
    aes(long, lat, map_id = region),
    color = "white", fill = "lightgray", size = 0.1
  ) + labs(x = "\nLongitude", y = "\nLatitude ") +
  geom_point(data = ANI_mat_long_meta_cazy_antismash_seqfu
             , aes(Longitude, Latitude,  
                                           color = S3320, 
                                           #color = as.factor(Sub.sample), 
                                           size = Spot_S.aur.
                                           ), 
             #alpha = 0.5, 
             #size=2
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
  scale_size(range = c(-5, 10)) + theme(#panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank(),
                                        axis.title = element_text(face = "bold"),
                                        axis.text = element_text(face = "bold"),
                                        #axis.title.x = element_text(margin = margin(t = -10)),
                                        legend.position = "bottom"
                                        ) 
world_plot



library(cowplot)
world_plot_automlst <- plot_grid(world_plot, automlst_dist_vs_latitude_plot, ncol = 2, align = "h")
world_plot_automlst_noSize <- plot_grid(world_plot, automlst_dist_vs_latitude_plot, ncol = 2, align = "h")
world_plot_automlst_Spot_S.aur <- plot_grid(world_plot, automlst_dist_vs_latitude_plot, ncol = 2, align = "h")


ggsave(world_plot_automlst, file = "../../world_plot_automlst.png", width = 9, height =5, limitsize = FALSE)
ggsave(world_plot_automlst_noSize, file = "../../world_plot_automlst_noSize.png", width = 9, height =5, limitsize = FALSE)
ggsave(world_plot_automlst_Spot_S.aur, file = "../../world_plot_automlst_Spot_S.aur.png", width = 9, height =5, limitsize = FALSE)


library(cowplot)

# CAZymes superfamily distributions
# plots <- lapply(14:20, function(i) {
#   ggplot() +
#     geom_map(data = world_map, map = world_map,
#       aes(long, lat, map_id = region),
#       color = "white", fill = "lightgray", size = 0.1) + 
#  geom_point(data = ANI_mat_long_meta_cazy_antismash_seqfu, aes(Longitude, Latitude,  
#  color = as.numeric(cluster), size = (ANI_mat_long_meta_cazy_antismash_seqfu[[i]])),
# alpha = 0.5) + scale_colour_steps2(low = "black",
#                                                   mid = "blue",
#                                                   high = "red",
#                                                   midpoint = 3) + 
#     scale_size(range = c(-5, 10)) +
#     labs(size = colnames(ANI_mat_long_meta_cazy_antismash_seqfu)[i]) # set the plot title to the column name
# })
# 
# 
# CAZy_summary_plots <- plot_grid(plotlist = plots, nrow = 2, ncol = 3)
# 
# ggsave(CAZy_summary_plots, file = "../../CAZy_summary_plots.png", width = 20, height =10, limitsize = FALSE)

# Individual hmm families

# plots <- lapply(which(colnames(ANI_mat_long_meta_cazy_antismash_seqfu) %in% df_num_t_test_names)
# , function(i) {
#   ggplot() +
#     geom_map(data = world_map, map = world_map,
#              aes(long, lat, map_id = region),
#              color = "white", fill = "lightgray", size = 0.1) + 
#     geom_point(data = ANI_mat_long_meta_cazy, aes(Longitude, Latitude,  
#                                                   color = as.numeric(cluster), size = ANI_mat_long_meta_cazy_antismash_seqfu[[i]]),
#                alpha = 0.5) + scale_colour_steps2(low = "black",
#                                                   mid = "blue",
#                                                   high = "red",
#                                                   midpoint = 3) + 
#     scale_size(range = c(0, 10)) +
#     labs(size = colnames(ANI_mat_long_meta_cazy_antismash_seqfu)[i]) # set the plot title to the column name
# })


# plots <- lapply(grep("^PL", colnames(ANI_mat_long_meta_cazy)), function(i) {
#   ggplot() +
#     geom_map(data = world_map, map = world_map,
#              aes(long, lat, map_id = region),
#              color = "white", fill = "lightgray", size = 0.1) + 
#     geom_point(data = ANI_mat_long_meta_cazy, aes(Longitude, Latitude,  
#                                                   color = as.numeric(cluster), size = ANI_mat_long_meta_cazy[[i]]),
#                alpha = 0.5) + scale_colour_steps2(low = "black",
#                                                   mid = "blue",
#                                                   high = "red",
#                                                   midpoint = 3) + 
#     scale_size(range = c(0, 10)) +
#     labs(size = colnames(ANI_mat_long_meta_cazy)[i]) # set the plot title to the column name
# })

length(plots)

# combine the plots into a single figure
Sig_latitude_plots <- plot_grid(plotlist = plots, nrow = 5, ncol = 3)

ggsave(Sig_latitude_plots, file = "../../Sig_latitude_plots.png", width = 28, height =18, limitsize = FALSE)




# Metadata file for anvio######
write.table(Vibrio_clusters_10, "../Metadata/Vibrio_clusters_10.txt", col.names = TRUE, row.names = F, quote = FALSE, sep = "\t")


# Specific BGCs####
df_regions_antismash_6.1.1 <- read.table(file = "../bgcflow/tables/df_regions_antismash_6.1.1_corrected.txt",  header = T , sep="\t")

df_regions_antismash_6.1.1$region <- gsub("\\_", ".", df_regions_antismash_6.1.1$region )
df_regions_antismash_6.1.1$similarity <- as.numeric(gsub("\\_", ".", df_regions_antismash_6.1.1$similarity))
df_regions_antismash_6.1.1$most_similar_known_cluster_type
df_regions_antismash_6.1.1$product_1
colnames(df_regions_antismash_6.1.1)[2] <- "samples"

#Summarize described BGCs
# df_regions_antismash_6.1.1 %>% filter(most_similar_known_cluster_description!="") %>%
#   group_by(genome_id, most_similar_known_cluster_type, most_similar_known_cluster_description) %>%
#   dplyr::summarise(sum=n()) 


df_regions_antismash_6.1.1_meta <- df_regions_antismash_6.1.1 %>% #filter(contig_edge == FALSE) %>%
  group_by(samples, product_1) %>%
  dplyr::summarise(count=n()) %>% left_join(., ANI_mat_long_meta_cazy_antismash_seqfu)

df_regions_antismash_6.1.1_meta_id <- df_regions_antismash_6.1.1 %>% #filter(contig_edge == FALSE) %>%
  group_by(samples, most_similar_known_cluster_description) %>%
  dplyr::summarise(count=n()) %>% left_join(., ANI_mat_long_meta_cazy_antismash_seqfu)



df_regions_antismash_6.1.1_meta %>% group_by(product_1) %>% filter(n() > 8) %>% 
  ggplot(aes(fill=product_1, x=reorder(samples,S3320), y=count)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_fill_brewer(palette = "Spectral")

df_regions_antismash_6.1.1_meta_id %>% group_by(most_similar_known_cluster_description) %>% filter(n() > 5) %>% 
  ggplot(aes(fill=most_similar_known_cluster_description, x=reorder(samples,S3320), y=count)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_fill_brewer(palette = "Spectral")

df_regions_antismash_6.1.1_meta_id %>% group_by(most_similar_known_cluster_id) %>% filter(n() > 5) %>% 
  ggplot(aes(fill=most_similar_known_cluster_id, x=reorder(samples,S3320), y=count)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_fill_brewer(palette = "Spectral")


bar_bgcs_cosmo_product_1 <- df_regions_antismash_6.1.1_meta %>% filter(S3320 <= 0.1) %>% group_by(product_1) %>% filter(n() > 8) %>% 
  ggplot(aes(fill=product_1, x=reorder(samples,S3320), y=count)) + 
  geom_bar(position="stack", stat="identity") + ylim(0,15) +
  scale_fill_brewer(palette = "Spectral") + theme_bw(base_size = 10) + 
  theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    #axis.title.x = element_text(margin = margin(t = -10)))
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "none"
  ) 
bar_bgcs_cosmo_product_1

bar_bgcs_tropical_product_1 <- df_regions_antismash_6.1.1_meta %>% filter(S3320 >= 0.1) %>% group_by(product_1) %>% filter(n() > 8) %>% 
  ggplot(aes(fill=product_1, x=reorder(samples,S3320), y=count)) + 
  geom_bar(position="stack", stat="identity") + ylim(0,15) +
  scale_fill_brewer(palette = "Spectral") + theme_bw(base_size = 10) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        #axis.title.x = element_text(margin = margin(t = -10)))
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        #legend.position = "none"
  ) 
bar_bgcs_tropical_product_1

bar_bgcs_cosmo_tropical <- plot_grid(bar_bgcs_cosmo_product_1, bar_bgcs_tropical_product_1, ncol = 2, align = "h")

ggsave(bar_bgcs_cosmo_tropical, file = "../../bar_bgcs_cosmo_tropical.png", width = 10, height =4, limitsize = FALSE)

###CAZY overall

bar_bgcs_cosmo_Cazyme_fam <-  ANI_mat_long_meta_cazy_antismash_seqfu %>% 
  select(samples, S3320, AA, CBM, CE, GH, GT, PL, cohesin, dockerin, SLH) %>%
pivot_longer(cols = -c(samples,S3320), names_to = "Cazyme_fam", values_to = "count") %>% filter(S3320 <= 0.1) %>% 
  ggplot(aes(fill=Cazyme_fam, x=reorder(samples,S3320), y=log(count))) + 
  geom_bar(position="stack", stat="identity") + #ylim(0,15) +
  scale_fill_brewer(palette = "Spectral") + theme_bw(base_size = 10) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        #axis.title.x = element_text(margin = margin(t = -10)))
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none"
  ) 


bar_bgcs_tropical_Cazyme_fam  <- ANI_mat_long_meta_cazy_antismash_seqfu %>% 
  select(samples, S3320, AA, CBM, CE, GH, GT, PL, cohesin, dockerin, SLH) %>%
  pivot_longer(cols = -c(samples,S3320), names_to = "Cazyme_fam", values_to = "count") %>% filter(S3320 >= 0.1) %>% 
  ggplot(aes(fill=Cazyme_fam, x=reorder(samples,S3320), y=log(count))) + 
  geom_bar(position="stack", stat="identity") + #ylim(0,15) +
  scale_fill_brewer(palette = "Spectral") + theme_bw(base_size = 10) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        #axis.title.x = element_text(margin = margin(t = -10)))
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        #legend.position = "none"
  ) 

bar_bgcs_cosmo_tropical_cazy <- plot_grid(bar_bgcs_cosmo_Cazyme_fam, bar_bgcs_tropical_Cazyme_fam, ncol = 2, align = "h")

ggsave(bar_bgcs_cosmo_tropical_cazy, file = "../../bar_bgcs_cosmo_tropical_cazy.png", width = 10, height =4, limitsize = FALSE)


###GHs

bar_bgcs_cosmo_GH <-  ANI_mat_long_meta_cazy_antismash_seqfu %>% 
  select(samples, S3320, starts_with("GH")) %>% select(-GH) %>% 
  pivot_longer(cols = -c(samples,S3320), names_to = "Cazyme_fam", values_to = "count") %>% filter(S3320 <= 0.1) %>% 
  group_by(Cazyme_fam) %>% filter(n() > 8) %>% 
  ggplot(aes(fill=Cazyme_fam, x=reorder(samples,S3320), y=log(count))) + 
  geom_bar(position="stack", stat="identity") + # ylim(0,15) +
  scale_fill_brewer(palette = "Spectral") + theme_bw(base_size = 10) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        #axis.title.x = element_text(margin = margin(t = -10)))
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none"
  ) 


bar_bgcs_tropical_Cazyme_fam  <- ANI_mat_long_meta_cazy_antismash_seqfu %>% 
  select(samples, S3320, AA, CBM, CE, GH, GT, PL, cohesin, dockerin, SLH) %>%
  pivot_longer(cols = -c(samples,S3320), names_to = "Cazyme_fam", values_to = "count") %>% filter(S3320 >= 0.1) %>% 
  ggplot(aes(fill=Cazyme_fam, x=reorder(samples,S3320), y=log(count))) + 
  geom_bar(position="stack", stat="identity") + #ylim(0,15) +
  scale_fill_brewer(palette = "Spectral") + theme_bw(base_size = 10) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        #axis.title.x = element_text(margin = margin(t = -10)))
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        #legend.position = "none"
  ) 

bar_bgcs_cosmo_tropical_cazy <- plot_grid(bar_bgcs_cosmo_Cazyme_fam, bar_bgcs_tropical_Cazyme_fam, ncol = 2, align = "h")
bar_bgcs_cosmo_tropical_cazy
ggsave(bar_bgcs_cosmo_tropical_cazy, file = "../../bar_bgcs_cosmo_tropical_cazy.png", width = 10, height =4, limitsize = FALSE)


bar_bgcs_cosmo_tropical_cazy_richness <- ANI_mat_long_meta_cazy_antismash_seqfu %>% 
  select(samples, S3320, AA, CBM, CE, GH, GT, PL, cohesin, dockerin, SLH) %>%
  pivot_longer(cols = -c(samples,S3320), names_to = "Cazyme_fam", values_to = "cazy_count") %>% 
  mutate(Group = ifelse(S3320 > 0.1, "Group 2 ~ Non-polar", "Group 1")) %>% na.omit()  %>%
ggplot(aes(x=Cazyme_fam, y=cazy_count, fill =Group)) + 
  geom_boxplot() + theme_bw(base_size = 10) + 
  scale_fill_brewer(palette = "Spectral") + theme_bw(base_size = 10) +
theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(face = "bold"),
      #axis.title.x = element_text(margin = margin(t = -10)))
      #axis.title.x = element_blank(),
      #axis.text.x = element_blank(),
      #legend.position = "none"
) 


ANI_mat_long_meta_cazy_antismash_seqfu %>% 
  select(samples, S3320, AA, CBM, CE, GH, GT, PL, cohesin, dockerin, SLH) %>%
  pivot_longer(cols = -c(samples,S3320), names_to = "Cazyme_fam", values_to = "cazy_count") %>% 
  mutate(Group = ifelse(S3320 > 0.1, "Group 2 ~ Non-polar", "Group 1")) %>% na.omit() %>% 
  group_by(Group) %>%
  distinct(samples) %>%
  summarise(n = n())

ggsave(bar_bgcs_cosmo_tropical_cazy_richness, file = "../../bar_bgcs_cosmo_tropical_cazy_richness.png", width = 10, height =4, limitsize = FALSE)


# df_regions_antismash_6.1.1_meta_id$group <- ifelse(df_regions_antismash_6.1.1_meta_id$S3320 < 0.1, "Group1", "Group2")
# 
# 
# bar_bgcs_id <- df_regions_antismash_6.1.1_meta_id %>% group_by(most_similar_known_cluster_id) %>% filter(n() > 8) %>% 
#   ggplot(aes(fill=most_similar_known_cluster_id, x=reorder(samples,S3320), y=count)) + 
#   geom_bar(position="stack", stat="identity") + ylim(0,5) +
#   scale_fill_brewer(palette = "Spectral") + theme_bw(base_size = 10) + 
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.title = element_text(face = "bold"),
#         axis.text = element_text(face = "bold"),
#         #axis.title.x = element_text(margin = margin(t = -10)))
#         axis.title.x = element_blank(),
#         #axis.text.x = element_blank(),
#         #legend.position = "none"
#   ) + facet_wrap(~group, scales="free_x")
# bar_bgcs_id






# df_regions_antismash_6.1.1_meta %>% group_by(product_1) %>% filter(n() > 8) %>%
# group_by(samples) %>% 
# ggplot(aes(x=S3320, fill=product_1)) +
#   geom_histogram(aes(y=..count../sum(..count..)), col="black", binwidth=0.005) + 
#   #geom_histogram(aes(y=..count../sum(..count..)),binwidth=0.01) + 
#   facet_wrap(~ product_1, ncol = 1, scales="free_y") + 
#   scale_fill_brewer(palette = "Spectral")

# df_regions_antismash_6.1.1_meta %>% group_by(product_1) %>% filter(n() > 4) %>% 
#   ggplot(aes(y=product_1, x=reorder(samples,S3320), fill=count)) +
#   geom_tile()

# df_regions_antismash_6.1.1_meta_id$Group <- ifelse(df_regions_antismash_6.1.1_meta_id$S3320 <= 0.1, "Group 1", "Group2")
# df_regions_antismash_6.1.1_meta_id %>% group_by(most_similar_known_cluster_description) %>% filter(most_similar_known_cluster_description!="",n() > 4) %>%
#   ggplot(aes(y=most_similar_known_cluster_description, x=as.factor(S3320), fill=count)) +
#   geom_tile() + facet_wrap(~Group,scales="free_x")
# # 
# # df_regions_antismash_6.1.1_meta %>% group_by(product_1) %>% filter(n() > 8) %>% 
# #   ggplot(aes(x=S3320, fill=product_1)) + 
# #   geom_density(..scaled..) + 
# #   facet_wrap(~ product_1, ncol = 1) + 
# #   scale_fill_brewer(palette = "Spectral")
# df_regions_antismash_6.1.1$most_similar_known_cluster_id
# df_regions_antismash_6.1.1 %>% group_by(most_similar_known_cluster_id) %>% dplyr::summarise(n()) %>% arrange(desc(`n()`))
#  #Core_bgcs





# Betalactone #####
betalactone_genomes <- df_regions_antismash_6.1.1_meta %>% filter(product_1 == "betalactone") %>%
  select(samples) %>% as.data.frame(.)

unique(betalactone_genomes$samples) #208 genomes contain betalactone BGC

df_regions_antismash_6.1.1 %>% filter(product_1 == "betalactone")
