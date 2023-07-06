library(tidyverse)
library(dplyr)
library(stringr)
library(ggpubr)
library(ggrepel)

ANIb_dist <- read.table(file = "../ANIb_output/ANIb_percentage_identity.tab", header = TRUE, row.names = 1,  sep = "\t")
str(ANIb_dist)

#clean up table
colnames(ANIb_dist) <-  gsub("\\_contigs.fixed", "", colnames(ANIb_dist))
rownames(ANIb_dist) <-  gsub("\\_contigs.fixed", "", rownames(ANIb_dist))

## remove photobacterium spp. "S3704" "S3856" "S3662" "S2753" "S2754" "S2541"
#remove <- c("S3704", "S3856", "S3662", "S2753", "S2754", "S2541")

#ANIb_dist <- ANIb_dist[-which(rownames(ANIb_dist) %in% remove),]
#ANIb_dist <- ANIb_dist[,-which(colnames(ANIb_dist) %in% remove)]





#convert to long format
ANIb_dist_long <- ANIb_dist %>%
  mutate(genome_id = rownames(.)) %>% 
  pivot_longer(cols = -genome_id, names_to = "quary", values_to = "ANI_similarity")

#### Add clades to the comparisons ######

#removing diagonal combinations
ANIb_dist_long <- ANIb_dist_long[ANIb_dist_long$quary!=ANIb_dist_long$genome_id,]

# Remove rows with duplicate combinations
ANIb_dist_long <- ANIb_dist_long[!duplicated(apply(ANIb_dist_long[, c("genome_id", "quary")], 1, function(x) paste(sort(x), collapse = ""))), ]



#Add metadata to quary and genome_id
ANIb_dist_long_meta <- df_Vibrio %>% select(genome_id, CTDTemperature, CTDSalinity, Taxonomy) %>% left_join(ANIb_dist_long, df_Vibrio, by="genome_id") 

colnames(ANIb_dist_long_meta)[2] <- "CTDTemperature_genome_id"
colnames(ANIb_dist_long_meta)[3] <- "CTDSalinity_genome_id"
colnames(ANIb_dist_long_meta)[4] <- "Taxonomy_genome_id"

colnames(ANIb_dist_long_meta)[1] <- "genome_id_1"
colnames(ANIb_dist_long_meta)[5] <- "genome_id"

ANIb_dist_long_meta <- df_Vibrio %>% select(genome_id, CTDTemperature, CTDSalinity, Taxonomy) %>% left_join(ANIb_dist_long_meta, df_Vibrio, by="genome_id") 

ANIb_dist_long_meta$temp_diff <- abs(c(ANIb_dist_long_meta$CTDTemperature - ANIb_dist_long_meta$CTDTemperature_genome_id))
ANIb_dist_long_meta$sal_diff <- abs(c(ANIb_dist_long_meta$CTDSalinity - ANIb_dist_long_meta$CTDSalinity_genome_id))
ANIb_dist_long_meta$shared.group <- ifelse(ANIb_dist_long_meta$Taxonomy_genome_id==ANIb_dist_long_meta$Taxonomy, 
                                           ANIb_dist_long_meta$Taxonomy_genome, "")


df_Vibrio %>% filter(S3320 < 0.001) %>% group_by(Taxonomy) %>% summarise(n())

ANIb_dist_long_meta %>% filter(shared.group=="echinoideorum") %>%
  ggplot(aes(temp_diff, 1-ANI_similarity, col=shared.group)) +
  #ggplot(aes(abs(Latitude), as.numeric(bgcs_count), col=as.numeric(cluster))) +
  #geom_point(aes(size=bgcs_count)) + scale_size(range = c(-5, 10)) 
  geom_point(alpha=0.3) + stat_smooth(method = "lm")


ANIb_dist_long_meta %>% filter(shared.group=="echinoideorum") %>%
  ggplot(aes(CTDTemperature,CTDTemperature_genome_id, col=shared.group)) +
  #ggplot(aes(abs(Latitude), as.numeric(bgcs_count), col=as.numeric(cluster))) +
  #geom_point(aes(size=bgcs_count)) + scale_size(range = c(-5, 10)) 
  geom_point(alpha=0.3) + stat_smooth(method = "lm")


ANIb_dist_long_meta %>% 
  ggplot(aes((temp_diff), 1-ANI_similarity)) +
  #ggplot(aes(abs(Latitude), as.numeric(bgcs_count), col=as.numeric(cluster))) +
  #geom_point(aes(size=bgcs_count)) + scale_size(range = c(-5, 10)) 
  geom_point(alpha=0.3) +  stat_smooth(method = "lm", col = "red")

summary(lm(ANIb_dist_long_meta$temp_diff ~ (1-ANIb_dist_long_meta$ANI_similarity)))$adj.r.squared 


#### 
#Create matrix from long format for distances


df <- ANIb_dist_long_meta %>% select(genome_id, genome_id_1, ANI_similarity) %>% as.data.frame(.)
str(df)
df <- df[-13469,]
# 
# 
# df[13469,]
# which(is.na(df$ANI_similarity))

library(igraph)

# Make undirected so that graph matrix will be symmetric
g <- graph.data.frame(df, directed=FALSE)

# add ANI_similarity as a weight attribute
distance_ANI_similarity <- get.adjacency(g, attr="ANI_similarity", sparse=FALSE)


df_2 <- ANIb_dist_long_meta %>% select(genome_id, genome_id_1, temp_diff) %>% as.data.frame(.)

df_2 <- df_2[-13469,]


library(igraph)

# Make undirected so that graph matrix will be symmetric
g <- graph.data.frame(df_2, directed=FALSE)

# add ANI_similarity as a weight attribute
distance_temp_diff <- get.adjacency(g, attr="temp_diff", sparse=FALSE)



ANIb_TempMCA <- mpmcorrelogram(as.dist(distance_ANI_similarity), as.dist(distance_temp_diff), method="spearman", permutations=999)


ANIb_test <- mpmcorrelogram(as.dist(distance_ANI_similarity), as.dist(distance_ANI_similarity), method="spearman", permutations=999, breaks=5)

summary(ANIb_test)



