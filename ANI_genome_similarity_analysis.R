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

## Isolate the big clade:
Big_clade <- df_Vibrio %>% filter(S3320 < 0.04) %>% select(genome_id) %>% as.data.frame(.)

ANIb_dist_cosmopolitans <- ANIb_dist[Big_clade$genome_id, Big_clade$genome_id]

##NMDS
library(vegan)


#NMDSQ <- metaMDS(unifracQ, distance = "bray", k = 3, trymax = 500, autotransform = T)
#nMMDS=metaMDS(as.dist(ANIb_dist),  trymax = 50 )
 
nNMDS=metaMDS(ANIb_dist, distance = "euclidian")
nNMDS=metaMDS((ANIb_dist_cosmopolitans ), distance = "euclidian", k=3)



nMDS_1_2 <- as.data.frame(scores(nNMDS, display = "sites"))

nMDS_1_2$genome_id <- rownames(nMDS_1_2)

###tree_cluster cluster def
tree_cluster <- read.table(file = "../treecluster/phylogenomic-tree_noPhotobacterium_out", header = TRUE, sep = "\t")

colnames(tree_cluster)[1] <- "genome_id"

nMDS_1_2_meta <- left_join(nMDS_1_2,df_Vibrio )
nMDS_1_2_meta<- left_join(nMDS_1_2_meta, tree_cluster)


colnames(nMDS_1_2_meta)

nMDS_1_2_meta %>% filter(genome_id=="S4053")

nMDS_1_2_meta$ClusterNumber <- as.factor(nMDS_1_2_meta$ClusterNumber)

nMDS_1_2_meta %>% filter(S3320 < 0.04) %>%#filter(ClusterNumber %in% c("2", "1",  "5", "6","8")) %>% 
  ggplot(aes(NMDS1, NMDS3, col=CTDTemperature)) +
  #ggplot(aes(abs(Latitude), as.numeric(bgcs_count), col=as.numeric(cluster))) +
  geom_point(aes(size = bgcs_count)) #+ scale_size(range = c(-5, 10)) 
  #geom_point(size=4) 
#geom_point(aes(size = bgcs_count)) + scale_size(range = c(-5, 10)) 
nMDS_1_2_meta %>% filter(NMDS1<0.01, NMDS2<0.01)

### PERMANOVA ####

ANIb_dist_meta <- ANIb_dist %>% mutate(genome_id = rownames(.)) %>%
left_join(., df_Vibrio, by="genome_id")  %>% 
  left_join(., tree_cluster) %>%
  select(genome_id, Latitude, ClusterNumber, Longitude, Type)

ANIb_dist_meta$Latitude <- round(abs(ANIb_dist_meta$Latitude), digits = 0)
ANIb_dist_meta$Latitude <- as.factor(ANIb_dist_meta$Latitude)
ANIb_dist_meta$Longitude <- round(ANIb_dist_meta$Longitude, digits = 0)
ANIb_dist_meta$Longitude <-  as.factor(ANIb_dist_meta$Longitude)
ANIb_dist_meta$Type <-  as.factor(ANIb_dist_meta$Type)


ANIb_dist_meta$ClusterNumber <- as.factor(ANIb_dist_meta$ClusterNumber)

str(ANIb_dist_meta)

### permanova to test for site effects
dist.mat <- vegdist(ANIb_dist, method = "bray")

perm.log <- adonis2(dist.mat ~ Type, data=ANIb_dist_meta, permutations = 999, method = "bray")
perm.log

#Cut-off treecluster 0.4 
#Df SumOfSqs      R2        F Pr(>F)    
#ClusterNumber           25 0.185621 0.96300 6808.709  0.001 ***
#  Latitude                67 0.005024 0.02607   68.769  0.001 ***
#  ClusterNumber:Latitude  33 0.002023 0.01049   56.212  0.001 ***
#  Residual                78 0.000085 0.00044                    
#Total                  203 0.192754 1.00000                    
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#####
dist.mat <- vegdist(ANIb_dist, method = "bray")
perm.log <- adonis2(dist.mat ~ Latitude, data=ANIb_dist_meta, permutations = 999, method = "bray")
perm.log

perm.log <- adonis2(dist.mat ~ Latitude, data=ANIb_dist_meta, permutations = 999, method = "bray")
perm.log

#######



#convert to long format
ANIb_dist_long <- ANIb_dist %>%
  mutate(genome_id = rownames(.)) %>% 
  pivot_longer(cols = -genome_id, names_to = "quary", values_to = "ANI_similarity")


#ANIb_dist_long$comb <- paste(ANIb_dist_long$quary, ANIb_dist_long$genome_id )

#ANIb_dist_long_meta <- left_join(ANIb_dist_long, df_Vibrio) %>% select(quary,genome_id,comb,ANI_similarity,Latitude, Longitude, bgcs_count)

#removing diagonal combinations
ANIb_dist_long <- ANIb_dist_long[ANIb_dist_long$quary!=ANIb_dist_long$genome_id,]

# Remove rows with duplicate combinations
ANIb_dist_long <- ANIb_dist_long[!duplicated(apply(ANIb_dist_long[, c("genome_id", "quary")], 1, function(x) paste(sort(x), collapse = ""))), ]







#Add metadata to quary and genome_id
ANIb_dist_long_meta <- df_Vibrio %>% select(genome_id, Latitude, Longitude) %>% left_join(ANIb_dist_long, df_Vibrio, by="genome_id") 
colnames(ANIb_dist_long_meta) <- c("genome_id_1", "Latitude_genome_id","Longitude_genome_id", "genome_id", "ANI_similarity")
ANIb_dist_long_meta <- df_Vibrio %>% select(genome_id, Latitude, Longitude) %>% left_join(ANIb_dist_long_meta, df_Vibrio, by="genome_id") 
colnames(ANIb_dist_long_meta)[2] <- "Latitude_genome_q"
colnames(ANIb_dist_long_meta)[3] <- "Longitude_genome_q"


#Calculate geospatial distance between two points (lat,long)

ANIb_dist_long_meta$Distance <- NA  # Create an empty column for distances

library(geosphere)

for (i in 1:nrow(ANIb_dist_long_meta)) {
  coords_q <- c(ANIb_dist_long_meta$Longitude_genome_q[i], ANIb_dist_long_meta$Latitude_genome_q[i])
  coords_id <- c(ANIb_dist_long_meta$Longitude_genome_id[i], ANIb_dist_long_meta$Latitude_genome_id[i])
  
  distance <- distGeo(coords_q, coords_id)
  
  ANIb_dist_long_meta$Distance[i] <- distance  # Assign the calculated distance to the new column
}



ANIb_dist_long_meta %>% 
  ggplot(aes(Distance, ANI_similarity)) +
  #ggplot(aes(abs(Latitude), as.numeric(bgcs_count), col=as.numeric(cluster))) +
  #geom_point(aes(size=bgcs_count)) + scale_size(range = c(-5, 10)) 
  geom_point(alpha=0.3) +  stat_smooth(method = "lm", col = "red")

summary(lm(ANIb_dist_long_meta$Distance ~ ANIb_dist_long_meta$ANI_similarity))$adj.r.squared 

library(mpmcorrelogram)


df_Vibrio_sort <- df_Vibrio %>% filter(genome_id %in% rownames(ANIb_dist)) %>% arrange(match(genome_id, rownames(ANIb_dist)))


Geodist <- distm(cbind(df_Vibrio_sort$Longitude, df_Vibrio_sort$Latitude), fun=distGeo)
rownames(Geodist) <- df_Vibrio_sort$genome_id
colnames(Geodist) <- df_Vibrio_sort$genome_id
# Geodist <- as.dist(Geodist)

rownames(ANIb_dist) == rownames(Geodist) 

ANIb_GeoMCA <- mpmcorrelogram(as.dist(ANIb_dist), as.dist(Geodist), method="spearman", permutations=999)

summary(ANIb_GeoMCA)




