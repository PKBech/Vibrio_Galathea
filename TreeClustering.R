library(tidyverse)
library(dplyr)
library(stringr)
library(ggpubr)
library(ggrepel)


#Import Phylogenomic distances
#Import autoMLST distance matrix and use S3320
#anvio_dist <- read.table(file = "anvio_tree160523.csv", header = FALSE, sep = ";")
anvio_dist <- read.table(file = "anvio_tree230523_noRef_noPhotobacterium.csv", header = FALSE, sep = ";")
rownames(anvio_dist) <- anvio_dist$V1
anvio_dist <- anvio_dist[,-1]
colnames(anvio_dist) <- rownames(anvio_dist)


PCoA=ape::pcoa(as.dist(anvio_dist))
str(PCoA)

PCOA_1_2 <- as.data.frame(PCoA$vectors[,1:2])

PCOA_1_2$genome_id <- rownames(PCOA_1_2)

PCOA_1_2_meta <- left_join(PCOA_1_2,df_Vibrio )

colnames(PCOA_1_2_meta)

PCOA_1_2_meta %>% 
  ggplot(aes(Axis.1, Axis.2, col=Taxonomy)) +
  #ggplot(aes(abs(Latitude), as.numeric(bgcs_count), col=as.numeric(cluster))) +
  #geom_point(aes(size = bgcs_count)) + ylim(-80,80) +
  geom_point(aes(size = bgcs_count)) 


PCOA_1_2_meta %>%   
  #filter(cluster!=0)  %>%
  ggplot(aes(Axis.1, Axis.2, col=as.factor(cluster))) +
  #ggplot(aes(abs(Latitude), as.numeric(bgcs_count), col=as.numeric(cluster))) +
  #geom_point(aes(size = bgcs_count)) + ylim(-80,80) +
  #geom_point(aes(shape = clade), size=4) 
  geom_point(aes(size = bgcs_count)) + scale_size(range = c(-5, 10)) 


##NMDS
library(vegan)

nMMDS=metaMDS(as.dist(anvio_dist))

scores(nMMDS, display = "sites")
nMDS_1_2 <- as.data.frame(scores(nMMDS, display = "sites"))

nMDS_1_2$genome_id <- rownames(nMDS_1_2)

nMDS_1_2_meta <- left_join(nMDS_1_2,df_Vibrio )

colnames(nMDS_1_2_meta)


nMDS_1_2_meta %>% 
  ggplot(aes(NMDS1, NMDS2, col=clade)) +
  #ggplot(aes(abs(Latitude), as.numeric(bgcs_count), col=as.numeric(cluster))) +
  #geom_point(aes(size = bgcs_count)) + ylim(-80,80) +
  geom_point(aes(shape = as.factor(clade), ), size=4) 
  #geom_point(aes(size = bgcs_count)) + scale_size(range = c(-5, 10)) 

nMDS_1_2_meta %>% 
  filter(cluster!=0)  %>%
  ggplot(aes(NMDS1, NMDS2, col=as.factor(cluster))) +
  #ggplot(aes(abs(Latitude), as.numeric(bgcs_count), col=as.numeric(cluster))) +
  geom_point(aes(size = bgcs_count))
  #geom_point(aes(shape = clade), size=4) 
#geom_point(aes(size = bgcs_count)) + scale_size(range = c(-5, 10)) 

###unsupervised clading dbscan
library(dbscan)

DB=dbscan::dbscan(as.dist(anvio_dist), eps = 0.001,minPts = 2)


DB$cluster

nMDS_1_2_meta$cluster <- DB$cluster


PCOA_1_2_meta$cluster <- DB$cluster

findEpsi=function(L, minRange=0,maxRange=3, steps=0.1, maxY=200, minP=10) {
  plot(0,0, col=0,xlim=c(minRange,maxRange), ylim=c(0,maxY))
  legend("topright", legend = c("nClust","nOutliers"), col=1:2, pch=16)
  for(i in seq(minRange,maxRange,steps)) {
    
    DB=dbscan::dbscan(L,i,minP)
    
    points(i,length(unique(DB$cluster)), pch=16, col=1)
    points(i,length(which(DB$cluster==0)), pch=16,col=2)
    
  }
}
findEpsi(as.dist(anvio_dist), minRange = 0,maxRange = .01,steps = 0.0001,maxY = 210,minP = 2)


nMDS_1_2_meta %>%  group_by(ClusterNumber,ribS2_taxon) %>% summarise(n()) %>% print(., n=38)

tree_cluster <- read.table(file = "../treecluster/phylogenomic-tree_noPhotobacterium_out", header = TRUE, sep = "\t")
colnames(tree_cluster)[1] <- "genome_id"

nMDS_1_2_meta<- left_join(nMDS_1_2_meta, tree_cluster)

nMDS_1_2_meta %>% 
  filter(cluster!=0) %>% 
  ggplot(aes(S3320, Latitude, col=as.factor(ribS2_taxon))) +
  #ggplot(aes(abs(Latitude), as.numeric(bgcs_count), col=as.numeric(cluster))) +
  #geom_point(aes(size = bgcs_count)) + ylim(-80,80) +
  geom_point(aes(size = bgcs_count)) + ylim(-80,80) +
  #geom_point(size=2) + ylim(-80,80) +
  #geom_smooth() +
  #geom_smooth(method = lm, se=FALSE) +
  #stat_cor(method = "spearman", p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name="rho", size = 2)+
  #scale_size(range = c(0.1, 10)) +
  #labs(x = "\nCore genes phylogenetic distance", y = "\nLatitude ") +
  theme_bw(base_size = 10) +
  #scale_colour_gradient(low = "blue", high = "red", na.value = NA) +
  # scale_colour_steps2(low = "white",
  #                     mid = "blue",
  #                     high = "red",
  #                     midpoint = 3) + 
  #scale_size(range = c(-5, 10)) +
  #scale_fill_manual(values = c("#f25a13", "#580a82"))+
  #scale_color_manual(values = c("#f25a13", "#580a82"))+
   scale_size(range = c(-5, 10))
  #annotate('rect', xmin=0, xmax=0.015, ymin=-Inf, ymax=Inf, alpha=.5, fill='red') +
 # annotate('rect', xmin=0.016, xmax=0.06, ymin=-Inf, ymax=Inf, alpha=.5, fill='blue') +
  #annotate('rect', xmin=0.13, xmax=0.17, ymin=-Inf, ymax=Inf, alpha=.3, fill='green') +
  #annotate('rect', xmin=0.18, xmax=0.20, ymin=-Inf, ymax=Inf, alpha=.3, fill='purple')
#+ facet_grid(.~cluster)
#scale_y_continuous(limits = c(0,1000))
