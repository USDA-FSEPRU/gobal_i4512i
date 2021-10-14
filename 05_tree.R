library(tidyverse)
library(ggtree)
library(ape)
library(ggplot2)
library(ggnewscale)
library(treeio)
library(tidytree)
library(caper)

# TODO #
# NEED A METADATA DF FOR TREE REPS ONLY, MASTERMETA WONT CUT IT #
# add in cluster size #
# CLASSIFY GENOMES BY TREE CLADES
tr <- read.raxml('./Jan_21_update/core_genome_tree2/RAxML_bipartitionsBranchLabels.core_genome_tree_2') 
str(tr)


# 'GCA_008862265'
# weirdos
# GCA_007978025
# GCA_009227915


master_meta <- read_tsv('output/04_metaJAN_UPDATE.tsv')



ggtree(tr)
# tr_filt <- drop.tip(tr, bad_tips)

ggtree(tr, layout = 'circular') 

# tr_filt %>% as_tibble() %>% arrange(desc(branch.length))


# hist(tr@phylo$edge.length, breaks = 1000)


# which(tr@phylo$edge.length > .002)
# 
# tr@phylo

# str(tree_tib)

# tree_tib <- tidytree::as_tibble(tr)

# tree_tib %>% arrange(desc(branch.length))

# offspring(.data = tree_tib, .node = 1891)

# tree_tib_filt <- tree_tib %>%
  # tidytree::filter(!parent %in% c(1890, 1891)) %>%
  # tidytree::filter(!(node %in% c(1890, 1891)))

# str(tree_tib_filt)

# as.phylo(tree_tib_filt)

# ggtree(as.phylo(tree_tib_filt))
# tr@phylo <- ape::drop.tip(phy=tr@phylo ,tip = c('SRR10232035', 'SRR8698045'))
# TR_RED <- drop.tip(phy=as.phylo(tr) ,tip = c('SRR10232035', 'SRR8698045'))
# tmp <- tr@data
# ggtree(tr, layout = 'circular',branch.length = 'none') #+ geom_tiplab(size=2)

# starting on bootstraps...


# 654 is the node that has good support to separate the south, but leaves out
# some genomes too...
#652 captures the south well but has poor bs support

# p$data$label
# clade.members(x = 1887, phy = tr@phylo, tip.labels = TRUE)
p <- ggtree(tr, layout = 'circular') +
  # geom_point2(aes(subset=!isTip & (bootstrap > 70),
  #                 fill=cut(bootstrap, c(0, 70, 90, 100))),
  #             shape=21, size=4) +
  # geom_text2(aes(subset=(!isTip & (bootstrap > 70)), label=node), size=4)+
  theme_tree(legend.position=c(.21, .91)) +
  geom_hilight(node=1387, fill="steelblue", alpha=0.5) +
  geom_hilight(node=926, fill='purple', alpha=.5) +
  geom_hilight(node=1022, fill='green', alpha=.5) 


ggtree(tr) +
  # geom_point2(aes(subset=!isTip & (bootstrap > 70),
  #                 fill=cut(bootstrap, c(0, 70, 90, 100))),
  #             shape=21, size=4) +
  geom_text2(aes(subset=(!isTip & (bootstrap > 50)), label=bootstrap), size=4)+
  geom_text2(aes(subset=(!isTip & (bootstrap > 50)), label=node),nudge_y = 15,color='blue', size=4)





ggtree(tr)+
  geom_hilight(node=1387, fill="steelblue", alpha=0.5) +
  geom_hilight(node=926, fill='purple', alpha=.5) +
  geom_hilight(node=1022, fill='green', alpha=.5) 


# view clade doesnt work with layout=circular

viewClade(ggtree(tr) , node=1387) 
viewClade(ggtree(tr) , node=926) 
# 1022 or 1024?!?!?!?
# 1022 captures the deeply branching members, but 1024 is probably better representative of a distinct clade
viewClade(ggtree(tr) , node=1022) 



clade_one <- caper::clade.members(1387, phy = tr@phylo, tip.labels = TRUE)
clade_two <- caper::clade.members(926, phy = tr@phylo, tip.labels = TRUE)
### ALERT HERE!!! 1022 vs 1024 (1024 was original) 
clade_three <- caper::clade.members(1024, phy = tr@phylo, tip.labels = TRUE)

clad_mat <- caper::clade.matrix(tr@phylo)

tree_tips <- tr@phylo$tip.label
length(tree_tips)


### GO BACK AND CHANGE THIS IN TREECHOOSE SCRIPT
# tree_meta <- 
#   read_tsv('./output/04_JAN_TREE_META.tsv') %>% 
#   filter(tree_rep %in% tree_tips)

tree_meta <- 
  read_tsv('./output/04_JAN_TREE_META2.tsv') %>% 
  filter(tree_rep %in% tree_tips)


# master_meta <- read_tsv('output/04_meta.tsv')
# weirdos <- master_meta %>% filter(ID %in% c('SRR10232035', 'SRR8698045'))


tree_meta <- 
  tree_meta %>%
  mutate(tree_clade=
           case_when(
             tree_rep %in% clade_one     ~ 'clade_one', 
             tree_rep %in% clade_two     ~ 'clade_two', 
             tree_rep %in% clade_three   ~ 'clade_three',
             TRUE                        ~ 'other'))


tree_meta$tree_clade %>% table()
LOOK <- tree_meta %>% filter(tree_clade == 'other')
# 
# hist(clad_mat$edge.length, breaks=1000)
# 
# 
# 
# test <- clad_mat$edge.length
# tibble(edge_num=names(clad_mat$edge.length), 
#        edge=paste(clad_mat$edge[,1], clad_mat$edge[,2], sep = '_'))

# 
# phy_clust <- read_tsv('./data/cluster_phydelity_k2_sol0_RAxML_bipartitions.txt') %>% 
#   mutate(rep_genome=TAXA, 
#          phy_cluster=CLUSTER) %>% dplyr::select(rep_genome, phy_cluster)
# 
# table(phy_clust$phy_cluster)
# 
# 
# tree_meta <- tree_meta %>% left_join(phy_clust)





###### FIND CLOSEST REFERENCE GENOME ####
phylo_dist <- ape::cophenetic.phylo(tr@phylo) %>%
  as.data.frame() %>% 
  rownames_to_column(var='from') %>% 
  pivot_longer(cols=-from,names_to = 'to', values_to='phydist')



these <- grep('GCA', phylo_dist$from, invert = TRUE)

phylo_refdist <- phylo_dist[these,]
unique(phylo_refdist$from)


phylo_refdist <- 
  phylo_refdist %>% 
  group_by(to) %>% 
  summarise(closest_ref=from[which.min(phydist)]) %>% 
  mutate(tree_rep=to) %>%
  dplyr::select(-to)

tree_meta <-
  tree_meta %>%
  left_join(phylo_refdist)



### WRITE OUT MASTERMETA HERE

clust_helper <- tree_meta %>%
  dplyr::select(tree_cluster, tree_clade, closest_ref)

master_meta <- master_meta %>% left_join(clust_helper)

write_tsv(master_meta, './output/05_JAN_meta.tsv')

### overview tree construction #

TEST_MET <- tree_meta %>%
  filter(tree_rep %in% tr@phylo$tip.label) %>% 
  dplyr::select(tree_rep, everything())


p <- ggtree(tr, layout = 'circular') %<+% TEST_MET 


p0 <- p + 
  geom_hilight(node=1387, fill="steelblue", alpha=0.5) +
  geom_hilight(node=926, fill='purple', alpha=.5) +
  geom_hilight(node=1024, fill='green', alpha=.5) 

library(stringr)

ISO_df <- TEST_MET %>% 
  transmute(tree_rep=tree_rep, 
            Source=str_to_title(isolation), 
            Source=case_when(
              Source == 'Cow'   ~ 'Bovine', 
              TRUE              ~  Source
            )) %>% 
  column_to_rownames(var='tree_rep')
  

ISLAND_DF <- TEST_MET %>%
  transmute(tree_rep=tree_rep,
            SGI4=perc_full_SGI4,
            MDR=perc_full_MDR,
            sopE=perc_full_sopE_mtmV) %>%
  column_to_rownames(var='tree_rep')


#TODO# move this to cluster script
leid_DF <- TEST_MET %>% 
  dplyr::select(tree_rep, leiden_cluster) %>% 
  transmute(tree_rep=tree_rep,
            Accessory=leiden_cluster) %>% 
  column_to_rownames(var='tree_rep')

# year_df <- TEST_MET %>%
#   dplyr::select(tree_rep, year) %>% 
#   mutate(year=as.numeric(year)) %>% 
#   column_to_rownames(var='tree_rep')

cont_df <- TEST_MET %>%
  transmute(tree_rep=tree_rep,
            Continent=continent) %>% 
  column_to_rownames(var='tree_rep')

looks <- p$data
p0$data$label <- sub('_92-0392','',p0$data$label)
# p$data$x
# p$data$angle[p$data$isTip]
p1 <- p0 +
  # geom_tippoint(aes(color=continent), alpha=.75)+
  # geom_tiplab(align = TRUE, size=2)+
  geom_tiplab(aes(subset=label %in% c('LT2',
                                     'USDA15WA1',
                                     'SL1344',
                                     'Enteritidis', 
                                     'L-4234', 
                                     'SO4698-09', 
                                     'TW-Stm6'),
                  angle=angle, 
                  x=.00021,
                  geom='label'),
              # angle=c(-90, 90, 90, 90, 90,90,90),
              size=2.75)#+
  # geom_point2(aes(subset=(!isTip & (bootstrap > 90))),color='black', size=1)
  # theme(legend.position = "right") +
  # scale_color_brewer(palette = 'Dark2')

p1

# p1 + geom_tippoint(aes(color=phy_cluster), alpha=.75)


p2 <- gheatmap(p1, leid_DF, offset=.0000095, width=.1,
         colnames_angle=-90, 
         colnames_offset_y = .25,
         font.size = 2, 
         legend_title = 'Accessory Genome') +
  scale_fill_manual(name='Accessory Genome', values=c(A='skyblue',B='forestgreen',C='orange'), na.translate=F) +
  new_scale_fill()
  # scale_color_brewer(palette = 'Set1')
p2

# p3 <- p2 + new_scale_fill()

# this is for MDR SGI4 ETC
p3 <- gheatmap(p2, ISLAND_DF,
               offset = .000039,
               width = .2,
               colnames_angle = -90,
               color = 'grey',
               low = 'white', high = 'red', 
               font.size = 2, 
               legend_title = 'Percent Present') +
  new_scale_fill()

p3

# ISOLATION sources

unique(ISO_df$Source)

library(RColorBrewer)

SOURCE_COLORS <- brewer.pal(length(unique(ISO_df$Source)), name = 'Set1')
names(SOURCE_COLORS) <- levels(factor(ISO_df$Source))


p4 <- gheatmap(p3, ISO_df, 
               offset = .000095,
               width = .1,
               colnames_angle = -90,
               color = 'grey', font.size = 2, 
               legend_title = 'Source') +
  scale_fill_manual(name='Source', values = SOURCE_COLORS, na.translate=F) +
  new_scale_fill()


p4


# continent
p5 <- gheatmap(p4, cont_df,
            offset= .000129,
            width=.1,
            colnames_angle=-90,
            color='grey', font.size=2) + 
  scale_fill_brewer(name='Continent', palette = 'Dark2', na.translate=F)


p5


ggsave(filename = './output/OVERVIEW_TREE_JAN_NEW.jpeg', width = 10, height = 10, units = 'in')

