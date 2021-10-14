library(tidyverse)

# leid <- read_csv('output/02_leiden_clust_JAN_UPDATE.csv') %>%
#   transmute(genome=genome, leid_clust=leiden_cluster) %>% 
#   filter(leid_clust != 3)
# 
# 
# leid %>% count(leid_clust) 
# # louv <- read_csv('output/02_louvain_clust.csv') %>%
#   # transmute(genome=genome, louv_clust=louvain_cluster)
# # unique(master_meta$continent)
# 



# PROBLEM, WITH TREE REPS. though the tree only took 48 hours.... can prob run over the weekend.
CURRENT_TREE_REPS <- read_lines("./Jan_21_update/TREE_REPS2.txt")


master_meta <- 
  read_tsv('output/03_meta_JAN.tsv',guess_max = 10000) 


master_meta %>% 
  group_by(leiden_cluster) %>% 
  tally()

# these are suspicious tips identified from the first iteration
bad_tips <-c('GCA_007978025','GCA_009227915', "GCA_008862265", "GCA_007603015", 'GCA_007029205')



master_meta <- 
  master_meta %>% 
  filter(leiden_cluster != 'none') %>% 
  filter(!(genome %in% bad_tips)) %>% 
  group_by(continent, year, ISOLATION, PDS_acc, leiden_cluster) %>% 
  mutate(tree_cluster=cur_group_id()) %>% 
  ungroup()

# tree_groups <- 
#   master_meta %>% 
#   summarise(tree_cluster=tree_cluster[1], 
#             num_genomes=length(genome), 
#             rep_genome=genome[1], 
#             MDR_types=paste(unique(MDR_class),sep = '_', collapse = '_'),
#             SGI4_types=paste(unique(SGI4_class),sep = '_', collapse = '_'),
#             av_MDR_cov=mean(MDR_length), 
#             av_SGI4_cov=mean(SGI4_length), 
#             av_MDR_perc=(av_MDR_cov/28186)*100, 
#             av_SGI4_perc=(av_SGI4_cov/80687)*100, 
#             n_full_SGI4=sum(SGI4_class == 'full'), 
#             n_full_MDR =sum(MDR_class == 'full'), 
#             n_full_SGI4=sum(SGI4_class == 'full')/length(SGI4_class) * 100, 
#             n_full_MDR =sum(MDR_class == 'full')/length(MDR_class) * 100) %>%
#   arrange(desc(num_genomes)) %>%
#   ungroup() %>% 
#   write_tsv('./output/tree_dereplication2.tsv')

# tree_groups2 <- read_tsv('./output/tree_dereplication.tsv')
# writing 2nd just to be sure not overwriting somethign imp
# tree_groups %>% write_tsv('./output/tree_dereplication.tsv')


tree_groups <- 
  master_meta %>%
  group_by(tree_cluster) %>% 
  summarise(tree_cluster=tree_cluster[1], 
            num_genomes=length(genome), 
            rep_genome=genome[1])





tree_reps <- 
  tree_groups %>%
  pull(rep_genome) #%>%
  # write_lines("./Jan_21_update/TREE_REPS2.txt")
all(CURRENT_TREE_REPS %in% tree_reps)

# apparently 3 of the genomes in TREE_REPS2 have allready been filtered from the dataset?
# oh well.  We're just going to an overview, not going to sweat it too much

# master_meta %>%
#   group_by(`SNP cluster`, leid_clust) %>%
#   tally() %>%
#   arrange(desc(n))
master_meta %>%
  group_by(leiden_cluster) %>%
  tally() %>%
  arrange(desc(n))

# master_meta %>% filter(is.na(leid_clust))
# master_meta[grep('t',master_meta$ID),]

hist(log(tree_groups$num_genomes), breaks=100)

master_meta <- 
  master_meta %>%
  left_join(tree_groups) %>% 
  mutate(tree_rep = genome %in% tree_reps) %>% 
  write_tsv('./output/04_metaJAN_UPDATE.tsv')



TREE_META <- 
  master_meta %>%
  group_by(tree_cluster) %>% 
  summarise(tree_rep=genome[tree_rep],
            all_genomes=paste(genome, sep = ';', collapse = ';'),
            continent = unique(continent), 
            perc_full_SGI4 = sum(targ_presence_SGI4 == 'full') / length(targ_presence_SGI4) * 100, 
            perc_full_MDR  = sum(targ_presence_MDR == 'full') / length(targ_presence_MDR) * 100,
            perc_full_thrW18kb= sum(targ_presence_thrW18kb_GI == 'full') / length(targ_presence_thrW18kb_GI) * 100,
            perc_full_sopE_mtmV=sum(targ_presence_SO4698_sopE_mTmV == 'full') / length(targ_presence_SO4698_sopE_mTmV) * 100,
            perc_full_USDA15WA1_PRO = sum(targ_presence_USDA15WA1_581ProFrag == 'full') / length(targ_presence_USDA15WA1_581ProFrag) * 100,
            perc_full_NOsopE_PRO = sum(targ_presence_SO4698_NOsopE_PRO == 'full') / length(targ_presence_SO4698_NOsopE_PRO) * 100,
            year = unique(year), 
            num_genome = length(genome), 
            isolation = unique(ISOLATION), 
            leiden_cluster = unique(leiden_cluster))

## CHANGED THIS TO 2 BECAUSE I THINK SOME THINGS ARE FUCKY
write_tsv(TREE_META, './output/04_JAN_TREE_META2.tsv')





master_meta$targ_presence_SO4698_NOsopE_PRO
master_meta$targ_presence_USDA15WA1_581ProFrag
###
"A pangenome containing all 5365 i4512i- genomes passing QC as well as the reference genomes,
LT2, Enteriditis, and SL1344 was constructed with ppanggolin.  ppangolin grouped the genes in the 
pan genome into 3 partitions, persistent (X genes), shell (X genes) and cloud (X genes).
To identify genomes that share similar accessory genes (cloud partition) a graph was contstructed
in which the nodes are genomes and the edges join nodes according to the number of cloud genes 
they share, that is, the edges of this graph are weighted according to the number of shared genes.
Once this graph was constructed we applied the leiden community detection algorithm.  This 
This algorithm assigned genomes into 3 clusters, genomes within the same cluster tend to share more cloud 
genes.  

To assess the phylogenetic relatedness of the isolates we sought to construct a phylogenetic tree based of 
an alignment of the core genome.  To ease the computational burden as well as assist in visualization
genomes were clustered and representative genomes chosen according to the following scheme:
Genomes were grouped by 1) country of origin, 2) year of isolation, 3) `Isolation type` (environmental or clinical), 
4) `SNP cluster` (as assigned by NCBI pathogen detection pipeline), 5) accessory genome cluster (from ppanggolin pangenome), and 6) Associated host.  This produced
produced 652 groups and a representative genome was chosen at random from within each.  The core genes from 
these 652 genomes plus the 4 reference strains were aligned using the pipeline in Roary, mafft
A maximum likelyhood phylogenetic tree was produced from this alignment of the core genes with raxml and 
visualized with ggtree."



### Now a section for ag animal isolates only ###
# this s getting hairy... need to reduce data in a way that doesnt bias the transmission results

#split into 
# North and South America
# Europe Asia Oceania 

# remove PDS_acc  with low # total isolates?  Not sure about this one...
#' 
#' master_meta %>%
#'   filter(global_area == 'Americas') %>%
#'   count(PDS_acc) %>% 
#'   filter(n > 3)
#' 
#' 
#' master_meta %>% group_by(continent) %>% tally()
#' 
#' master_meta <- 
#'   master_meta %>% 
#'   mutate(global_area=case_when(
#'     continent %in% c('North America', 'South America')      ~ 'Americas',
#'     continent %in% c('Europe', 'Asia', 'Oceania', 'Africa') ~ 'Eurasia',
#'     TRUE                                                    ~ 'ERROR'))
#' 
#' 
#' master_meta <- master_meta %>% ungroup()
#' Americas_non_human <- 
#'   master_meta %>% 
#'   filter(ISOLATION != 'human') %>%
#'   filter(global_area == 'Americas')
#' 
#' Americas_human <- 
#'   master_meta %>%
#'   filter(ISOLATION == 'human' & global_area == 'Americas') %>% 
#'   filter(!(country == 'USA' & state == 'unknown'))
#' 
#' 
#' Americas_non_human %>% count(country)
#' Americas_non_human %>% count(ISOLATION)
#' 
#' Americas_human %>% 
#'   group_by(year, PDS_acc, country, state, leid_clust) %>% 
#'   tally() %>%
#'   arrange(desc(n))
#' 
#' Americas_non_human %>% 
#'   group_by(year, PDS_acc, country, state, ISOLATION, leid_clust) %>% 
#'   tally() %>% 
#'   arrange(desc(n))
#' 
#' 
#' Americas_non_human <- 
#'   Americas_non_human %>%
#'   filter(!(country == 'USA' & state == 'unknown' & ISOLATION == 'other'))
#' 
#' 
#' Americas_human %>% 
#' 
#' Americas_non_human %>% count(ISOLATION)
#' 
#' Americas_non_human %>% 
#'   group_by(year, PDS_acc, country, state, ISOLATION, leid_clust) %>% 
#'   tally() %>% 
#'   arrange(desc(n))
#' 
#' 
#' 
#' 
#' 
#' Eurasia_non_human <-
#'   master_meta %>%
#'   filter(ISOLATION != 'human') %>% 
#'   filter(global_area == 'Eurasia')
#' 
#' Eurasia_non_human %>% count(ISOLATION)
#' 
#' 
#' unique(master_meta$ISOLATION)
#' 
#' #'Felis catus' = 'Cat'
#' #'Canis' = 'Dog'
#' #'Ovis aries' = 'goat'
#' 
#' # 
#' # master_meta <- master_meta %>% 
#' #   mutate(iso_source = fct_recode(iso_source,
#' #                                  Cat='Felis catus',
#' #                                  Dog='Canis', 
#' #                                  goat='Ovis aries'))
#' 
#' 
#' hums <- master_meta %>% filter(ISOLATION == 'human')
#' 
#' 
#' hums <- 
#'   hums %>%
#'   group_by(country,state, year, ISOLATION, PDS_acc, leid_clust) %>% 
#'   mutate(hum_cluster=cur_group_id())
#' 
#' 
#' humsum <- 
#'   hums %>% 
#'   summarise(tree_cluster=tree_cluster[1], 
#'           num_genomes=length(genome), 
#'           rep_genome=genome[1], 
#'           MDR_types=paste(unique(MDR_class),sep = '_', collapse = '_'),
#'           SGI4_types=paste(unique(SGI4_class),sep = '_', collapse = '_'),
#'           av_MDR_cov=mean(MDR_length), 
#'           av_SGI4_cov=mean(SGI4_length), 
#'           av_MDR_perc=(av_MDR_cov/28186)*100, 
#'           av_SGI4_perc=(av_SGI4_cov/80687)*100, 
#'           n_full_SGI4=sum(SGI4_class == 'full'), 
#'           n_full_MDR =sum(MDR_class == 'full'), 
#'           n_full_SGI4=sum(SGI4_class == 'full')/length(SGI4_class) * 100, 
#'           n_full_MDR =sum(MDR_class == 'full')/length(MDR_class) * 100) %>%
#'   arrange(desc(num_genomes)) %>%
#'   ungroup()
#' 
#' hums4tree <- humsum %>% pull(rep_genome)
#' nonhum4tree <- non_human %>% pull(ID)
#' 
#' transmission_tree <- c(hums4tree, nonhum4tree)
#' 
#' 
#' transmission_tree %>% write_lines('HUM_NONHUM_TRANSMISSION_TREE.txt')
