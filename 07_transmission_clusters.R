library(tidyverse)
library(ape)
library(cowplot)

theme_set(theme_cowplot(font_family = 'arial'))
theme_update(panel.grid.major=element_line(color='grey80'), 
             panel.border = element_rect(color='black'),
             axis.title=element_text(face = 'bold'))



mm <- read_tsv('output/06_master_meta.tsv')

good_PDS <- mm %>% count(PDS_acc) %>% filter(n > 5) %>% arrange(desc(n))

good_PDS

mm$collected_by %>% grepl(., 'FSIS')

mm %>% group_by(continent, country) %>% tally() %>% arrange(desc(n))
mm %>% group_by(continent, country, tree_clade) %>% tally() %>% write_tsv('CONTINENT_COUNTRY_TREE_CLADE.tsv')
mm %>% group_by(continent, country, tree_clade, closest_ref) %>% tally() %>% write_tsv('CONTINENT_COUNTRY_TREE-CLADE_CLOSEST-REF.tsv')

# tr <- read.tree('./Jan_21_update/SNP_trees/PDS000004605.21.newick_tree.newick')
# 
# 
# write.tree(tr, file = 'TEST_TREE.newick')
# 
# tree_path <- './Jan_21_update/SNP_trees/PDS000004605.21.newick_tree.newick'
# 
# sub('./Jan_21_update/SNP_trees/(PDS[0-9]+.[0-9]+).newick_tree.newick','\\1',tp)


# For some reason phydelity doesnt like the trees directly from ncbi, 
# I have to read them in with ape and then write back out and they are fine...
if (!(dir.exists('FIXED_SNP_TREES'))) {
  
  system('mkdir FIXED_SNP_TREES')
  
  fix_trees <- 
    function(tree_path){
      tr <- read.tree(tree_path)
      SNP_name <- sub('./Jan_21_update/SNP_trees/(PDS[0-9]+.[0-9]+).newick_tree.newick','\\1',tree_path)
      new_path <- paste0('./FIXED_SNP_TREES/', SNP_name, '.newick')
      write.tree(tr, new_path)
    }
  
  
  tree_files <- list.files('./Jan_21_update/SNP_trees', pattern = 'newick', full.names = T)
  
  fix_trees(tree_files[1])
  
  lapply(tree_files, fix_trees)
  
  
  
  
}


trans_clust_files <- 
  list.files(path = './FIXED_SNP_TREES/', pattern = 'txt', full.names = T)


read_phydelity <- 
  function(phy_txt_path){
    SNP_clust <- 
      sub('./FIXED_SNP_TREES//cluster_phydelity_k.*_(PDS[0-9]+.[0-9]+).txt',
          '\\1',
          phy_txt_path)
    read_tsv(phy_txt_path) %>% 
      transmute(target_acc=TAXA, 
                PDS_acc = SNP_clust, 
                phy_clust=CLUSTER)
    
  
}


phy_res <- lapply(trans_clust_files, read_phydelity)

Trans_clusts <- 
  bind_rows(phy_res) %>%
  group_by(phy_clust, PDS_acc) %>% 
  mutate(PHY_CLUST = cur_group_id()) %>% 
  ungroup() %>% 
  select(target_acc, PDS_acc, PHY_CLUST) %>% 
  arrange(PHY_CLUST)


# clean this one up?
master_meta <- 
  Trans_clusts %>%
  right_join(mm) %>% 
  write_csv('./output/07_meta.csv')

### Trans cluster summary ###
# num trans clusters
# num genomes in trans clusters
# biggest trans cluster

# FOR THESE INTERESTING CLUSTERS I WANT TO KNOW WHAT % OF THE SNP CLUSTER IS 
# REPRESENTED IN THAT PHY CLUST

master_meta %>%
  group_by(PDS_acc) %>% 
  summarise(num_PHY=length(unique(PHY_CLUST)))


PDS_tots <- master_meta %>% count(PDS_acc, name='PDS_total')

PHY_tots <- master_meta %>% count(PHY_CLUST)

singleton_transmission_clusters <- 
  PHY_tots %>% 
  arrange((n)) %>% 
  filter(n==1) %>% 
  pull(PHY_CLUST) %>%
  as.character()

master_meta %>% 
  group_by(PDS_acc) %>% 
  summarise(num_PHY = length(unique(PHY_CLUST))) %>% 
  arrange(desc(num_PHY))



master_meta %>%
  group_by(PHY_CLUST, leiden_cluster) %>%
  tally() %>% ungroup() %>% 
  group_by(PHY_CLUST) %>% 
  summarise(num_leid=length(leiden_cluster)) %>% filter(num_leid > 1) %>% 
  arrange(desc(num_leid))
#



### multihost trans clust ###

transmission_summary <- 
  master_meta %>% 
  filter(!(PHY_CLUST %in% singleton_transmission_clusters)) %>%
  filter(!is.na(PHY_CLUST)) %>% 
  group_by(PHY_CLUST) %>%
  summarise(num_isolates=n(),
            PDS_acc=unique(PDS_acc), 
            tree_clade = paste(unique(tree_clade), collapse = '_'),
            countries=paste(unique(country), collapse = '_'), 
            n_countries=length(sort(unique(country))),
            continents=paste(unique(continent), collapse = '_'),
            ISOLATIONS=paste(sort(unique(ISOLATION)), collapse = '_'), 
            num_ISO=length(unique(ISOLATION)),
            genomes=paste(unique(genome), collapse = '_'), 
            acc_genome = paste(leiden_cluster, collapse = '_'),
            SGI4 = paste(targ_presence_SGI4, collapse = '_'), 
            MDR  = paste(targ_presence_MDR , collapse = '_'), 
            sopE_mTmV = paste(targ_presence_SO4698_sopE_mTmV, collapse = '_'), 
            thrW18kb_GI = paste(targ_presence_thrW18kb_GI, collapse = '_'), 
            raw_iso = paste(isolation_source, collapse = '_')
            
            ) %>% 
  arrange(desc(num_ISO)) %>% 
  left_join(PDS_tots) %>% 
  mutate(percent_PDS=num_isolates/PDS_total) %>% 
  write_csv('./data/transmission_summary.csv')


library(cowplot)

cont_order <- c('North America', 'Europe','South America', 'Asia',
                'North America and Europe', 'North America and South America',
                'Europe and North America', 'Asia and North America', 'Asia and Europe')


transmission_summary <- 
  read_csv('./data/transmission_summary.csv') %>% 
  mutate(ISOLATIONS=gsub('_',' / ',ISOLATIONS), 
         continents = gsub('_', ' and ', continents), 
         continents = factor(continents, levels=cont_order), 
         tree_clade=factor(tree_clade, levels = c('clade_one', 'clade_two', 'clade_three')), 
         `Tree Clade`=fct_recode(tree_clade, `Clade 1`='clade_one', `Clade 2`='clade_two', `Clade 3`='clade_three')) 

transmission_summary$`Tree Clade`

# human only transmission clusters
transmission_summary %>%
  filter(ISOLATIONS == 'Human') %>% 
  # mutate(ISOLATIONS=fct_rev(fct_infreq(f = ISOLATIONS,ordered = T ))) %>%
  # filter(ISOLATIONS %in% c('human', 'human_other', 'other', 'human_swine', 'swine', 'chicken_human', 'cow_human', 'chicken', 'human_other_swine', 'other_swine')) %>% 
  ggplot(aes(x=tree_clade, fill=continents, )) + 
  geom_bar(width = .75, color='black') +
  coord_flip() + 
  ylab('number of transmission clusters') + 
  xlab('Host combinations') + 
  ggtitle('Number of transmission clusters only associated with humans') + 
  theme_cowplot() 



"human_swine	103
swine	61
chicken_human	36

cow_human	27
human_other_swine	24
chicken	20
other_swine	18
cow	12
human_turkey	15
turkey	10
cow_swine	10"

# top 9 transmission host-transmission patterns
# Brad asked for this one
transmission_summary %>%
  filter(percent_PDS !=1) %>% 
  mutate(ISOLATIONS=fct_rev(fct_infreq(f = ISOLATIONS,ordered = T )), 
         tree_clade=factor(tree_clade, levels = c('clade_one', 'clade_two', 'clade_three'))) %>%
  filter(ISOLATIONS %in% c('Turkey', 'Bovine / Swine', 'Human / Turkey', 'Human / Swine', 'Swine', 'Chicken / Human', 'Bovine / Human', 'Chicken', 'Human / Other / Swine', 'Other / Swine', 'Bovine')) %>% 
  # filter(ISOLATIONS %in% c('Human / Swine', 'Swine', 'Chicken / Human', 'Cow / Human', 'Chicken', 'Human / Other / Swine', 'Other / Swine', 'Cow', 'Human / Turkey', 'Turkey', 'Cow / Swine')) %>% 
  mutate(ISOLATIONS=fct_rev(fct_infreq(f = ISOLATIONS,ordered = T ))) %>%
  ggplot(aes(x=ISOLATIONS, fill=continents)) + 
  geom_bar(width = .75, color = 'black') +
  coord_flip() + 
  ylab('Number of Transmission Clusters') + 
  xlab('Host Combinations') + 
  facet_wrap(~tree_clade, scales = 'free_x') + 
  scale_fill_brewer(palette = 'Set1', direction = 1) + 
  theme(strip.background = element_blank(), 
        strip.text = element_blank(), 
        legend.position = 'bottom', 
        legend.title = element_blank()) + 
  guides(fill=guide_legend(nrow=3))

ggsave('Figure8.jpeg', width = 9, height = 5, units = 'in', 
       bg='white')


# checking brads numbers #
transmission_summary %>% colnames()

transmission_summary %>%  group_by(ISOLATIONS) %>% tally() %>% arrange(desc(n))
master_meta
###
# human_swine
# swine
# chicken_human
# cow_human
# human_other_swine
# chicken
# other_swine
# #
# cow
# human_turkey
# turkey
# cow_swine


###

# other host transmission patterns
transmission_summary %>%
  mutate(ISOLATIONS=fct_rev(fct_infreq(f = ISOLATIONS,ordered = T ))) %>%
  filter(!(ISOLATIONS %in% c('Human', 'Human_Other', 'Other', 'Human_Swine', 'Swine', 'Chicken_Human', 'Cow_Human', 'Chicken', 'Human_Other_Swine', 'Other_Swine'))) %>% 
  ggplot(aes(x=ISOLATIONS, fill=continents)) + 
  geom_bar(width = .75, color= 'black') +
  coord_flip() + 
  ylab('number of transmission clusters') + 
  xlab('Host combinations') + 
  ggtitle('Number of transmission clusters associated with various host patterns', 
          'other host combinations') + 
  theme_cowplot() + 
  theme(panel.grid.major = element_line(colour = 'grey'))+
  facet_wrap(~tree_clade, scales = 'free_x')




transmission_summary %>%
  filter(percent_PDS != 1) %>% 
  arrange(desc(num_ISO), desc(num_isolates))


### multi ag host trans clust ###


master_meta %>%  
  select(genome, ISOLATION, PDS_acc, PHY_CLUST) %>% 
  ungroup() %>% 
  filter(!is.na(PHY_CLUST)) %>% 
  filter(ISOLATION != 'human') %>%
  group_by(PHY_CLUST) %>%
  mutate(num_iso=length(unique(ISOLATION))) %>% 
  filter(num_iso >1)

### multi country trans clust ###


LOOK <- 
  master_meta %>%  
  select(genome, country, PDS_acc, PHY_CLUST) %>% 
  ungroup() %>% 
  filter(!is.na(PHY_CLUST)) %>% 
  group_by(PHY_CLUST) %>%
  mutate(num_iso=length(unique(country))) %>% 
  filter(num_iso >1) %>% 
  arrange(desc(num_iso))







# mm %>% filter(PDS_acc == 'PDS000076482.55') %>% count(minsame)
# 
# library(castor)
# 
# # pretty sure this is what I want, but need way to figure out which 'cluster' original tips belong to.
# big_guy_2 <- collapse_tree_at_resolution(tree = big_guy, resolution = 0)
# 
# 
# # along these lines, maybe worthwhile to look into pangenome dereplication?
# 
# sum(table(big_guy_2$collapsed_nodes) != 1)
# #### the big PDS is taking forever.... ###
# 
# # PDS000076482.55.newick
# # can I split it in a reasonable way?
# library(ggtree)
# ggtree(big_guy, layout = 'circular')
# 
# hist(log(big_guy$edge.length))
# 
# 
# big_guy <- read.tree('./FIXED_SNP_TREES/PDS000076482.55.newick')
# 
# big_guy$tip.label <- 1:length(big_guy$tip.label)
# plot(big_guy)
# 
# big_PDS_mm <- mm %>% filter(PDS_acc %in% 'PDS000076482.55')
# 
# big_PDS_mm %>% group_by(targ_presence_MDR) %>% tally()
# 
# big_PDS_mm %>% group_by(continent) %>% tally()
# 
# 
# 
# 
