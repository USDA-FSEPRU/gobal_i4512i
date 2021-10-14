library(tidyverse)


read_blast_results <- function(BLAST_RES_PATH){
  all_hits <- read_tsv(BLAST_RES_PATH,
                       col_names = c('qseqid', 'sseqid', 'pident',
                                     'length', 'mismatch', 'gapopen',
                                     'qstart', 'qend', 'sstart', 'send',
                                     'slen','evalue', 'bitscore', 'qcovs', 'qcovhsp')) %>% 
    mutate(genome=sub('(.*)_[0-9]+', '\\1', qseqid))
  
  
  
}

calc_blast_cov1 <- function(all_hits, top_n_hits=5, feature_length, PREFIX){
  
  all_hits <- all_hits %>%
    filter(length > 499) %>% 
    group_by(genome) %>% 
    top_n(top_n_hits, length) %>%
    ungroup()
  
  all_nuc <- all_hits %>% mutate(start=ifelse(sstart>send, send, sstart), 
                                 end  =ifelse(sstart > send, sstart, send)) %>% 
    select(genome, qseqid, start, end, slen) %>% unique()
  
  
  all_nuc <- as.data.frame(all_nuc)
  
  
  res <- list()
  for (gen in unique(all_nuc$genome)){
    OG <- rep(F,feature_length)
    # print(gen)
    dat <- all_nuc[all_nuc$genome == gen,]
    
    for (row_num in 1:nrow(dat)){
      # print(row_num)
      covered <- c(rep(F,dat[row_num,3]-0), rep(T, dat[row_num,4]-dat[row_num,3]), rep(F, dat[row_num,5] - dat[row_num,4]))
      OG <- OG | covered
    }
    res[[gen]] <- sum(OG)
  }
  # browser()
  MDR_cov <- data.frame(genome=names(res), 
                        tot_cov=unlist(res)) %>% 
    mutate(perc_cov=tot_cov/feature_length)
  
  colnames(MDR_cov)[2:3] <- paste(PREFIX, colnames(MDR_cov)[2:3], sep = '_')
  
  return(MDR_cov)
  # hist(MDR_cov$tot_cov)
  
}


#######

master_meta <-
  read_tsv('./output/02_meta_Jan.tsv', guess_max = 14000)

all_genomes <- master_meta %>% select(genome)




# SGI4

USDA15WA1_SGI4 <- read_blast_results('./Jan_21_update/blast_results/USDA15WA-1_SGI4.results')
USDA15WA1_SGI4_cov <- calc_blast_cov1(all_hits = USDA15WA1_SGI4,
                            top_n_hits = 5, 
                            feature_length = 80661, 
                            PREFIX = 'SGI4')

# MDR



USDA15WA1_MDR <- read_blast_results('./Jan_21_update/blast_results/USDA15WA-1_MDR.results')
USDA15WA1_MDR_cov <- calc_blast_cov1(all_hits = USDA15WA1_MDR,
                                      top_n_hits = 5, 
                                      feature_length = 28188, 
                                      PREFIX = 'MDR')



thrW18kb_GI <- read_blast_results('./Jan_21_update/blast_results/USDA15WA-1_18kb_GI_thrW.results')
thrW18kb_cov <- calc_blast_cov1(all_hits = thrW18kb_GI,
                            top_n_hits = 5, 
                            feature_length = 18391, 
                            PREFIX = 'thrW18kb_GI')



#

USDA15WA1_PROfrag_581bp <- read_blast_results(BLAST_RES_PATH = './Jan_21_update/blast_results/USDA15WA-1_581bp_39kbPRO_thrW.results')
USDA15WA1_PROfrag_cov <- 
  calc_blast_cov1(all_hits = USDA15WA1_PROfrag_581bp,
                  top_n_hits = 3,
                  feature_length = 581, 
                  PREFIX = 'USDA15WA1_581ProFrag')



SO4698_sopE_mTmV <- read_blast_results(BLAST_RES_PATH = './Jan_21_update/blast_results/SO4698-09_sopE_mTmV_phage.results')
SO4698_sopE_mTmV_cov <- calc_blast_cov1(all_hits = SO4698_sopE_mTmV, 
                                        top_n_hits = 3, 
                                        feature_length =723, 
                                        PREFIX = 'SO4698_sopE_mTmV')

SO4698_NOsopE_PRO <- read_blast_results(BLAST_RES_PATH = './Jan_21_update/blast_results/SO4698-09_PRO_separate_from_sopE.results')
SO4698_NOsopE_PRO_cov <- calc_blast_cov1(all_hits = SO4698_NOsopE_PRO, top_n_hits = 3, feature_length = 1605, PREFIX = 'SO4698_NOsopE_PRO')

### INSERTION TO DOUBLE CHECK BLAST RESULTS ###
# 
# 
# SGI4_NEW <- read_blast_results(BLAST_RES_PATH = './Jan_21_update/BLAST_DOUBLE_CHECK/CP040686_SGI4.results')
# SGI4_NEW_cov <- calc_blast_cov1(all_hits = SGI4_NEW, top_n_hits = 5, feature_length = 80687, PREFIX = 'CP040686_SGI4')
# 
# 
# 
# 
# MDR_NEW <- read_blast_results(BLAST_RES_PATH = './Jan_21_update/BLAST_DOUBLE_CHECK/CP040686_MDR.results')
# MDR_NEW_cov <- calc_blast_cov1(all_hits = MDR_NEW, top_n_hits = 5, feature_length = 28194, PREFIX = 'CP040686_MDR')
# 
# 
# 
# 
# All_blast_covs <-
#   list(USDA15WA1_SGI4_cov,
#        USDA15WA1_MDR_cov,
#        thrW18kb_cov,
#        USDA15WA1_PROfrag_cov,
#        SO4698_sopE_mTmV_cov,
#        SO4698_NOsopE_PRO_cov,
#        SGI4_NEW_cov,
#        MDR_NEW_cov) %>%
#   reduce(full_join) %>%
#   right_join(all_genomes)
# 
# 
# All_blast_long <-
#   All_blast_covs %>%
#   select(genome, ends_with('perc_cov')) %>%
#   pivot_longer(cols = ends_with('perc_cov'),
#                names_to='target',
#                values_to='perc_cov') %>%
#   mutate(target=sub('_perc_cov','',target),
#          perc_cov=ifelse(is.na(perc_cov), 0, perc_cov),
#          targ_presence=case_when(
#            perc_cov >= .9 ~ 'full',
#            perc_cov >= .5 ~ 'partial',
#            TRUE           ~ 'absent'
#          ))
# 
# 
# blast_hists <-
#   All_blast_long %>%
#   ggplot(aes(x=perc_cov))+
#   geom_vline(xintercept = .9, color='purple')+
#   geom_vline(xintercept = .5, color='orange')+
#   geom_histogram() +
#   facet_wrap(~target)  +
#   ggtitle('Histograms of percent coverage of genomic features of interest',
#           'purple line indicates cutoff for "full" (90%)\norange line indicates cutoff for "partial" (50%)') +
#   xlab('percent coverage')
# 
# ggsave(plot=blast_hists, filename = './output/blast_histograms.jpeg')
# 
# All_blast_long %>% filter(grepl('SGI4|MDR',target ))
# 
# All_blast_res <- All_blast_long %>%
#   filter(grepl('SGI4|MDR',target ))%>%
#   pivot_wider(names_from = target, values_from=c(perc_cov, targ_presence)) %>%
#   select(genome, starts_with('targ_pres')) %>%
#   select(genome, contains('SGI4'), contains('MDR'))
# 
# all(All_blast_res[,2] == All_blast_res[,3])
# all(All_blast_res[,4] == All_blast_res[,5])
# all the same.

### END INSERTION ###


All_blast_covs <- 
  list(USDA15WA1_SGI4_cov,
     USDA15WA1_MDR_cov,
     thrW18kb_cov,
     USDA15WA1_PROfrag_cov,
     SO4698_sopE_mTmV_cov,
     SO4698_NOsopE_PRO_cov) %>%
  reduce(full_join) %>% 
  right_join(all_genomes)


All_blast_long <- 
  All_blast_covs %>%
  select(genome, ends_with('perc_cov')) %>% 
  pivot_longer(cols = ends_with('perc_cov'),
               names_to='target', 
               values_to='perc_cov') %>% 
  mutate(target=sub('_perc_cov','',target), 
         perc_cov=ifelse(is.na(perc_cov), 0, perc_cov), 
         targ_presence=case_when(
           perc_cov >= .9 ~ 'full', 
           perc_cov >= .5 ~ 'partial', 
           TRUE           ~ 'absent'
         ))


blast_hists <- 
  All_blast_long %>%
  ggplot(aes(x=perc_cov))+
  geom_vline(xintercept = .9, color='purple')+
  geom_vline(xintercept = .5, color='orange')+
  geom_histogram() +
  facet_wrap(~target)  +
  ggtitle('Histograms of percent coverage of genomic features of interest', 
          'purple line indicates cutoff for "full" (90%)\norange line indicates cutoff for "partial" (50%)') + 
  xlab('percent coverage')

ggsave(plot=blast_hists, filename = './output/blast_histograms.jpeg')


#### SPREAD, MERGE INTO MASTER_META ####
All_blast_res <- All_blast_long %>%
  pivot_wider(names_from = target, values_from=c(perc_cov, targ_presence)) %>% 
  select(genome, starts_with('targ_pres'))



####### merge in blast features to master meta

master_meta <-
  read_tsv('./output/02_meta_Jan.tsv', guess_max = 14000) %>%
  left_join(All_blast_res) 

master_meta %>% write_tsv('./output/03_meta_JAN.tsv')

# master_meta <- 
#   master_meta %>% 
#   mutate(across(.cols = starts_with('targ_presence'), 
#          ~ifelse(is.na(.x), 'absent', .x)))
# 




## merging island identification with mastermeta ##
# 
# all_islands <- merge(METAL_hits, ABX_hits, by = 'ID', all=TRUE)
# 
# colnames(all_islands) <- c('ID', 'SGI4_pident', 'SGI4_length', 'SGI4_class', 'MDR_pident', 'MDR_length', 'MDR_class')
# 
# 
# master_meta <- master_meta %>% left_join(all_islands)
# # master_meta <- merge(master_meta, all_islands, by='ID', all=TRUE)
# 
# # ONLY SRA GENOMES FROM HERE ON
# master_meta <- master_meta %>% filter(!grepl('GCF', ID))
# 
# master_meta$SGI4_class[is.na(master_meta$SGI4_class)] <- 'absent'
# master_meta$MDR_class[is.na(master_meta$MDR_class)] <- 'absent'
# 
# master_meta$MDR_length[is.na(master_meta$MDR_length)] <- 0
# master_meta$SGI4_length[is.na(master_meta$SGI4_length)] <- 0
# 
# # if putting fragment back in these need to be adjusted
# master_meta$SGI4_class <- factor(master_meta$SGI4_class, levels = c('absent', 'partial', 'full'))
# master_meta$MDR_class <- factor(master_meta$MDR_class, levels = c('absent', 'partial', 'full'))
# 
# master_meta$both_class <- paste(master_meta$SGI4_class , master_meta$MDR_class, sep=' / ')
# 
# #### sequence typing #####
# # used tseemanns mlst program. uses classic 7 gene typing scheme
# 
# seqtypes <- read_tsv('./data/ALL.seqtypes', col_names = c('file', 'scheme', 'ST', 'aroC', 'dnaN', 'hemD', 'hisD', 'purE', 'sucA', 'thrA'))
# seqtypes <- seqtypes %>%
#   mutate(ID=sub('.fasta','',file)) %>%
#   select(ID, ST)
# 
# master_meta <- master_meta %>% left_join(seqtypes)
# 
# 
# seqtypes %>% group_by(ST) %>% tally() %>% arrange(desc(n))
# #####################
# ### ADD IN ABRICATE
# 
# 
# 
# abricate_cols <- 
#   c('FILE', 'SEQUENCE', 'START', 'END', 'STRAND', 'GENE', 'COVERAGE', 'COVERAGE_MAP', 
#     'GAPS', 'P_COVERAGE', 'P_ID', 'DATABASE', 'ACCESSION', 'PRODUCT', 'RESISTANCE')
# 
# megares <-
#   read_tsv('./data/abricate_res/MEGARES.ALL', col_names = abricate_cols, comment = '#')
# 
# ncbi <- 
#   read_tsv('./data/abricate_res/NCBI.ALL', col_names = abricate_cols, comment = '#')
# 
# plasmidfinder <- 
#   read_tsv('./data/abricate_res/PLASMID.ALL', col_names = abricate_cols, comment = '#')
# 
# vfdb <- 
#   read_tsv('./data/abricate_res/VFDB.ALL', col_names = abricate_cols, comment = '#')
# 
# viroseqs <-
#   read_tsv('./data/abricate_res/VIRO.ALL', col_names = abricate_cols, comment = '#')
# 
# all_abricate <- bind_rows(megares, ncbi, plasmidfinder, vfdb, viroseqs)
# 
# 
# ncbi <- ncbi %>% mutate(genome=sub('./(.*).fasta','\\1',FILE)) %>% select(genome, everything(), -FILE)
# 
# 
# res_types <- 
#   ncbi %>%
#   separate_rows(RESISTANCE, sep = ';') %>% 
#   group_by(genome) %>%
#   arrange(RESISTANCE) %>% 
#   summarise(res_type=paste(unique(RESISTANCE), sep = '~', collapse = '~'), 
#             num_res =length(RESISTANCE)) %>% 
#   write_tsv('./output/04_restypes.tsv')
# 
# res_types %>% group_by(res_type, num_res) %>% tally() %>% arrange(desc(n))
# 
# res_types_summary <- res_types %>%
#   group_by(res_type) %>%
#   tally() %>%
#   arrange(desc(n))
# 
# 
# 
# 
# ### megares stuff ###
# 
# megares <- megares %>% mutate(TYPE=sub('(.*)\\|(.*)','\\1',PRODUCT))
# unique(megares$TYPE)
# 
# check <- megares %>% filter(TYPE == 'Multi-compound')
# ### colistin ###
# 
# grep('COLISTIN', ncbi$RESISTANCE)
# 
# 
# colistin_res <- ncbi[grep('COLISTIN', ncbi$RESISTANCE),]
# 
# qnrs <- ncbi[grep('qnr', ncbi$GENE),]
# 
# 
# #####
# master_meta <- master_meta %>%
#   mutate(continent=case_when(
#     country %in% c('USA', 'Mexico', 'Canada')  ~ 'North America',
#     country %in% c('Thailand', 'Taiwan', 'China')  ~ 'Asia', 
#     country %in% c('Brazil', 'Trinidad and Tobago', 'Peru') ~ 'South America', 
#     country %in% c('United Kingdom', 'Germany', 'Denmark')  ~ 'Europe', 
#     TRUE ~ 'unknown'
#   ))
# 
# ###
# ######### TSNE / UMAP #####
# 
# all_dist <- funfuns::read_mash_tri('./data/mash_files/passing_refs.dist')
# 
# # ALL_list <- read_tsv('./mash_files/ALL.dist', col_names = c('genome1', 'genome2', 'distance', 'pval', 'kmers'))
# # 
# # all_dist <- ALL_list %>% 
# #   select(genome1, genome2, distance) %>% 
# #   mutate(genome1=sub('.fasta','',genome1), 
# #          genome2=sub('.fasta','',genome2)) %>% 
# #   filter(genome1 %in% master_meta$ID) %>% 
# #   filter(genome2 %in% master_meta$ID) %>% 
# #   spread(key=genome2, value=distance) %>%
# #   column_to_rownames(var = 'genome1') %>%
# #   as.dist()
# library(Rtsne)
# all_tsne <- Rtsne(is_distance=TRUE, all_dist)
# 
# 
# tsne_tib <- tibble(ID=attributes(all_dist)$Labels, T1=all_tsne$Y[,1], T2=all_tsne$Y[,2]) %>% 
#   mutate(ID=sub('(.*).fasta','\\1',ID))
# 
# 
# 
# master_meta <- master_meta %>% full_join(tsne_tib)
# 
# 
# 
# #####
# 
# library(uwot)
# 
# umap_test <- umap(all_dist)
# plot(umap_test)
# 
# umap_tib <- tibble(ID=attributes(all_dist)$Labels, 
#                    U1=umap_test[,1], 
#                    U2=umap_test[,2]) %>% 
#   mutate(ID=sub('.fasta','',ID))
# 
# 
# 
# 
# master_meta <- master_meta %>% full_join(umap_tib)
# 
# # # clustering info #
# # 
# # 
# # ### SHOW THAT MASH BASED DREP IS NOT GOOD! ###
# # clust_info <- read_csv('passingQC_drep/data_tables/Cdb.csv')%>%
# #   mutate(ID=sub('.fasta','',genome), 
# #          drep_cluster=primary_cluster) %>% 
# #   select(ID, drep_cluster)
# # ### can prob remove this
# # 
# # master_meta <- master_meta %>% left_join(clust_info)
# # 
# 
# master_meta %>% write_tsv('./output/01_clean.tsv')
# 

# blast hits broken by contig ends #

# For calculating broken blast hits #
# contig_lens <- read_tsv('./data/ALL_CONTIG_LENGTHS.tsv', col_names=c('contig', 'len'))
# 
# 
# test <- all_hits %>% group_by(genome) %>% top_n(10, length)
# 
# test <- test %>%ungroup() %>% group_by(genome) %>%  summarise(sum_hits=sum(length))
# 
# 
# longest_hits <- all_hits %>% group_by(genome) %>%
#   slice(which.max(length)) %>%
#   left_join(test, by = 'genome') %>% 
#   ungroup()
# 
# 
# longest_hits %>% ggplot(aes(x=length, y=sum_hits)) + geom_point()
# 
# longest_hits <- longest_hits %>% 
#   left_join(contig_lens) %>%
#   mutate(hit_contig_diff = qend - contig_len, 
#          broken_assembly = ifelse(hit_contig_diff == 0, 'probably', 'no'), 
#          broken_assembly = ifelse(qstart == 1, 'probably', broken_assembly))
# 
# longest_hits %>% filter(length > 2000) %>% group_by(broken_assembly) %>% tally()
# 
# longest_hits %>% filter(length > 2000) %>% group_by(broken_assembly) %>%
#   summarise(mean_hitlen=mean(length), 
#             maxhitlen = max(length), 
#             min_hitlen=min(length))
# 
#### USE THESE WORDS ###

"The presence of the MDR module often breaks short read assemblies in this data.
4460 genomes in this dataset had blast hits longer than 2000 bp corresponding to the MDR module.  
3839 of these hits abutted contig ends, meaning the assembly was broken in this module.
621 genomes had hits that did not abut contig ends, but the mean length of these hits 
was 10482 meaning many were incomplete versions of the MDR module."

