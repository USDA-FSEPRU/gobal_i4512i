library(tidyverse)



# make this script a pangenome analysis
# table with stats for the pangenome
# GO enrich for persistent shell cloud
# clusters, GO enrich for clusters


# three clustering iterations
# first with only shell genes = shell
# second with shell and cloud = sheclo
# third with only cloud = cloud

# matr <- read_csv('ppanggolin_gbk_defrag/matrix.csv')
# matr2 <- read_tsv('ppanggolin_gbk_defrag/gene_presence_absence.Rtab')


matr2 <- read_tsv('Jan_21_update/ppangg_gbk_feb21/gene_presence_absence.Rtab')


ctyps <- c('cccdddcccccddd', rep('c', ncol(matr2) - 1) %>%
             paste(collapse = '')) %>% 
  paste(collapse = '')

matr <- read_csv('Jan_21_update/ppangg_gbk_feb21/matrix.csv',
                 col_types = ctyps, quote = '"',
                 quoted_na = T,
                 trim_ws = T)

# matr$GCA_005688305[2]
# sheclo <- matr[matr[,2] != 'persistent',] %>% pull(Gene)
# shells <- matr[matr[,2] == 'shell',] %>% pull(Gene)
clouds <- matr[matr[,2] == 'cloud',] %>% pull(Gene)
shell <- matr[matr[,2] == 'shell',] %>% pull(Gene)
persistent <- matr[matr[,2] == 'persistent',] %>% pull(Gene)


#### ADD IN PAN SUMMARY HERE


length(clouds)       # 17546
length(shell)        # 344
length(persistent)   # 4271

# persistent histogram : number of persistent genes per genome
matr2[1:10,1:10]

persist_matrix <- matr2 %>% filter(Gene %in% persistent) %>% column_to_rownames(var='Gene') %>% as.matrix()
shell_matrix <- matr2 %>% filter(Gene %in% shell) %>% column_to_rownames(var='Gene')%>% as.matrix()
cloud_matrix <- matr2 %>% filter(Gene %in% clouds) %>% column_to_rownames(var='Gene')%>% as.matrix()


hist(colSums(persist_matrix))
hist(colSums(shell_matrix))
hist(colSums(cloud_matrix))
# matr[,2]


# rm(matr)


cloud_mat <- matr2 %>% 
  filter(Gene %in% clouds) %>% 
  column_to_rownames(var = 'Gene') %>% 
  as.matrix() 


cloud_colsums <- colSums(cloud_mat)
hist(cloud_colsums, breaks = 100)  # is this the number of cloud genes each genome has?

sum((cloud_colsums) == 0)

cloud_mat <- cloud_mat[,colSums(cloud_mat) > 0]

# rntmp <- rownames(cloud_mat)
# cloud_mat <- apply(cloud_mat, 2, as.numeric)
# rownames(cloud_mat) <- rntmp
# rownames(cloud_mat)
dim(cloud_mat)
str(cloud_mat)
####### CONSIDER FILTERING GENES HERE #########
# KEEP IN MIND PERSISTENT GENES ARE ALREADY REMOVED #
cloud_mat_transpose <- t(cloud_mat)


# remove singleton cloud genes

clouds_filt <- cloud_mat_transpose[,colSums(cloud_mat_transpose) > 1]

# rm(clouds)
# rm(cloud_mat)
# rm(matr2)

hist(colSums(clouds_filt), breaks = 100) # number of genomes each cloud gene is detected in

sum(colSums(clouds_filt) == 2) # 2213 genes only found in 2 genomes.  
max(colSums(clouds_filt))      # 5613

print('constructing graph: Nodes are genomes, edges are number of shared genes')

# only run this chunk if the leiden clusters dont exist yet
# will need to update the file name as we move away from the "JAN_UPDATE" scheme
if(!(file.exists('./output/02_leiden_clust_JAN_UPDATE.csv'))){
  #cross product
  co_mat <- clouds_filt %*% t(clouds_filt)
  
  #co_mat is a co-occurance matrix, both rows and columns are genomes
  # numbers indicate how many genes those genomes share with eachother
  # set diagonal to 0
  diag(co_mat) <- 0
  #co_mat[1:10,1:10]
  library(igraph)
  
  g <- graph_from_adjacency_matrix(co_mat, mode = "upper", weighted = TRUE)
  clouv <- cluster_louvain(g)
  # welp there is a very frustrating issue with rstudio, reticulate, some c libraries and LD_LIB_PATH or something.
  # I'm just going to run it in python.
  
  clust_info <- tibble(genome = names(membership(clouv)),
                       louvain_cluster = membership(clouv)) 
  
  # clust_info %>% 
  #   write_csv('./output/02_louvain_clust_JAN_UPDATE.csv')
  # 
  
  clust_info %>% group_by(louvain_cluster) %>% tally() %>% arrange(desc(n))
  
  
  write.graph(graph = g, format = 'gml', file = "./output/JAN_UPDATE_cloud_co.gml")
  
  # then run leiden_clust.py
  
  # this will only run on the fsep11 server because my conda path is hardcoded
  system('./scripts/leiden_wrapper.sh')

}



leid <-  
  read_csv('./output/02_leiden_clust_JAN_UPDATE.csv') %>% 
  dplyr::select(genome, leiden_cluster)

master_meta <- read_tsv('./output/01_meta_JAN_UPDATE.tsv', guess_max = 14000) %>% 
  left_join(leid)

master_meta %>% group_by(leiden_cluster) %>% tally()

master_meta <- 
  master_meta %>%
  mutate(leiden_cluster=case_when(
    leiden_cluster == 0   ~ 'A',
    leiden_cluster == 1   ~ 'B',
    leiden_cluster == 2   ~ 'C',
    leiden_cluster == 3   ~ 'none',
    is.na(leiden_cluster) ~ 'none'
  ))

master_meta %>% group_by(leiden_cluster) %>% tally()

master_meta %>% write_tsv('./output/02_meta_Jan.tsv')

# 
# c0 <- leid$genome[leid$leiden_cluster == 0]
# c1 <- leid$genome[leid$leiden_cluster == 1]
# c2 <- leid$genome[leid$leiden_cluster == 2]
# c3 <- leid$genome[leid$leiden_cluster == 3]

# colnames(interpro1)


# 
# 
# c0_mat <- clouds_filt[rownames(clouds_filt) %in% c0,]
# c1_mat <- clouds_filt[rownames(clouds_filt) %in% c1,]
# c2_mat <- clouds_filt[rownames(clouds_filt) %in% c2,]
# c3_mat <- clouds_filt[rownames(clouds_filt) %in% c3,]

# 
# 
# clust_pres <- 
#   tibble(Gene=names(((colSums(c0_mat) / nrow(c0_mat)) * 100)), 
#          c0=((colSums(c0_mat) / nrow(c0_mat)) * 100), 
#          c1=((colSums(c1_mat) / nrow(c1_mat)) * 100), 
#          c2=((colSums(c2_mat) / nrow(c2_mat)) * 100))
# 
# 
# 
# clust_pres %>%
#   mutate(G50_anywhere=
#            case_when(c0 > 50 ~ TRUE, 
#                      c1 > 50 ~ TRUE, 
#                      c2 > 50 ~ TRUE, 
#                      TRUE ~ FALSE), 
#   ) 

# 
# 
# gene_clust_enich <- 
#   clust_pres %>%
#   pivot_longer(cols = starts_with('c'),
#                names_to = 'cluster',
#                values_to = 'abundance') %>% 
#   group_by(Gene) %>%
#   arrange(abundance) %>% 
#   summarise(max_clust=cluster[which.max(abundance)],
#             min_clust=cluster[which.min(abundance)],
#             ind1 = cluster[1], 
#             ind2=cluster[2], 
#             ind3=cluster[3], 
#             # ind4=cluster[4], 
#             max_abund=abundance[which.max(abundance)], 
#             min_abund=abundance[which.min(abundance)], 
#             max_m_min=max_abund - min_abund, 
#             max_m_2nd=abundance[3]-abundance[2], 
#             max_m_mean=abundance[3] - (abundance[2]+abundance[1])/2, 
#             max_d_mean=abundance[3] / ((abundance[2]+abundance[1])/2)+.1)
# 
# hist(gene_clust_enich$max_d_mean)
# 
# 
# # hist(tmp$max_m_2nd, breaks = 100)
# # hist(tmp$max_m_min, breaks = 100)
# 
# 
# sum(gene_clust_enich$max_m_2nd > 10)
# sum(gene_clust_enich$max_m_2nd > 25)
# sum(gene_clust_enich$max_m_2nd > 50)
# sum(gene_clust_enich$max_m_2nd > 75)

# gene was considered to be loosely associated with a cluster if the within cluster percent presence minus
# the percent presence of the next closest cluster was greater than 10%

# A gene was considered to be associated with a cluster if the within cluster percent presence minus
# the next closest cluster percent presence was 25%


# A gene was considered to be enriched within a cluster if the within cluster percent presence minus
# the next closest cluster percent presence was 50%


# A gene was considered to be strongly enriched within a cluster if the within cluster percent presence minus
# the next closest cluster percent presence was 75%
# 
# loose <- gene_clust_enich %>% filter(max_m_2nd >10)
# 
# assoc <- gene_clust_enich %>% filter(max_m_2nd >25)
# 
# enrich <- gene_clust_enich %>% filter(max_m_2nd >50)
# 
# stron_enrich <- gene_clust_enich %>% filter(max_m_2nd >75)
# 
# 
# hist(stron_enrich$max_abund)
# 
# stron_enrich %>% group_by(max_clust) %>% tally()
# 
# 
# enrich %>% group_by(max_clust) %>% tally()
# 

# c1_strong <- stron_enrich %>%filter(max_clust =='c1')%>% transmute(accno=Gene)
# c3_strong <- stron_enrich %>%filter(max_clust =='c3') %>% transmute(accno=Gene)

##### NEED TO CHANGE THINGS HERE #####
# 
# 
# c0_enrich <- enrich %>%filter(max_clust =='c0') %>% pull(Gene)
# c1_enrich <- enrich %>%filter(max_clust =='c1') %>% pull(Gene)
# c2_enrich <- enrich %>%filter(max_clust =='c2') %>% pull(Gene)
# 
# 
# c0_asso <- assoc %>%filter(max_clust =='c0') %>% pull(Gene)
# c1_asso <- assoc %>%filter(max_clust =='c1') %>% pull(Gene)
# c2_asso <- assoc %>%filter(max_clust =='c2') %>% pull(Gene)
# 
# 
# 
# c0_loose <- loose %>%filter(max_clust =='c0') %>% pull(Gene)
# c1_loose <- loose %>%filter(max_clust =='c1') %>% pull(Gene)
# c2_loose <- loose %>%filter(max_clust =='c2') %>% pull(Gene)
# 





# 
# library(funfuns)
# 
# # library(topGO)
# 
# c1_GO <- 
#   list(topGO_wrapper(myInterestingGenes = c1_asso, mapping_file = './output/JAN_UPDATE_GO_gene_universe.tsv', ont = 'BP'),
#        topGO_wrapper(myInterestingGenes = c1_asso, mapping_file = './output/JAN_UPDATE_GO_gene_universe.tsv', ont = 'MF'),
#        topGO_wrapper(myInterestingGenes = c1_asso, mapping_file = './output/JAN_UPDATE_GO_gene_universe.tsv', ont = 'CC')) %>% 
#   bind_rows() %>%
#   filter(pval < .05) %>% 
#   write_tsv('./output/JAN_cloud_clust_1_GO.tsv')
# 
# 
# 
# c2_GO <- 
#   list(topGO_wrapper(myInterestingGenes = c2_asso,mapping_file = './output/JAN_UPDATE_GO_gene_universe.tsv', ont = 'BP'), 
#        topGO_wrapper(myInterestingGenes = c2_asso,mapping_file = './output/JAN_UPDATE_GO_gene_universe.tsv', ont = 'MF'),
#        topGO_wrapper(myInterestingGenes = c2_asso,mapping_file = './output/JAN_UPDATE_GO_gene_universe.tsv', ont = 'CC'))%>% 
#   bind_rows() %>% filter(pval < .05) %>% 
#   write_tsv('./output/JAN_cloud_clust_2_GO.tsv')
# 
