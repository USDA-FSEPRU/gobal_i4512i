library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(stringr)


theme_set(theme_cowplot(font_family = 'arial'))
theme_update(panel.grid.major=element_line(color='grey80'), 
             panel.border = element_rect(color='black'),
             axis.title=element_text(face = 'bold'))

master_meta <- read_tsv('./output/05_JAN_meta.tsv') %>% 
  filter(!is.na(tree_clade)) %>% 
  filter(tree_clade != 'other') %>%
  mutate(tree_clade=factor(tree_clade, levels=c('clade_one', 'clade_two', 'clade_three')),
         `Tree Clade`=fct_recode(tree_clade, `Clade 1`='clade_one', `Clade 2`='clade_two', `Clade 3`='clade_three'),
         ISOLATION=str_to_title(ISOLATION), 
         ISOLATION=case_when(
           ISOLATION == 'Cow'   ~  'Bovine', 
           TRUE              ~   ISOLATION), 
         Source=factor(ISOLATION, levels =c('Human','Bovine', 'Chicken', 'Other', 'Swine', 'Turkey')))


master_meta$Source
master_meta$`Tree Clade`
source_colors <- brewer.pal(length(levels(factor(master_meta$ISOLATION))), name = 'Set1')
names(source_colors) <- levels(factor(master_meta$ISOLATION))

# IF TURKEY COLOR PROBLEM ADDRESS HERE
# source_colors['Turkey'] <- '#00FFFF'
##

######## NEWFIGS ########
master_meta %>%
  filter(year >2003) %>% 
  count(year, targ_presence_SGI4) %>% 
  group_by(year) %>% 
  mutate(tot_genomes_year=sum(n), 
         perc_genomes_year=(n/tot_genomes_year)*100, 
         `SGI-4`=factor(str_to_title(targ_presence_SGI4), levels = c('Full', 'Partial', 'Absent'))) %>% 
  ggplot(aes(x=year, y=perc_genomes_year, fill=`SGI-4`)) + 
  geom_col(color='white') + 
  geom_text(aes(y=50, label=tot_genomes_year), size=3) + 
  ylab('Percent of Genomes Containing SGI-4') +
  xlab('Year') + 
  scale_fill_brewer(palette = 'Set1')

ggsave(filename = 'Figure4B.jpeg', width = 9, height = 5, units = 'in', bg='white')

master_meta %>%
  filter(year >2003) %>% 
  filter(tree_clade== 'clade_one') %>% 
  count(year, targ_presence_SGI4) %>% 
  group_by(year) %>% 
  mutate(tot_genomes_year=sum(n), 
         perc_genomes_year=(n/tot_genomes_year)*100, 
         `SGI-4`=factor(str_to_title(targ_presence_SGI4), levels = c('Full', 'Partial', 'Absent'))) %>% 
  ggplot(aes(x=year, y=perc_genomes_year, fill=`SGI-4`)) + 
  geom_col(color='white') + 
  geom_text(aes(y=50, label=tot_genomes_year), size=3) + 
  ylab('Percent of Clade 1 Genomes Containing SGI-4') +
  xlab('Year')+
  scale_fill_brewer(palette = 'Set1')


ggsave(filename = 'Figure4C.jpeg', width = 9, height = 5, units = 'in', bg='white')

master_meta %>%
  # filter(year >2003) %>% 
  # filter(tree_clade== 'clade_one') %>% 
  count(country, targ_presence_SGI4) %>% 
  group_by(country) %>% 
  mutate(tot_genomes_country=sum(n), 
         perc_genomes_country=(n/tot_genomes_country)*100, 
         `SGI-4`=factor(str_to_title(targ_presence_SGI4), levels = c('Full', 'Partial', 'Absent'))) %>% 
  filter(tot_genomes_country > 4) %>% 
  ggplot(aes(x=country, y=perc_genomes_country, fill=`SGI-4`)) + 
  geom_col(color='white') + 
  geom_text(aes(y=50, label=tot_genomes_country), size=3) + 
  ylab('Percent of Genomes Containing SGI-4') +
  xlab('Year') + 
  theme(axis.text.x = element_text(angle=-45, hjust=.1)) + 
  scale_fill_brewer(palette = 'Set1')


ggsave(filename = 'Figure4A.jpeg', width = 9, height = 5, units = 'in', bg='white')



# MDR NOW

master_meta %>%
  filter(year >2003) %>% 
  count(year, targ_presence_MDR) %>% 
  group_by(year) %>% 
  mutate(tot_genomes_year=sum(n), 
         perc_genomes_year=(n/tot_genomes_year)*100, 
         `MDR`=factor(str_to_title(targ_presence_MDR), levels = c('Full', 'Partial', 'Absent'))) %>% 
  ggplot(aes(x=year, y=perc_genomes_year, fill=`MDR`)) + 
  geom_col(color='white') + 
  geom_text(aes(y=50, label=tot_genomes_year), size=3) + 
  ylab('Percent of Genomes Containing MDR') +
  xlab('Year') + 
  scale_fill_brewer(palette = 'Set1')


ggsave(filename = 'Figure5B.jpeg', width = 9, height = 5, units = 'in', bg='white')

master_meta %>%
  filter(year >2003) %>% 
  filter(tree_clade== 'clade_one') %>% 
  count(year, targ_presence_MDR) %>% 
  group_by(year) %>% 
  mutate(tot_genomes_year=sum(n), 
         perc_genomes_year=(n/tot_genomes_year)*100, 
         `MDR`=factor(str_to_title(targ_presence_MDR), levels = c('Full', 'Partial', 'Absent'))) %>% 
  ggplot(aes(x=year, y=perc_genomes_year, fill=`MDR`)) + 
  geom_col(color='white') + 
  geom_text(aes(y=50, label=tot_genomes_year), size=3) + 
  ylab('Percent of Clade 1 Genomes Containing MDR') +
  xlab('Year') + 
  scale_fill_brewer(palette = 'Set1')


ggsave(filename = 'Figure5C.jpeg', width = 9, height = 5, units = 'in', bg='white')

master_meta %>%
  # filter(year >2003) %>% 
  # filter(tree_clade== 'clade_one') %>% 
  count(country, targ_presence_MDR) %>% 
  group_by(country) %>% 
  mutate(tot_genomes_country=sum(n), 
         perc_genomes_country=(n/tot_genomes_country)*100, 
         `MDR`=factor(str_to_title(targ_presence_MDR), levels = c('Full', 'Partial', 'Absent'))) %>% 
  filter(tot_genomes_country > 4) %>% 
  ggplot(aes(x=country, y=perc_genomes_country, fill=`MDR`)) + 
  geom_col(color='white') + 
  geom_text(aes(y=50, label=tot_genomes_country), size=3) + 
  ylab('Percent of Genomes Containing MDR') +
  xlab('Year') + 
  theme(axis.text.x = element_text(angle=-45, hjust=.1)) + 
  scale_fill_brewer(palette = 'Set1')


ggsave(filename = 'Figure5A.jpeg', width = 9, height = 5, units = 'in', bg='white')





# CO-OCCUR
# master_meta %>% 
#   # group_by(targ_presence_SGI4, targ_presence_MDR) %>% 
#   mutate(`SGI-4 / MDR` = paste(targ_presence_SGI4, targ_presence_MDR, sep = ' / ')) %>% 
#   group_by()



master_meta %>% 
  mutate(`SGI-4 / MDR` = paste(str_to_title(targ_presence_SGI4), str_to_title(targ_presence_MDR), sep = ' / '), 
         `SGI-4 / MDR` = factor(`SGI-4 / MDR`, levels = c(
           'Full / Full', 'Full / Partial', 'Full / Absent', 
           'Partial / Full', 'Partial / Partial', 'Partial / Absent', 
           'Absent / Full', 'Absent / Partial', 'Absent / Absent'))) %>% 
  filter(year >2003) %>% 
  filter(tree_clade== 'clade_one') %>% 
  count(year, `SGI-4 / MDR`) %>% 
  group_by(year) %>% 
  mutate(tot_genomes_year=sum(n), 
         perc_genomes_year=(n/tot_genomes_year)*100) %>%  
  ggplot(aes(x=year, y=perc_genomes_year, fill=`SGI-4 / MDR`)) + 
  geom_col(color='white') + 
  geom_text(aes(y=50, label=tot_genomes_year), size=3) + 
  ylab('Percent of Clade 1 Genomes Containing \nSGI-4 and MDR') +
  xlab('Year') + 
  scale_fill_brewer(palette = 'Set1')

ggsave(filename = 'Figure6B.jpeg', width = 9, height = 5, units = 'in', bg='white')

master_meta %>% 
  mutate(`SGI-4 / MDR` = paste(str_to_title(targ_presence_SGI4), str_to_title(targ_presence_MDR), sep = ' / '), 
         `SGI-4 / MDR` = factor(`SGI-4 / MDR`, levels = c(
           'Full / Full', 'Full / Partial', 'Full / Absent', 
           'Partial / Full', 'Partial / Partial', 'Partial / Absent', 
           'Absent / Full', 'Absent / Partial', 'Absent / Absent'))) %>% 
  filter(tree_clade == 'clade_one') %>%
  count(country, `SGI-4 / MDR`) %>% 
  group_by(country) %>% 
  mutate(tot_genomes_country=sum(n), 
         perc_genomes_country=(n/tot_genomes_country)*100, 
         `MDR`=factor(str_to_title(`SGI-4 / MDR`), levels = c('Full', 'Partial', 'Absent'))) %>% 
  filter(tot_genomes_country > 4) %>% 
  ggplot(aes(x=country, y=perc_genomes_country, fill=`SGI-4 / MDR`)) + 
  geom_col(color='white') + 
  geom_text(aes(y=50, label=tot_genomes_country), size=3) + 
  ylab('Percent of Clade 1 Genomes Containing \nSGI-4 and MDR') +
  xlab('') +
  theme(axis.text.x = element_text(angle=-45, hjust=.1)) + 
  scale_fill_brewer(palette = 'Set1')


ggsave(filename = 'Figure6A.jpeg', width = 9, height = 5, units = 'in', bg='white')







# FIG2
master_meta %>%
  filter(year>2003) %>% 
  count(`Tree Clade`, year) %>% 
  group_by(year) %>% 
  mutate(tot_genomes_year=sum(n), 
         perc_genomes_year=(n/tot_genomes_year) *100) %>% 
  ggplot(aes(x=year, y=perc_genomes_year, fill=`Tree Clade`))+
  geom_col(color='white') + 
  geom_text(aes(y=50, label=tot_genomes_year), size=3) + 
  ylab('Percent of Total Genomes per Year') + 
  scale_fill_manual(values = c('steelblue', 'purple', 'forestgreen'))

ggsave(filename = 'Figure2.jpeg', width = 9, height = 5, units = 'in', bg='white')

### CHANGED ISOLATION to SOURCE HERE
tree_clade_iso_totals <- 
  master_meta %>%
  count(tree_clade, Source, .drop=F,name="tot_n_genomes") 



tree_clade_totals <- 
  master_meta %>%
  group_by(tree_clade) %>%
  summarise(tot_n_genomes=n())


tree_clade_year_totals <- 
  master_meta %>%
  group_by(tree_clade, year) %>%
  summarise(tot_n_genomes=n())


P1_tree_clade_year_totals <- 
  tree_clade_year_totals %>% 
  filter(year < 2021) %>% 
  ggplot(aes(x=year, y=tot_n_genomes)) + 
  geom_line() +
  facet_wrap(~tree_clade)


# fig 3
# P1_tree_clade_year_totals + ylab('Total Number of Genomes') + 
#   theme(strip.background = element_blank(), 
#         strip.text = element_blank())
# 
# 
# ggsave('fig3_1.jpg', width = 9, height = 5, units = 'in', bg='white')
# 

# P1 alternate

p1A <- 
  master_meta %>%
  filter(year < 2021) %>% 
  count(year) %>%
  ggplot(aes(x=year, y=n)) +
  geom_line() + 
  ylab('Total Number of Genomes')

p1A
ggsave('FigureS1.jpg', width = 9, height = 5, units = 'in', bg='white')




###


tree_clade_year_iso_totals <-
  master_meta %>%
  filter(year < 2021) %>% 
  count(tree_clade, Source, year)

P2_tree_clade_year_iso_totals <- 
  tree_clade_year_iso_totals %>%
  ggplot(aes(x=year, y=n, color=Source)) + 
  geom_line() +
  facet_wrap(~tree_clade, scales = 'free_y') + 
  scale_color_manual(values = source_colors) + 
  ylab('Total Number of Genomes') + 
  theme(strip.background = element_blank(), 
        strip.text=element_blank())

P2_tree_clade_year_iso_totals

ggsave('figure3A.jpg', width = 9, height = 5, units = 'in', bg='white')

source_colors2 <- source_colors[-c(3,4)]

P3_tree_clade_year_iso_no_hum <- 
  tree_clade_year_iso_totals %>% 
  filter(!(Source %in% c('Human', 'Other'))) %>% 
  ggplot(aes(x=year, y=n, color=Source)) + 
  geom_line() +
  facet_wrap(~tree_clade, scales = 'free_y') + 
  scale_color_manual(values=source_colors2) + 
  ylab('Total Number of Genomes') + 
  theme(strip.background = element_blank(), 
        strip.text=element_blank())

P3_tree_clade_year_iso_no_hum

ggsave('figure3B.jpg', width = 9, height = 5, units = 'in', bg='white')

# 
# P4_tree_clade_iso_totals <- 
#   master_meta %>%
#   filter(!(ISOLATION %in% c('Human', 'Other'))) %>% 
#   count(continent, tree_clade, ISOLATION) %>% 
#   ggplot(aes(x=continent, y=n, fill=ISOLATION)) +
#   geom_col() +
#   facet_wrap(~tree_clade, scales = 'free_y') + 
#   theme(axis.text.x = element_text(angle=-45, hjust = .1)) + 
#   scale_fill_brewer(palette = 'Set1') + 
#   ylab('number of genomes')
# 
# P4_tree_clade_iso_totals

# ggsave('P4.jpg', width = 9, height = 5, units = 'in')



# master_meta %>% filter(ISOLATION == 'other') %>% select(isolation_source) %>% table()

######## AMR stress section #######
AMR_reference <- read_tsv('./Jan_21_update/ReferenceGeneCatalog.txt')

# 
# AMR_reference %>% 
#   pivot_longer(cols = c(allele, gene_family), names_to='NAMES_TO', values_to="VALUES_TO")



GENE_2_CLASS_DICT <- 
  AMR_reference %>%
  dplyr::select(allele, gene_family, scope:subclass) %>% 
  pivot_longer(cols = c(allele, gene_family),
               names_to='ORIGIN',
               values_to="GENE_ID") %>% 
  filter(!is.na(GENE_ID)) %>%
  unique()
  


AMR_finder <- 
  master_meta %>% 
  dplyr::select(genome, AMR_genotypes, stress_genotypes, virulence_genotype) %>% 
  transmute(genome=genome, 
            AMR_ID=gsub('"','',AMR_genotypes), 
            stress_ID=gsub('"','',stress_genotypes),
            virulence_ID=gsub('"','',virulence_genotype))
  


# separate_rows(ends_with('ID'), sep = ',') %>% 
  # filter(!grepl('PARTIAL', AMR_ID)) %>% 
  # filter(!grepl('PARTIAL', virulence_ID)) %>% 
  # filter(!grepl('PARTIAL', stress_ID))

# COUNTS OF AMR CLASSES HERE
AMR <- 
  AMR_finder %>%
  transmute(genome=genome, 
            GENE_ID=AMR_ID) %>% 
  separate_rows(GENE_ID, sep = ',') %>% 
  filter(!grepl('PARTIAL', GENE_ID)) %>% 
  mutate(GENE_ID=sub('(.*)=.*', '\\1', GENE_ID)) %>% 
  left_join(GENE_2_CLASS_DICT) 

AMR <- master_meta %>%
  dplyr::select(genome, tree_clade ,Source, leiden_cluster, year, country) %>%
  right_join(AMR)  %>% 
  separate_rows(class, sep = '/') %>%
  unique()



unique(AMR$class)

# AMR %>%
#   select(genome, class) %>% 
#   separate_rows(class, sep = '/') %>% 
#   unique() %>% 
#   mutate(PRESENT = 1) %>% 
#   spread(key = class, value=PRESENT, fill = 0)



AMR_sum <- 
  AMR %>%
  separate_rows(class, sep = '/') %>% 
  group_by(genome) %>%
  arrange(class) %>%
  summarise(all_AMR_genes=paste(unique(GENE_ID), collapse = '_'), 
            AMR_classes=paste(unique(class[!is.na(class)]), collapse='_'), 
            num_AMR_classes=length(unique(class[!is.na(class)])))



VIR <- 
  AMR_finder %>%
  transmute(genome=genome, 
            GENE_ID=virulence_ID) %>% 
  separate_rows(GENE_ID, sep = ',') %>% 
  filter(!grepl('PARTIAL', GENE_ID)) %>% 
  mutate(GENE_ID=sub('(.*)=.*', '\\1', GENE_ID)) %>% 
  left_join(GENE_2_CLASS_DICT) 

VIR_sum <- 
  VIR %>% 
  group_by(genome) %>% 
  arrange(class) %>% 
  summarise(all_vir_genes=paste(unique(GENE_ID), collapse = '_'), 
            vir_types = paste(unique(type[!is.na(type)]), collapse='_'), 
            vir_subtypes=paste(unique(subtype[!is.na(subtype)]), collapse='_'), 
            vir_classes=paste(unique(class[!is.na(class)]), collapse='_'), 
            vir_subclasses=paste(unique(subclass[!is.na(subclass)]), collapse='_'))

STRESS <- 
  AMR_finder %>%
  transmute(genome=genome, 
            GENE_ID=stress_ID) %>% 
  separate_rows(GENE_ID, sep = ',') %>% 
  filter(!grepl('PARTIAL', GENE_ID)) %>% 
  mutate(GENE_ID=sub('(.*)=.*', '\\1', GENE_ID)) %>% 
  left_join(GENE_2_CLASS_DICT) %>% 
  mutate(class=ifelse(is.na(class), subtype, class))


### All genomes have copper and gold resisitance
# maybe quantify copper resistance with # copper res genes?
STRESS %>% filter(class == 'COPPER') %>% group_by(GENE_ID) %>% tally()
STRESS %>% filter(class == 'GOLD') %>% group_by(GENE_ID) %>% tally()

STRESS_sum <- 
  STRESS %>% 
  separate_rows(class, sep = '/') %>% 
  group_by(genome) %>% 
  arrange(class) %>% 
  summarise(all_stress_genes=paste(unique(GENE_ID), collapse = '_'), 
            stress_classes=paste(unique(class[!is.na(class)]), collapse='_'), 
            num_stress_classes = length(unique(class[!is.na(class)]))) %>% 
  arrange(desc(num_stress_classes))


STRESS <- master_meta %>%
  dplyr::select(genome, tree_clade ,Source, leiden_cluster, year, country) %>%
  right_join(STRESS)  %>% 
  separate_rows(class, sep = '/') %>%
  unique()



master_meta <- 
  master_meta %>%
  left_join(STRESS_sum) %>% 
  left_join(AMR_sum) %>% 
  filter(tree_clade %in% c('clade_one', 'clade_two', 'clade_three')) %>% 
  mutate(resist_type=paste(AMR_classes, stress_classes, sep = '_'), 
         tree_clade=factor(tree_clade, levels=c('clade_one', 'clade_two', 'clade_three')), 
         AMR_classes = factor(AMR_classes), 
         stress_classes = factor(stress_classes))




### ADD MLST DATA ###
MLST <- 
  read_tsv('output/MLST_JAN.tsv',
           col_names = c('ID', 'species','MLST',
                         'aroC','dnaN','hemD',
                         'hisD','purE','sucA','thrA')) %>% 
  mutate(genome=sub('.fna.gz','',ID)) %>% 
  dplyr::select(genome, MLST)




### maybe clean this up a little
master_meta <- 
  master_meta %>%
  left_join(MLST) %>% 
  write_tsv('./output/06_master_meta.tsv')


master_meta %>% group_by(tree_clade, MLST) %>% tally() %>% arrange(desc(n))

# overview, num resistances per clade

# by clade by host
  # percent resistant to individual AMR classes 
  # percent resistant to individual stress classes

# # fig 2
# P5_AMR_classes_violin <- master_meta %>%
#   # filter(tree_clade %in% c('clade_one', 'clade_two', 'clade_three')) %>% 
#   ggplot(aes(x=tree_clade, y=num_AMR_classes)) +
#   geom_violin() +
#   ggtitle('number AMR classes per genome')
# 
# P5_AMR_classes_violin

# ggsave('P5.jpg', width = 9, height = 5, units = 'in')


# # tree clade AVG Num AMR over time
# P6_AMR_year <-
#   master_meta %>% group_by(tree_clade, year) %>% 
#   summarise(avg_num_AMR_classes=mean(num_AMR_classes), 
#             se_num_AMR = sd(num_AMR_classes)/sqrt(n())) %>%
#   ggplot(aes(x=year, y=avg_num_AMR_classes, color=tree_clade)) + 
#   geom_line() + 
#   ggtitle('number AMR classes over time')
# 
# P6_AMR_year
# # ggsave('P6.jpg', width = 9, height = 5, units = 'in')



# tree clade AVG num Stress over time
# P7_stress_year <- master_meta %>%
#   group_by(tree_clade, year) %>% 
#   summarise(avg_num_stress_classes=mean(num_stress_classes), 
#             se_num_AMR = sd(num_stress_classes)/sqrt(n())) %>%
#   ggplot(aes(x=year, y=avg_num_stress_classes, color=tree_clade)) + 
#   geom_line()+
#   ggtitle('number stress res classes over time')
# 
# P7_stress_year
# 
# # ggsave('P7.jpg', width = 9, height = 5, units = 'in')
# 
# 
# 
# 
# P8_stress_violin <- 
#   master_meta %>%
#   # filter(tree_clade %in% c('clade_one', 'clade_two', 'clade_three')) %>% 
#   ggplot(aes(x=tree_clade, y=num_stress_classes)) +
#   geom_violin() + 
#   ggtitle('number stress classes over time')
# 
# P8_stress_violin
# 
# # ggsave('P8.jpg', width = 9, height = 5, units = 'in')
# 
# 
# P9_both_violin <- 
#   master_meta %>%
#   filter(tree_clade %in% c('clade_one', 'clade_two', 'clade_three')) %>% 
#   ggplot(aes(x=tree_clade, y=num_stress_classes + num_AMR_classes)) +
#   geom_violin() + ggtitle('number AMR + Stress classes')
# 
# P9_both_violin
# 
### MOST COMMON AMR GENOTYPES BY TREE CLADE ###
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(12)

mycolors <- c(mycolors, 'white')

# P10_AMR_genotypes <-
#   master_meta %>%
#   count(tree_clade, AMR_classes) %>%
#   left_join(tree_clade_totals) %>% 
#   mutate(ptot=n/tot_n_genomes *100) %>% 
#   arrange(desc(ptot)) %>% 
#   write_tsv('./output/AMR_class_genotypes_by_clade.tsv') %>% 
#   mutate(classes=fct_lump(AMR_classes, n = 12, w = ptot))  %>% 
#   ggplot(aes(x=tree_clade, y=ptot, fill=classes)) +
#   geom_col(color='black') +
#   scale_fill_manual(values=mycolors) + 
#   theme(legend.text = element_text(size=6)) + 
#   ggtitle('AMR genotypes by clade')
# 
# P10_AMR_genotypes  
# 
# ggsave('P10.jpg', width = 12, height = 5, units = 'in')


# MOST COMMON STRESS GENOTYPES BY TREE CLADE
# P11_stress_genotypes <- 
#   master_meta %>%
#   count(tree_clade, stress_classes) %>%
#   left_join(tree_clade_totals) %>% 
#   mutate(ptot=n/tot_n_genomes *100) %>% 
#   arrange(desc(ptot)) %>% 
#   write_tsv('./output/stress_class_genotypes_by_clade.tsv') %>% 
#   mutate(classes=fct_lump(stress_classes, n = 12, w = ptot))  %>% 
#   ggplot(aes(x=tree_clade, y=ptot, fill=classes)) +
#   geom_col(color='black') + 
#   scale_fill_manual(values=mycolors)+
#   theme(legend.text = element_text(size=6)) + 
#   ggtitle('stress resistance genotypes by clade')
# P11_stress_genotypes

# ggsave('P11.jpg', width = 12, height = 5, units = 'in')


 #### These are good #####
# percent resist to AMR for isolation host by clade
Percent_resist_ISO_clade <- 
  AMR %>%
  dplyr::select(genome, Source, class, tree_clade) %>%
  unique() %>%
  mutate(class=factor(class), 
         tree_clade=factor(tree_clade)) %>% 
  count(class, Source, tree_clade, .drop=F) %>%
  left_join(tree_clade_iso_totals) %>% 
  filter(class != 'EFFLUX') %>%
  mutate(prop_genomes=n/tot_n_genomes, 
         perc_genomes=prop_genomes *100) 

Percent_resist_ISO_clade %>% 
  arrange(tree_clade) %>% 
  write_tsv('./output/Percent_AMR_resist_isolation_clade.tsv')


mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(16)


names(mycolors) <- unique(Percent_resist_ISO_clade$class)

P12_percent_resist_ISO_clade_one <-
  Percent_resist_ISO_clade %>% 
  filter(tree_clade == 'clade_one') %>%
  ggplot(aes(x=class, y=perc_genomes, fill=class)) +
  geom_col() +
  geom_text(aes(x=class, y=50, label=round(perc_genomes, 1))) +
  scale_fill_manual(values = mycolors) +
  facet_wrap(~Source+ tot_n_genomes, nrow=1) +
  scale_y_continuous(limits = c(0,100),breaks = c(0,50,100))+  coord_flip()  +
  theme(axis.text.x=element_text(size=8), 
        legend.position = 'none') + 
  ylab('Percent of Genomes') + 
  xlab('')

P12_percent_resist_ISO_clade_one

ggsave('figure7A.jpg', width = 12, height = 5, units = 'in', bg='white')


P13_percent_resist_ISO_clade_two <- 
  Percent_resist_ISO_clade %>%
  filter(tree_clade == 'clade_two') %>%
  ggplot(aes(x=class, y=perc_genomes, fill=class)) +
  geom_col() +
  geom_text(aes(x=class, y=50, label=round(perc_genomes, 1))) +
  scale_fill_manual(values = mycolors) +
  facet_wrap(~Source+ tot_n_genomes, nrow=1) +
  scale_y_continuous(limits = c(0,100),breaks = c(0,50,100))+
  coord_flip()  +
  theme(axis.text.x=element_text(size=8), 
        legend.position = 'none') + 
  ylab('Percent of Genomes') + 
  xlab('')
P13_percent_resist_ISO_clade_two

ggsave('figure7B.jpg', width = 12, height = 5, units = 'in', bg='white')



P14_percent_resist_ISO_clade_three <- Percent_resist_ISO_clade %>%
  filter(tree_clade == 'clade_three') %>%
  ggplot(aes(x=class, y=perc_genomes, fill=class)) +
  geom_col() +
  geom_text(aes(x=class, y=50, label=round(perc_genomes, 1))) +
  scale_fill_manual(values = mycolors) +
  facet_wrap(~Source+ tot_n_genomes, nrow=1) +
  scale_y_continuous(limits = c(0,100),breaks = c(0,50,100))+  
  coord_flip()  +
  theme(axis.text.x=element_text(size=8), 
        legend.position = 'none') + 
  ylab('Percent of Genomes') + 
  xlab('')

P14_percent_resist_ISO_clade_three
ggsave('figure7C.jpg', width = 12, height = 5, units = 'in', bg='white')

#

##### DO ANY OF THESE CLASSES OCCUR AT HIGHER THAN EXPECTED RATES
##### ASSUMING THERE IS NO DEPENDENCE BETWEEN ISO AND CLASS
#### COLISTIN INVESTIGATION ######

country_iso_tree_tots <- 
  master_meta %>%
  count(country, ISOLATION, tree_clade, name='tot_n_genomes', .drop=F)

Percent_resist_ISO_clade_country_colistin <- 
  AMR %>%
  # filter(country %in% c('USA', 'United Kingdom')) %>% 
  dplyr::select(genome,country, Source, class, tree_clade) %>%
  unique() %>%
  mutate(Source=factor(Source), 
         class=factor(class), 
         tree_clade=factor(tree_clade), 
         country=factor(country)) %>% 
  count(class, country, Source, tree_clade, .drop=F, name='number_resistant') %>%
  left_join(country_iso_tree_tots) %>% 
  mutate(prop_genomes_resistant=number_resistant/tot_n_genomes) %>%
  filter(class=='COLISTIN') %>% 
  filter(!is.na(prop_genomes_resistant)) 

master_meta %>% filter(tree_clade == 'clade_two') %>% count(country, ISOLATION) %>% filter(ISOLATION == 'Swine')
master_meta %>% filter(tree_clade == 'clade_one') %>% count(country, ISOLATION) %>% filter(ISOLATION == 'Swine')



P15_percent_resist_colistin_swine <- 
  Percent_resist_ISO_clade_country_colistin %>%
  # filter(tree_clade == 'clade_one') %>% 
  filter(Source == 'Swine') %>% 
  filter(prop_genomes_resistant > 0) %>% 
  ggplot(aes(x=country, y=prop_genomes_resistant, fill=tree_clade)) +
  geom_col() + 
  geom_text(aes(label=tot_n_genomes)) +
  ggtitle('proportion of swine associated isolates with predicted resistance to colistin', 
          'countries with no detected colistin resistance in swine associated isolates not shown')
P15_percent_resist_colistin_swine
ggsave('P14.jpg', width = 7, height = 5, units = 'in')


######


# Percent_resist_ISO_clade %>%
#   filter(ISOLATION == 'swine') %>%
#   ggplot(aes(x=class, y=prop_genomes, fill=class)) +
#   geom_col() +
#   geom_text(aes(x=class, y=.5, label=round(prop_genomes, 3) * 100)) +
#   scale_fill_manual(values=mycolors) +
#   facet_wrap(~ISOLATION+tree_clade+tot_n_genomes ) +
#   coord_flip()
# 
# 
# Percent_resist_ISO_clade %>%
#   filter(ISOLATION == 'turkey') %>%
#   ggplot(aes(x=class, y=prop_genomes, fill=class)) +
#   geom_col() +
#   geom_text(aes(x=class, y=.5, label=round(prop_genomes, 3) * 100)) +
#   scale_fill_manual(values=mycolors) +
#   facet_wrap(~ISOLATION+tree_clade + tot_n_genomes) +
#   coord_flip()
# 
# 
# Percent_resist_ISO_clade %>%
#   filter(ISOLATION == 'chicken') %>%
#   ggplot(aes(x=class, y=prop_genomes, fill=class)) +
#   geom_col() +
#   geom_text(aes(x=class, y=.5, label=round(prop_genomes, 3) * 100)) +
#   scale_fill_manual(values=mycolors) +
#   facet_wrap(~ISOLATION+tree_clade + tot_n_genomes) +
#   coord_flip()
# 
# master_meta %>% group_by(ISOLATION) %>% tally()
# 
# 
# Percent_resist_ISO_clade %>%
#   filter(ISOLATION == 'human') %>%
#   ggplot(aes(x=class, y=prop_genomes, fill=class)) +
#   geom_col() +
#   geom_text(aes(x=class, y=.5, label=round(prop_genomes, 3) * 100)) +
#   scale_fill_manual(values=mycolors) +
#   facet_wrap(~ISOLATION+tree_clade + tot_n_genomes) +
#   coord_flip()


######## Percent of genomes with STRESS genotypes #######
### USE THESE TOO
Percent_stressRes_ISO_clade <- 
  STRESS %>%
  dplyr::select(genome, Source, class, tree_clade) %>%
  unique() %>%
  mutate(Source=factor(Source), 
         class=factor(class), 
         tree_clade=factor(tree_clade)) %>% 
  count(class, Source, tree_clade, .drop=F) %>%
  left_join(tree_clade_iso_totals) %>% 
  filter(class != 'EFFLUX') %>%
  mutate(prop_genomes=n/tot_n_genomes) %>% 
  arrange(tree_clade)


Percent_stressRes_ISO_clade %>% 
  write_tsv('./output/Percent_stress_resist_isolation_clade.tsv')


num_cols <- length(levels(Percent_stressRes_ISO_clade$class))

mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(num_cols)

names(mycolors) <-levels(Percent_stressRes_ISO_clade$class)




P16_percent_stress_ISO_clade_one <- 
  Percent_stressRes_ISO_clade %>% 
  filter(tree_clade == 'clade_one') %>%
  ggplot(aes(x=class, y=prop_genomes, fill=class)) +
  geom_col() +
  geom_text(aes(x=class, y=.5, label=round(prop_genomes, 3) * 100)) +
  scale_fill_manual(values = mycolors) +
  facet_wrap(~Source+tree_clade + tot_n_genomes, nrow=1) +
  coord_flip() +
  theme(axis.text.x=element_text(angle=-45, hjust = .1))

P16_percent_stress_ISO_clade_one

ggsave('P16.jpg', width = 12, height = 5, units = 'in')




P17_percent_stress_ISO_clade_two <- 
  Percent_stressRes_ISO_clade %>%
  filter(tree_clade == 'clade_two') %>%
  ggplot(aes(x=class, y=prop_genomes, fill=class)) +
  geom_col() +
  geom_text(aes(x=class, y=.5, label=round(prop_genomes, 3) * 100)) +
  scale_fill_manual(values=mycolors) +
  facet_wrap(~Source+tree_clade +tot_n_genomes, nrow=1) +
  coord_flip() + 
  theme(axis.text.x=element_text(angle=-45, hjust = .1))

P17_percent_stress_ISO_clade_two

ggsave('P16.jpg', width = 12, height = 5, units = 'in')


P18_percent_stress_ISO_clade_three <- 
  Percent_stressRes_ISO_clade %>%
  filter(tree_clade == 'clade_three') %>%
  ggplot(aes(x=class, y=prop_genomes, fill=class)) +
  geom_col() +
  geom_text(aes(x=class, y=.5, label=round(prop_genomes, 3) * 100)) +
  scale_fill_manual(values=mycolors) +
  facet_wrap(~Source+tree_clade + tot_n_genomes, nrow=1) +
  coord_flip() + 
  ylim(0,1) + 
  theme(axis.text.x=element_text(angle=-45, hjust = .1))

P18_percent_stress_ISO_clade_three
ggsave('P16.jpg', width = 12, height = 5, units = 'in')


###### num copp res genes
# USE THIS

P19_num_coppres_genes <- 
  STRESS %>%
  count(Source, genome, class, tree_clade) %>% 
  group_by(tree_clade, Source, class) %>% 
  summarise(mean_num_res=mean(n)) %>% 
  filter(class == 'COPPER') %>% 
  ggplot(aes(x=Source, y=mean_num_res, fill=class)) +
  geom_col() + facet_wrap(~tree_clade) + 
  ggtitle('average number of genes predicted to confer copper resistance')

P19_num_coppres_genes

##############
#############
#########
#####
#

# Closest REFs

# 
# master_meta %>%
#   group_by(tree_clade, closest_ref) %>%
#   tally() %>% 
#   filter(!is.na(tree_clade)) %>% 
#   ggplot(aes(x=closest_ref, y=n, fill=tree_clade)) + 
#   geom_col()
# 
# 
# master_meta %>%
#   group_by(tree_clade, closest_ref) %>%
#   tally() %>% 
#   filter(!is.na(tree_clade)) %>% 
#   ggplot(aes(x=tree_clade, y=n, fill=closest_ref)) + 
#   geom_col()

####### THIS ONE
P20_closest_reference_continent_clade <- 
  master_meta %>%
  group_by(continent, tree_clade, closest_ref) %>%
  tally() %>% 
  filter(!is.na(tree_clade)) %>% 
  ggplot(aes(x=continent, y=n, fill=closest_ref)) + 
  geom_col() + facet_wrap(~tree_clade) + 
  theme(axis.text.x = element_text(angle=-45, hjust = .1))
P20_closest_reference_continent_clade

master_meta %>%
  group_by(year, continent, tree_clade, closest_ref) %>%
  tally() %>% 
  filter(tree_clade == 'clade_one') %>% 
  ggplot(aes(x=year, y=n, color=closest_ref)) + 
  geom_line() + facet_wrap(~continent, scales = 'free') + 
  theme(axis.text.x = element_text(angle=-45, hjust = .1))+ 
  ggtitle('Clade_one closest reference genomes over time')

#####


master_meta %>%
  group_by(year, continent, tree_clade, closest_ref) %>%
  tally() %>%
  filter(continent == 'North America') %>% 
  filter(tree_clade == 'clade_one') %>% 
  ggplot(aes(x=year, y=n, color=closest_ref)) + 
  geom_line() + facet_wrap(~continent, scales = 'free') + 
  theme(axis.text.x = element_text(angle=-45, hjust = .1))+ 
  ggtitle('Clade_one closest reference genomes over time')

# # rando figs # # SHOULD HAVE FIGURES SCRIPT SEPARATE FROM CLEAN SCRIPT
# 
# master_meta %>%
#   mutate(clust=fct_lump(as.factor(Isolation), 7)) %>%
#   ggplot(aes(x=T1, y=T2, color=clust)) +
#   geom_point(alpha=.5) +
#   scale_color_brewer(palette = 'Set2') + 
#   geom_point(data=filter(master_meta, tree_rep == TRUE), color='black', alpha=.5) + 
#   ggtitle('TSNE showing host species', 'black dots are tree representatives')
# 
# 
# 
# master_meta %>% 
#   mutate(clust=fct_lump(as.factor(Isolation), 7)) %>%
#   filter(!clust %in% c('human', 'other')) %>% 
#   ggplot(aes(x=T1, y=T2, color=clust)) +
#   geom_point(alpha=.75) +
#   scale_color_brewer(palette = 'Set1') + 
#   ggtitle('TSNE, limited to only ag hosts')
# 
# 
# master_meta %>% 
#   mutate(clust=fct_lump(as.factor(state), 7)) %>%
#   filter(!clust %in% c('unknown', 'other', 'Other')) %>% 
#   ggplot(aes(x=T1, y=T2, color=clust)) +
#   geom_point(alpha=.75) +
#   scale_color_brewer(palette = 'Set1') + 
#   ggtitle('tsne, showing states of isolation')
# 
# 
# master_meta %>%
#   ggplot(aes(x=T1, y=T2, color=factor(leid_clust))) +
#   geom_point(alpha=.2) + scale_color_brewer(palette = 'Set2') + 
#   ggtitle('TSNE showing accessory genome cluster')

# master_meta %>% ggplot(aes(x=U1, y=U2, color=country)) + geom_point(alpha=.2)
# master_meta %>% ggplot(aes(x=T1, y=U2, color=country)) + geom_point(alpha=.2)
# 
# master_meta %>% ggplot(aes(x=T2, y=U2, color=country)) + geom_point(alpha=.2)
# master_meta %>% ggplot(aes(x=T1, y=U1, color=country)) + geom_point(alpha=.2)
# 
# master_meta %>% ggplot(aes(x=T1+U1, y=T2+U2, color=country)) + geom_point(alpha=.2)

# 
# top_countries <- master_meta %>% group_by(country) %>% tally() %>% top_n(n = 5, wt = n) %>%
#   select(country) %>% unlist()
# 
# master_meta <- master_meta %>% 
#   mutate(Country=case_when(
#     country %in% top_countries ~ country, 
#     TRUE                       ~ 'other'), 
#     Country = factor(Country, levels = c('USA', 'United Kingdom', 'Denmark', 'Germany', 'unknown', 'other')))
# 
# 
# # TSNE by country
# master_meta %>% ggplot(aes(x=T1, y=T2, color=Country)) + geom_point(alpha=.3)+ 
#   guides(color = guide_legend(override.aes = list(alpha=1)))  + 
#   ggtitle('t-SNE visualization of genome similarities', '5634 total genomes') + scale_color_brewer(palette = 'Set1')#+
# # geom_point(data = ANNOTATE, aes(x=T1, y=T2), show.legend = FALSE, color='black', size=3, shape=24, fill='white', stroke=1)
# 
# 
# master_meta %>% ggplot(aes(x=U1, y=U2, color=Country)) + geom_point(alpha=.3)+ 
#   guides(color = guide_legend(override.aes = list(alpha=1)))  + 
#   ggtitle('umap visualization of genome similarities', '5634 total genomes') + scale_color_brewer(palette = 'Set1')#+


# TSNE by SGI4
# master_meta %>% ggplot(aes(x=T1, y=T2, color=SGI4_class)) + geom_point(alpha=.3)+ 
#   guides(color = guide_legend(override.aes = list(alpha=1)))  +
#   ggtitle('t-SNE visualization of genome similarities', '5634 total genomes') +
#   scale_color_brewer(palette = 'Set1')+ 
#   ggtitle("TSNE showing SGI4")
# geom_point(data = ANNOTATE, aes(x=T1, y=T2), show.legend = FALSE, color='black', size=3, shape=24, fill='white', stroke=1)
# master_meta %>% ggplot(aes(x=U1, y=U2, color=SGI4_class)) + geom_point(alpha=.3)+ 
#   guides(color = guide_legend(override.aes = list(alpha=1)))  +
#   ggtitle('t-SNE visualization of genome similarities', '5634 total genomes') +
#   scale_color_brewer(palette = 'Set1')#+
# geom_point(data = ANNOTATE, aes(x=T1, y=T2), show.legend = FALSE, color='black', size=3, shape=24, fill='white', stroke=1)

# master_meta %>% ggplot(aes(x=T1+U1, y=T2+U2, color=SGI4_class)) + geom_point(alpha=.3)+ 
#   guides(color = guide_legend(override.aes = list(alpha=1)))  +
#   ggtitle('t-SNE visualization of genome similarities', '5634 total genomes') +
#   scale_color_brewer(palette = 'Set1')#+
# geom_point(data = ANNOTATE, aes(x=T1, y=T2), show.legend = FALSE, color='black', size=3, shape=24, fill='white', stroke=1)



# 
# # TSNE by MDR
# master_meta %>% ggplot(aes(x=T1, y=T2, color=MDR_class)) + geom_point(alpha=.3)+ 
#   guides(color = guide_legend(override.aes = list(alpha=1)))  + 
#   ggtitle('TSNE showing MDR module') + scale_color_brewer(palette = 'Set1')# + 
# # geom_text(data = ANNOTATE, aes(x=T1, y=T2, label=lab), show.legend = FALSE, color='black') + 
# # geom_point(data = ANNOTATE, aes(x=T1, y=T2), show.legend = FALSE, color='black', size=3, shape=24, fill='white', stroke=1)


# TSNE by MDR
# master_meta %>% ggplot(aes(x=U1, y=U2, color=MDR_class)) + geom_point(alpha=.3)+ 
#   guides(color = guide_legend(override.aes = list(alpha=1)))  + 
#   ggtitle('t-SNE visualization of genome similarities', '5634 total genomes') + scale_color_brewer(palette = 'Set1')# + 
# geom_text(data = ANNOTATE, aes(x=T1, y=T2, label=lab), show.legend = FALSE, color='black') + 
# geom_point(data = ANNOTATE, aes(x=T1, y=T2), show.legend = FALSE, color='black', size=3, shape=24, fill='white', stroke=1)


# grep('USDA',master_meta$ID, value = TRUE)
# ANNOTATE <- master_meta[master_meta$ID == 'USDA15WA1.fna',]
# 
# ANNOTATE$lab <- 'USDAWA15-1'
# 
# ANNOTATE$T1
# overall SGI4
# combine these two

SGI4_tot <- 
  master_meta %>%
  group_by(targ_presence_SGI4, tree_clade) %>%
  tally() %>%
  ungroup() %>% 
  mutate(class=targ_presence_SGI4 , 
         type='SGI4') %>% 
  dplyr::select(class, n, type, tree_clade)

MDR_tot <- master_meta %>%
  group_by(targ_presence_MDR, tree_clade) %>%
  tally() %>% 
  ungroup() %>% 
  mutate(class=targ_presence_MDR, 
         type='MDR') %>% 
  dplyr::select(class, n, type, tree_clade)

thrW18kb_tot <- master_meta %>%
  group_by(targ_presence_thrW18kb_GI, tree_clade) %>%
  tally() %>% 
  ungroup() %>% 
  mutate(class=targ_presence_thrW18kb_GI, 
         type='thrW18kb') %>% 
  dplyr::select(class, n, type, tree_clade)


sopE_mTmV_tot <- master_meta %>%
  group_by(targ_presence_SO4698_sopE_mTmV, tree_clade) %>%
  tally() %>% 
  ungroup() %>% 
  mutate(class=targ_presence_SO4698_sopE_mTmV, 
         type='sopE_mTmV') %>% 
  dplyr::select(class, n, type, tree_clade)










all_tot <- rbind(MDR_tot, SGI4_tot, thrW18kb_tot, sopE_mTmV_tot) %>%
  # mutate(type=factor(type, levels=c('SGI4', 'MDR'))) %>% 
  filter(!is.na(class))

P21_genomic_feats_by_tree_clade <- all_tot %>%
  filter(tree_clade != 'other') %>% 
  ggplot(aes(x=0, y=n, fill=class)) +
  geom_col(color='black')+
  facet_wrap(~type) +
  xlim(-.75,.75) + 
  # theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.length.x = unit(0, units = 'mm'))  + 
  ggtitle('Overall presence of Interesting genomic features')+
  scale_fill_brewer(palette = 'Set1') + ylab('Number of genomes') + 
  facet_wrap(~tree_clade + type, ncol = 4, scales = 'free')

P21_genomic_feats_by_tree_clade
# ggsave(P_genomic_feats_by_tree_clade, filename = 'output/genomic_feature_frequency_by_clade.jpg',width = 7, height = 5, units = 'in')


# Mosaic plot here #


# master_meta %>%
#   group_by(SGI4_class) %>%
#   tally() %>%
#   ggplot(aes(x=0, y=n, fill=SGI4_class)) +
#   geom_col(color='black')+
#   xlim(-1,1) + 
#   theme_bw() +
#   theme(axis.title.x = element_blank(), 
#         axis.text.x = element_blank())  + 
#   ggtitle('Overall presence of SGI-4')+
#   scale_fill_brewer(palette = 'Set1')
# 
# master_meta %>%
#   group_by(MDR_class) %>%
#   tally() %>%
#   ggplot(aes(x=0, y=n, fill=MDR_class)) +
#   geom_col(color='black')+
#   xlim(-1,1) + 
#   theme_bw() +
#   theme(axis.title.x = element_blank(), 
#         axis.text.x = element_blank())  + 
#   ggtitle('Overall presence of MDR module')+
#   scale_fill_brewer(palette = 'Set1')
# 
# 



### start insert

# #### FIG SGI4 by country by tree_clade
# master_meta %>%
#   filter(!is.na(targ_presence_SGI4)) %>% 
#   filter(tree_clade != 'other') %>% 
#   group_by(country, targ_presence_SGI4, tree_clade) %>%
#   tally() %>%
#   ggplot(aes(x=country, y=n, fill=targ_presence_SGI4)) + geom_col(color='black') +
#   scale_fill_brewer(palette = 'Set1') +theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
#   ggtitle('Presence of SGI-4 in 4,5,12,i:- NCBI genomes') +
#   ylab('number of genomes') + 
#   facet_wrap(~tree_clade, ncol=1)
# 
# #### MDR by country by tree clade
# master_meta %>%
#   filter(!is.na(targ_presence_MDR)) %>% 
#   filter(tree_clade != 'other') %>% 
#   group_by(country, targ_presence_MDR, tree_clade) %>%
#   tally()%>%
#   ggplot(aes(x=country, y=n, fill=targ_presence_MDR)) + geom_col(color='black') +
#   scale_fill_brewer(palette = 'Set1') + theme(axis.text.x = element_text(angle = -45, hjust = 0))+
#   ggtitle('Presence of MDR module in 4,5,12,i:- NCBI genomes') +
#   ylab('number of genomes') + 
#   facet_wrap(~tree_clade, ncol=1)
# 
# 
# master_meta %>%
#   filter(!is.na(targ_presence_thrW18kb_GI)) %>% 
#   filter(tree_clade != 'other') %>% 
#   group_by(country, targ_presence_thrW18kb_GI, tree_clade) %>%
#   tally()%>%
#   ggplot(aes(x=country, y=n, fill=targ_presence_thrW18kb_GI)) + geom_col(color='black') +
#   scale_fill_brewer(palette = 'Set1') + theme(axis.text.x = element_text(angle = -45, hjust = 0))+
#   ggtitle('Presence of MDR module in 4,5,12,i:- NCBI genomes') +
#   ylab('number of genomes') + 
#   facet_wrap(~tree_clade, ncol=1)
# 
# 
# master_meta %>%
#   filter(!is.na(targ_presence_SO4698_sopE_mTmV)) %>% 
#   filter(tree_clade != 'other') %>% 
#   group_by(country, targ_presence_SO4698_sopE_mTmV, tree_clade) %>%
#   tally()%>%
#   ggplot(aes(x=country, y=n, fill=targ_presence_SO4698_sopE_mTmV)) + geom_col(color='black') +
#   scale_fill_brewer(palette = 'Set1') + theme(axis.text.x = element_text(angle = -45, hjust = 0))+
#   ggtitle('Presence of MDR module in 4,5,12,i:- NCBI genomes') +
#   ylab('number of genomes') + 
#   facet_wrap(~tree_clade, ncol=1)
# 
# 
# 



#### Metal island % genomes

# master_meta %>%
#   filter(!is.na(targ_presence_SGI4)) %>% 
#   group_by(country, targ_presence_SGI4, tree_clade) %>%
#   summarise(number_genomes=n()) %>%
#   group_by(country) %>%
#   mutate(tot_genomes_country=sum(number_genomes),
#          percent_tot_genomes=(number_genomes/tot_genomes_country)*100) %>%
#   ggplot(aes(x=country, y=percent_tot_genomes, fill=targ_presence_SGI4)) + 
#   geom_col(color='black') +
#   scale_fill_brewer(palette = 'Set1') +
#   theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
#   ggtitle('Presence of SGI-4 in 4,[5],12:i:- genomes') +
#   geom_label(aes(label=tot_genomes_country, y=40), fill='white') +
#   ylab('Percent of total genomes in each country') + 
#   facet_wrap(~tree_clade , ncol=1)


##### MDR % genomes

# 
# master_meta %>% 
#   filter(!is.na(targ_presence_MDR)) %>% 
#   group_by(country, targ_presence_MDR, tree_clade) %>%
#   summarise(number_genomes=n()) %>%
#   group_by(country) %>%
#   mutate(tot_genomes_country=sum(number_genomes),
#          percent_tot_genomes=(number_genomes/tot_genomes_country)*100) %>%
#   ggplot(aes(x=country, y=percent_tot_genomes, fill=targ_presence_MDR)) +
#   geom_col(color='black') +
#   scale_fill_brewer(palette = 'Set1') +
#   theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
#   ggtitle('Presence of MDR module in 4,[5],12:i:- genomes') +
#   geom_label(aes(label=tot_genomes_country, y=35), fill='white') + 
#   ylab('Percent of total genomes in each country') + 
#   facet_wrap(~tree_clade, ncol=1)



##### over time #####



# tmp <- master_meta %>% group_by(year) %>% tally()

# master_meta$year_num <- as.numeric(master_meta$year)

# master_meta %>%filter(year_num >2012) %>% 
#   ggplot(aes(x=T1, y=T2, color=year_num)) + geom_point(alpha=.3)+ 
#   # guides(color = guide_legend(override.aes = list(alpha=1)))  + theme_bw() + 
#   ggtitle('t-SNE visualization of genome similarities', '5634 total genomes') + 
#   scale_color_viridis_c() + theme_bw()
# 
# check <- master_meta %>% group_by(year) %>% tally()

# TEXT HERE 
# IT looks liek there is a relationship with the MDR module increasing
# over time, but its really just the south clade becomming more common

### USE ME ###
#MDR OVER TIME #
# master_meta %>%
#   group_by(year, targ_presence_MDR) %>%
#   summarise(number_genomes=n()) %>%
#   group_by(year) %>%
#   mutate(tot_genomes_country=sum(number_genomes),
#          percent_tot_genomes=(number_genomes/tot_genomes_country)*100) %>%
#   ungroup() %>% 
#   filter(year != 'unknown') %>%
#   # filter(as.numeric(year) >2005) %>% 
#   ggplot(aes(x=year, y=percent_tot_genomes, fill=targ_presence_MDR)) + geom_col() +
#   scale_fill_brewer(palette = 'Set1')+
#   theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
#   ggtitle('occurrence of MDR module by year') +
#   geom_text(aes(y=50, label=tot_genomes_country)) +
#   ylab('Percent total genomes')



#MDR OVER TIME BUT SPLIT BY TREE_CLADE
P22_MDR_time_clade <- master_meta %>% 
  filter(tree_clade != 'other') %>% 
  group_by(year, targ_presence_MDR, tree_clade) %>%
  summarise(number_genomes=n()) %>%
  group_by(year) %>% 
  mutate(tot_genomes_country=sum(number_genomes),
         percent_tot_genomes=(number_genomes/tot_genomes_country)*100) %>%
  ungroup() %>% 
  filter(year != 'unknown') %>%
  filter(as.numeric(year) >2005) %>% 
  ggplot(aes(x=year, y=percent_tot_genomes, fill=targ_presence_MDR)) + geom_col() +
  scale_fill_brewer(palette = 'Set1')+
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle('occurrence of MDR module by year') +
  geom_text(aes(y=50, label=tot_genomes_country)) +
  ylab('Percent total genomes (in each year)') + 
  facet_wrap(~tree_clade, ncol=1)

P22_MDR_time_clade


### USE ME ###
# SGI4 overtime
# master_meta %>%
#   group_by(year,targ_presence_SGI4) %>% 
#   summarise(number_genomes=n()) %>%
#   ungroup() %>% 
#   group_by(year) %>%
#   mutate(tot_genomes_country=sum(number_genomes),
#          percent_tot_genomes=(number_genomes/tot_genomes_country)*100) %>%
#   ungroup() %>% 
#   filter(year != 'unknown') %>%
#   filter(as.numeric(year) >2005) %>% 
#   ggplot(aes(x=year, y=percent_tot_genomes, fill=targ_presence_SGI4)) + geom_col() +
#   scale_fill_brewer(palette = 'Set1')+
#   theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
#   ggtitle('occurrence of SGI-4 by year') +
#   geom_text(aes(y=50, label=tot_genomes_country)) + 
#   ylab('Percent total genomes in each year')




# SGI4 OVER TIME SPLIT BY TREE_CLADE
P23_SGI4_time_clade <-
  master_meta %>%
  filter(tree_clade != 'other') %>% 
  group_by(year, targ_presence_SGI4, tree_clade) %>%
  summarise(number_genomes=n()) %>%
  ungroup() %>% 
  group_by(year) %>%
  mutate(tot_genomes_country=sum(number_genomes),
         percent_tot_genomes=(number_genomes/tot_genomes_country)*100) %>%
  ungroup() %>% 
  filter(year != 'unknown') %>%
  filter(as.numeric(year) >2005) %>% 
  ggplot(aes(x=year, y=percent_tot_genomes, fill=targ_presence_SGI4)) + geom_col() +
  scale_fill_brewer(palette = 'Set1')+
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle('occurrence of SGI-4 by year') +
  geom_text(aes(y=50, label=tot_genomes_country)) + 
  ylab('Percent total genomes in each year') + 
  facet_wrap(~tree_clade, ncol=1)

P23_SGI4_time_clade


# One genome in clade3 has full SGI4, also has partial MDR
# Canada, 2012, 
clade_3_SGI4 <- 
  master_meta %>% filter(tree_clade == 'clade_three') %>% 
  filter(targ_presence_SGI4 == 'full') %>% 
  dplyr::select(genome, country, year, ISOLATION, starts_with('targ_'))

clade_3_SGI4

#####

# thrW18kb OVER TIME SPLIT BY TREE_CLADE
P24_thrW18kb_GI_time_clade <-
  master_meta %>%
  filter(tree_clade != 'other') %>% 
  group_by(year, targ_presence_thrW18kb_GI, tree_clade) %>%
  summarise(number_genomes=n()) %>%
  ungroup() %>% 
  group_by(year) %>%
  mutate(tot_genomes_country=sum(number_genomes),
         percent_tot_genomes=(number_genomes/tot_genomes_country)*100) %>%
  ungroup() %>% 
  filter(year != 'unknown') %>%
  filter(as.numeric(year) >2005) %>% 
  ggplot(aes(x=year, y=percent_tot_genomes, fill=targ_presence_thrW18kb_GI)) + geom_col() +
  scale_fill_brewer(palette = 'Set1')+
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle('occurrence of thrW18kb_GI by year') +
  geom_text(aes(y=50, label=tot_genomes_country)) + 
  ylab('Percent total genomes in each year') + 
  facet_wrap(~tree_clade, ncol=1)

P24_thrW18kb_GI_time_clade


########
P25_SO4698_sopE_mTmV_time_clade <-
  master_meta %>%
  filter(tree_clade != 'other') %>% 
  group_by(year, targ_presence_SO4698_sopE_mTmV, tree_clade) %>%
  summarise(number_genomes=n()) %>%
  ungroup() %>% 
  group_by(year) %>%
  mutate(tot_genomes_country=sum(number_genomes),
         percent_tot_genomes=(number_genomes/tot_genomes_country)*100) %>%
  ungroup() %>% 
  filter(year != 'unknown') %>%
  filter(as.numeric(year) >2005) %>% 
  ggplot(aes(x=year, y=percent_tot_genomes, fill=targ_presence_SO4698_sopE_mTmV)) + geom_col() +
  scale_fill_brewer(palette = 'Set1')+
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle('occurrence of thrW18kb_GI by year') +
  geom_text(aes(y=50, label=tot_genomes_country)) + 
  ylab('Percent total genomes in each year') + 
  facet_wrap(~tree_clade, ncol=1)

P25_SO4698_sopE_mTmV_time_clade


SOPE_countries <- 
  master_meta %>% 
  filter(tree_clade =='clade_one') %>% 
  group_by(country, targ_presence_SO4698_sopE_mTmV) %>% tally() %>% 
  filter(targ_presence_SO4698_sopE_mTmV == 'full') %>% filter(n>2) %>% 
  pull(country)



SOPE_countries2 <- 
  master_meta %>% 
  filter(tree_clade =='clade_one') %>% 
  group_by(country, targ_presence_SO4698_sopE_mTmV) %>% tally() %>% 
  filter(targ_presence_SO4698_sopE_mTmV == 'full') %>% filter(n>50) %>% 
  pull(country)





P26_SO4698_sopE_mTmV_time_clade <-
  master_meta %>% 
  filter(country %in% SOPE_countries) %>% 
  filter(tree_clade == 'clade_one') %>% 
  group_by(year, targ_presence_SO4698_sopE_mTmV, country) %>%
  summarise(number_genomes=n()) %>%
  ungroup() %>% 
  group_by(year, country) %>%
  mutate(tot_genomes_country_year=sum(number_genomes),
         percent_tot_genomes=(number_genomes/tot_genomes_country_year)*100) %>%
  ungroup() %>% 
  filter(year != 'unknown') %>%
  filter(as.numeric(year) >2005) %>% 
  ggplot(aes(x=year, y=percent_tot_genomes, fill=targ_presence_SO4698_sopE_mTmV)) +
  geom_col() +
  scale_fill_brewer(palette = 'Set1')+
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle('occurrence of sopE_mTmV by year', 
          'only countries with at least 2 full sopE genomes') +
  # geom_text(aes(y=50, label=tot_genomes_country)) + 
  ylab('Percent total genomes in each year per country') + 
  facet_wrap(~country, ncol=11) + 
  theme(legend.position = 'bottom')

P26_SO4698_sopE_mTmV_time_clade


P27_SO4698_sopE_mTmV_time_clade2 <- 
  master_meta %>% 
  filter(country %in% SOPE_countries2) %>% 
  filter(tree_clade == 'clade_one') %>% 
  group_by(year, targ_presence_SO4698_sopE_mTmV, country) %>%
  summarise(number_genomes=n()) %>%
  ungroup() %>% 
  group_by(year, country) %>%
  mutate(tot_genomes_country_year=sum(number_genomes),
         percent_tot_genomes=(number_genomes/tot_genomes_country_year)*100) %>%
  ungroup() %>% 
  filter(year != 'unknown') %>%
  filter(as.numeric(year) >2005) %>% 
  ggplot(aes(x=year, y=percent_tot_genomes, fill=targ_presence_SO4698_sopE_mTmV)) +
  geom_col() +
  scale_fill_brewer(palette = 'Set1')+
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle('occurrence of sopE_mTmV by year', 
          'only countries with at least 50 sopE positive genomes') +
  # geom_text(aes(y=50, label=tot_genomes_country)) + 
  ylab('Percent total genomes in each year per country') + 
  facet_wrap(~country, ncol=20) + 
  theme(legend.position = 'bottom')

P27_SO4698_sopE_mTmV_time_clade2

######



######## TREE CLADES OVER TIME ####
# north clade is more abundant in earlier years
# once sequencing really gets going, south clade dominates
# patterns seen in other data might be a artifact of this change.
# looks like mdr and sgi4 are becoming more prevalent over time but it really is
# just that the strains carrying MDR and SGI4 already are becoming more abundant than
# those that are not.  Not a whole lot of evidence that SGI4 and MDR are being introduced
# to strains that are in the north clade.  Maybe a good experiment could be to 
# measure the transcongugation freq. moving SGI4 from south clade genome to north vs
# south to south.



# master_meta %>% 
#   filter(tree_clade != 'other') %>% 
#   group_by(year, tree_clade) %>% 
#   summarise(number_genomes=n()) %>%
#   ungroup() %>% 
#   group_by(year) %>% 
#   mutate(tot_genomes_country=sum(number_genomes),
#          percent_tot_genomes=(number_genomes/tot_genomes_country)*100) %>%
#   ungroup() %>% 
#   filter(year != 'unknown') %>%
#   # filter(as.numeric(year) >2005) %>% 
#   ggplot(aes(x=year, y=percent_tot_genomes, fill=tree_clade)) + geom_col() +
#   scale_fill_brewer(palette = 'Set1')+
#   theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
#   ggtitle('occurrence of genomes in tree_clade by year') +
#   geom_text(aes(y=50, label=tot_genomes_country)) + 
#   ylab('Percent total genomes (in each year)')


###NOW ONLY SOUTH CLADE SGI4 ###

P28_clade_one_SGI4 <- 
  master_meta %>%
  filter(tree_clade == 'clade_one') %>% 
  group_by(year, targ_presence_SGI4) %>%
  summarise(number_genomes=n()) %>%
  ungroup() %>% 
  group_by(year) %>% mutate(tot_genomes_country=sum(number_genomes),
                            percent_tot_genomes=(number_genomes/tot_genomes_country)*100) %>%
  ungroup() %>% 
  filter(year != 'unknown') %>%
  filter(as.numeric(year) >2005) %>% 
  ggplot(aes(x=year, y=percent_tot_genomes, fill=targ_presence_SGI4)) + geom_col() +
  scale_fill_brewer(palette = 'Set1')+
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle('occurrence of SGI-4 by year, clade one only') +
  geom_text(aes(y=50, label=tot_genomes_country)) + 
  ylab('Percent south clade genomes in each year')

P28_clade_one_SGI4

### NOW ONLY SOUTH CLADE MDR #####

P29_clade_one_MDR <- 
  master_meta %>%
  filter(tree_clade == 'clade_one') %>% 
  group_by(year, targ_presence_MDR) %>%
  summarise(number_genomes=n()) %>%
  group_by(year) %>%
  mutate(tot_genomes_country=sum(number_genomes),
         percent_tot_genomes=(number_genomes/tot_genomes_country)*100) %>%
  ungroup() %>% 
  filter(year != 'unknown') %>%
  filter(as.numeric(year) >2005) %>% 
  ggplot(aes(x=year, y=percent_tot_genomes, fill=targ_presence_MDR)) + geom_col() +
  scale_fill_brewer(palette = 'Set1')+
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle('occurrence of MDR module by year, South clade only') +
  geom_text(aes(y=50, label=tot_genomes_country)) +
  ylab('Percent south clade genomes in each year')



P29_clade_one_MDR

#########



# tmp <- master_meta %>% select(-c('#Organism Group',
#                                  'Assembly',
#                                  'K-mer group',
#                                  'Lat/Lon',
#                                  'Library Layout',
#                                  'Platform',
#                                  'Level',
#                                  'Method',
#                                  'PFGE Primary Enzyme Pattern',
#                                  'PFGE Secondary Enzyme Pattern',
#                                  'WGS Prefix',
#                                  'WGS Accession')) %>% 
#   select(ID, `Isolation Source`, Host, starts_with('iso'), everything())
# 

# need to combine 'Chicken' and 'chicken'
# some cases where Isolation is 'environmental' but there is useful host info
# Equine and Bovine need attention
# Meleagris gallopavo to turkey
# Poricne
# Bos taurus
# Equus caballus	
# Gallus gallus	

####### PASTE IN


# with human
P30_MDR <- 
  master_meta %>%
  filter(!is.na(targ_presence_MDR)) %>%
  # filter(tree_clade != 'other') %>% 
  # filter(!(Isolation %in% c('unknown', 'other', 'Food', 'environmental'))) %>%
  group_by(ISOLATION, targ_presence_MDR, tree_clade) %>%
  summarise(number_genomes=n()) %>%
  ungroup() %>% 
  ggplot(aes(x=ISOLATION, y=number_genomes, fill=targ_presence_MDR)) + geom_col() +
  scale_fill_brewer(palette = 'Set1')+
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle('occurrence of MDR module')+
  ylab('Number of genomes') + 
  facet_wrap(~tree_clade)

P30_MDR

# absolute no human
P31_MDR <- 
  master_meta %>%
  filter(tree_clade != 'other') %>% 
  filter(!is.na(targ_presence_MDR)) %>% 
  filter(!(ISOLATION %in% c('human', 'unknown', 'other', 'Food', 'environmental'))) %>%
  group_by(ISOLATION, targ_presence_MDR, tree_clade) %>%
  summarise(number_genomes=n()) %>%
  ungroup() %>% 
  ggplot(aes(x=ISOLATION, y=number_genomes, fill=targ_presence_MDR)) + geom_col() +
  scale_fill_brewer(palette = 'Set1')+
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle('occurrence of MDR module by host', 'excluding human isolates') +
  ylab('Number of genomes') + xlab('') + 
  facet_wrap(~tree_clade)

P31_MDR

# percent of total
P32_MDR <- 
  master_meta %>%
  filter(!is.na(targ_presence_MDR)) %>%
  filter(tree_clade != 'other') %>% 
  group_by(ISOLATION, targ_presence_MDR, tree_clade) %>%
  summarise(number_genomes=n()) %>%
  ungroup() %>% 
  group_by(ISOLATION) %>%
  mutate(tot_genomes_country=sum(number_genomes),
         percent_tot_genomes=(number_genomes/tot_genomes_country)*100) %>%
  ungroup() %>% 
  filter(!(ISOLATION %in% c('human', 'unknown', 'other', 'Food', 'environmental'))) %>%
  ggplot(aes(x=ISOLATION, y=percent_tot_genomes, fill=targ_presence_MDR)) + geom_col() +
  scale_fill_brewer(palette = 'Set1')+
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle('occurrence of MDR module by host', 'excluding human isolates') +
  geom_text(aes(y=50, label=tot_genomes_country)) + 
  ylab('Percent of genomes, in each host') + xlab('') + 
  facet_wrap(~tree_clade)

P32_MDR

P33_MDR <- master_meta %>%
  filter(!is.na(targ_presence_MDR)) %>%
  filter(tree_clade == 'clade_one') %>% 
  # filter(tree_clade == 'south') %>% 
  group_by(ISOLATION, targ_presence_MDR, tree_clade) %>%
  summarise(number_genomes=n()) %>%
  ungroup() %>% 
  group_by(ISOLATION) %>%
  mutate(tot_genomes_country=sum(number_genomes),
         percent_tot_genomes=(number_genomes/tot_genomes_country)*100) %>%
  ungroup() %>% 
  filter(!(ISOLATION %in% c('human', 'unknown', 'other', 'Food', 'environmental'))) %>%
  ggplot(aes(x=ISOLATION, y=percent_tot_genomes, fill=targ_presence_MDR)) + geom_col() +
  scale_fill_brewer(palette = 'Set1')+
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle('occurrence of MDR module by host', 'South clade only, excluding human isolates') +
  geom_text(aes(y=50, label=tot_genomes_country)) + 
  ylab('Percent of genomes, in each host') + xlab('') + 
  facet_wrap(~tree_clade)


P33_MDR
### SAME THING NOW FOR THE SGI4 ###


P34_SGI4 <- 
  master_meta %>%
  # filter(!(Isolation %in% c('unknown', 'other', 'Food', 'environmental'))) %>%
  group_by(ISOLATION, targ_presence_SGI4, tree_clade) %>%
  summarise(number_genomes=n()) %>%
  ungroup() %>% 
  ggplot(aes(x=ISOLATION, y=number_genomes, fill=targ_presence_SGI4)) + geom_col() +
  scale_fill_brewer(palette = 'Set1')+
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle('occurrence of SGI4 module')+
  ylab('Number of genomes') + 
  facet_wrap(~tree_clade)

P34_SGI4

# absolute no human
P35_SGI4 <- master_meta %>%
  filter(!is.na(targ_presence_SGI4)) %>% 
  filter(!(ISOLATION %in% c('human', 'unknown', 'other', 'Food', 'environmental'))) %>%
  group_by(ISOLATION, targ_presence_SGI4, tree_clade) %>%
  summarise(number_genomes=n()) %>%
  ungroup() %>% 
  ggplot(aes(x=ISOLATION, y=number_genomes, fill=targ_presence_SGI4)) + geom_col() +
  scale_fill_brewer(palette = 'Set1')+
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle('occurrence of SGI4 module by host', 'excluding human isolates') +
  ylab('Number of genomes') + xlab('') + 
  facet_wrap(~tree_clade)

P35_SGI4

# percent of total
P36_SGI4 <- 
  master_meta %>%
  filter(!is.na(targ_presence_SGI4)) %>%
  group_by(ISOLATION, targ_presence_SGI4, tree_clade) %>%
  summarise(number_genomes=n()) %>%
  ungroup() %>% 
  group_by(ISOLATION) %>%
  mutate(tot_genomes_country=sum(number_genomes),
         percent_tot_genomes=(number_genomes/tot_genomes_country)*100) %>%
  ungroup() %>% 
  filter(!(ISOLATION %in% c('human', 'unknown', 'other', 'Food', 'environmental'))) %>%
  ggplot(aes(x=ISOLATION, y=percent_tot_genomes, fill=targ_presence_SGI4)) + geom_col() +
  scale_fill_brewer(palette = 'Set1')+
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle('occurrence of SGI4 module by host', 'excluding human isolates') +
  geom_text(aes(y=50, label=tot_genomes_country)) + 
  ylab('Percent of genomes, in each host') + xlab('') + 
  facet_wrap(~tree_clade)


P36_SGI4


P37_SGI4 <- master_meta %>%
  filter(!is.na(targ_presence_SGI4)) %>%
  filter(tree_clade == 'clade_one') %>% 
  group_by(ISOLATION, targ_presence_SGI4, tree_clade) %>%
  summarise(number_genomes=n()) %>%
  ungroup() %>% 
  group_by(ISOLATION) %>%
  mutate(tot_genomes_country=sum(number_genomes),
         percent_tot_genomes=(number_genomes/tot_genomes_country)*100) %>%
  ungroup() %>% 
  filter(!(ISOLATION %in% c('human', 'unknown', 'other', 'Food', 'environmental'))) %>%
  ggplot(aes(x=ISOLATION, y=percent_tot_genomes, fill=targ_presence_SGI4)) + geom_col() +
  scale_fill_brewer(palette = 'Set1')+
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle('occurrence of SGI4 module by host', 'South clade only, excluding human isolates') +
  geom_text(aes(y=50, label=tot_genomes_country)) + 
  ylab('Percent of genomes, in each host') + xlab('') + 
  facet_wrap(~tree_clade)



P37_SGI4

##########

#### USE ME ####
# master_meta %>%
#   filter(!(Isolation %in% c('human', 'unknown', 'other', 'Food', 'environmental'))) %>%
#   filter(!is.na(SGI4_class)) %>% 
#   group_by(Isolation, SGI4_class, tree_clade) %>%
#   summarise(number_genomes=n()) %>%
#   ungroup() %>% 
#   ggplot(aes(x=Isolation, y=number_genomes, fill=SGI4_class)) + geom_col() +
#   scale_fill_brewer(palette = 'Set1')+
#   theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
#   ggtitle('occurrence of SGI4 module by host', 'excluding human isolates') +
#   xlab('') + 
#   ylab('number of genomes') + 
#   facet_wrap(~tree_clade)
# 
# 
# # percent of total
# master_meta %>% filter(!is.na(SGI4_class)) %>%  
#   group_by(Isolation, SGI4_class, tree_clade) %>%
#   summarise(number_genomes=n()) %>%
#   ungroup() %>% 
#   group_by(Isolation) %>%
#   mutate(tot_genomes_country=sum(number_genomes),
#          percent_tot_genomes=(number_genomes/tot_genomes_country)*100) %>%
#   ungroup() %>% 
#   filter(!(Isolation %in% c('human', 'unknown', 'other', 'Food', 'environmental'))) %>%
#   ggplot(aes(x=Isolation, y=percent_tot_genomes, fill=SGI4_class)) + geom_col() +
#   scale_fill_brewer(palette = 'Set1')+
#   theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
#   ggtitle('occurrence of SGI4 module byt host', 'all tree clades, excluding human isolates') +
#   geom_text(aes(y=50, label=tot_genomes_country))+
#   xlab('') + 
#   ylab('Percent total genomes')
# 
# ##
# 
# master_meta %>% filter(!is.na(SGI4_class)) %>%  
#   group_by(Isolation, SGI4_class, tree_clade) %>%
#   summarise(number_genomes=n()) %>%
#   ungroup() %>% 
#   group_by(Isolation) %>%
#   mutate(tot_genomes_country=sum(number_genomes),
#          percent_tot_genomes=(number_genomes/tot_genomes_country)*100) %>%
#   ungroup() %>% 
#   filter(!(Isolation %in% c('human', 'unknown', 'other', 'Food', 'environmental'))) %>%
#   ggplot(aes(x=Isolation, y=percent_tot_genomes, fill=SGI4_class)) + geom_col() +
#   scale_fill_brewer(palette = 'Set1')+
#   theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
#   ggtitle('occurrence of SGI4 module byt host', 'excluding human isolates') +
#   geom_text(aes(y=50, label=tot_genomes_country))+
#   xlab('') + 
#   ylab('Percent total genomes (in each host)') + facet_wrap(~tree_clade)
# 
# 
# master_meta %>% filter(!is.na(SGI4_class)) %>%  
#   filter(tree_clade == 'south') %>% 
#   group_by(Isolation, SGI4_class, tree_clade) %>%
#   summarise(number_genomes=n()) %>%
#   ungroup() %>% 
#   group_by(Isolation) %>%
#   mutate(tot_genomes_country=sum(number_genomes),
#          percent_tot_genomes=(number_genomes/tot_genomes_country)*100) %>%
#   ungroup() %>% 
#   filter(!(Isolation %in% c('human', 'unknown', 'other', 'Food', 'environmental'))) %>%
#   ggplot(aes(x=Isolation, y=percent_tot_genomes, fill=SGI4_class)) + geom_col() +
#   scale_fill_brewer(palette = 'Set1')+
#   theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
#   ggtitle('occurrence of SGI4 module byt host', 'South clade only excluding human isolates') +
#   geom_text(aes(y=50, label=tot_genomes_country))+
#   xlab('') + 
#   ylab('Percent total genomes (in each host)') + facet_wrap(~tree_clade)
# 
# 

# devtools::install_github('haleyjeppson/ggmosaic')
# 
# 
# library(ggmosaic)
# 
# # ggmosaic::geom_mosaic()
# 
# 
# sgi4_mosaic <- 
#   master_meta %>%
#   filter(!is.na(targ_presence_SGI4)) %>%  
#   group_by(ISOLATION, targ_presence_SGI4, tree_clade) %>%
#   summarise(number_genomes=n()) %>%
#   ungroup() #%>% 
#   # group_by(Isolation) %>%
#   # mutate(tot_genomes_country=sum(number_genomes),
#   #        percent_tot_genomes=(number_genomes/tot_genomes_country)*100) %>%
#   # ungroup() %>% 
#   # filter(!(Isolation %in% c('human', 'unknown', 'other', 'Food', 'environmental')))
# 
# 
# 
# sgi4_mosaic %>%
#   filter(ISOLATION %in% c('swine', 'cow', 'chicken', 'turkey')) %>%
#   mutate(targ_presence_SGI4=fct_relevel(targ_presence_SGI4, 'absent','partial','full'),
#          Isolation =fct_relevel(targ_presence_SGI4, 'chicken', 'turkey', 'swine', 'cow')) %>%
#   ggplot() +
#   geom_mosaic(aes(x=product(targ_presence_SGI4, ISOLATION), fill=targ_presence_SGI4, weight=number_genomes)) +
#   scale_fill_brewer(palette = 'Set1') +# theme_mosaic() +
#   theme(axis.text.x=element_text(angle=-45, hjust=-.21))+
#   facet_wrap(~tree_clade)
# 
# 


# sgi4_mosaic %>% filter(tree_clade == 'south') %>% 
#   filter(Isolation %in% c('swine', 'cow', 'chicken', 'turkey')) %>% 
#   mutate(SGI4_class=fct_relevel(SGI4_class, 'absent','partial','full'), 
#          Isolation =fct_relevel(Isolation, 'chicken', 'turkey', 'swine', 'cow')) %>% 
#   ggplot() +
#   geom_mosaic(aes(x=product(SGI4_class, Isolation), fill=SGI4_class, weight=number_genomes)) +
#   scale_fill_brewer(palette = 'Set1') + theme_mosaic() +
#   theme(axis.text.x=element_text(angle=-45, hjust=-.21))+
#   facet_wrap(~tree_clade)

# trying to add in tree_clade
# 
# sgi4_mosaic %>% 
#   filter(Isolation %in% c('swine', 'cow', 'chicken', 'turkey')) %>% 
#   mutate(SGI4_class=fct_relevel(SGI4_class, 'absent','partial','full'), 
#          Isolation =fct_relevel(Isolation, 'chicken', 'turkey', 'swine', 'cow')) %>% 
#   ggplot() +
#   geom_mosaic(aes(x=product(Isolation, tree_clade), fill=SGI4_class, weight=number_genomes)) +
#   scale_fill_brewer(palette = 'Set1') + theme_mosaic()




# 
# 
# mdr_mosaic <- 
#   master_meta %>% filter(!is.na(targ_presence_MDR)) %>%  
#   group_by(ISOLATION, targ_presence_MDR, tree_clade) %>%
#   summarise(number_genomes=n()) %>%
#   ungroup() 

# mdr_mosaic %>% 
#   filter(Isolation %in% c('swine', 'cow', 'chicken', 'turkey')) %>% 
#   mutate(MDR_class=fct_relevel(MDR_class, 'absent','partial','full'), 
#          Isolation =fct_relevel(Isolation, 'chicken', 'turkey', 'swine', 'cow')) %>% 
#   ggplot() +
#   geom_mosaic(aes(x=product(MDR_class, Isolation), fill=MDR_class, weight=number_genomes)) +
#   scale_fill_brewer(palette = 'Set1') + theme_mosaic() + facet_wrap(~tree_clade) +
#   theme(axis.text.x=element_text(angle=-45, hjust=-.21))
  

# 
# mdr_mosaic %>% filter(tree_clade == 'clade_one') %>% 
#   filter(ISOLATION %in% c('swine', 'cow', 'chicken', 'turkey')) %>% 
#   mutate(targ_presence_MDR=fct_relevel(targ_presence_MDR, 'absent','partial','full'), 
#          ISOLATION =fct_relevel(ISOLATION, 'chicken', 'turkey', 'swine', 'cow')) %>% 
#   ggplot() +
#   geom_mosaic(aes(x=product(targ_presence_MDR, ISOLATION), fill=targ_presence_MDR, weight=number_genomes)) +
#   scale_fill_brewer(palette = 'Set1') +
#   # theme_mosaic() +
#   facet_wrap(~tree_clade) + 
#   theme(axis.text.x=element_text(angle=-45, hjust=-.21)) 
# 
# 

# ### Coocurrance
# 
# 
# 
# ####
# 
# # 
# # test <- master_meta2 %>% group_by(SGI4_class, MDR_class) %>% tally() %>%
# #   spread(key=MDR_class, value=n)
# # 
# # test <- as.data.frame(test)
# # rownames(test) <- test$SGI4_class
# # test <- as.matrix(test[,-1])
# # 
# # image(test)
# # 
# # heatmap(test, Rowv = NA, Colv = NA, scale = 'none')
# # 
# # 
# # # install.packages("gplots")
# # library("gplots")
# # heatmap.2(test, scale = "none", col = bluered(100),
# #           trace = "none", density.info = "none")
# 
# 
# # master_meta$both_class <- paste(master_meta2$SGI4_class , master_meta2$MDR_class, sep=' / ')
# 
# # master_meta %>% group_by(both_class) %>% tally() %>% ggplot(aes(x=both_class, y=n)) + geom_col() +
# #   theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
# #   ggtitle('Co-occurrence of Metal and MDR islands',
# #           subtitle = 'SGI4 1st, MDR 2nd') + ylab('number of genomes')
# 
# 
# 
# ####### by country #######
# #
# # master_meta2$both_class2 <- ifelse(master_meta2$both_class == 'absent_fragment', 'other',
# #                                    ifelse(master_meta2$both_class == 'fragment_absent', 'other',
# #                                           ifelse(master_meta2$both_class =='fragment_fragment', 'other',
# #                                                  ifelse(master_meta2$both_class == 'fragment_full', 'other',
# #                                                         ifelse(master_meta2$both_class =='fragment_partial', 'other',
# #                                                                ifelse(master_meta2$both_class =='full_fragment', 'other',
# #                                                                       ifelse(master_meta2$both_class =='partial_absent', 'other',
# #                                                                              ifelse(master_meta2$both_class =='partial_fragment', 'other',master_meta2$both_class))))))))
# 
# 
# # master_meta2$both_class2 <- ifelse(master_meta2$both_class %in% c('full_full', 'absent_absent'), master_meta2$both_class, 'other')
# master_meta$both_class3 <- ifelse(master_meta$both_class %in% c('full / full', 'absent / absent', 'full / partial', 'partial / full'), master_meta$both_class, 'other')
# 
# master_meta$both_class3 <- factor(master_meta$both_class3, levels = c('absent / absent', 'full / full', 'full / partial', 'partial / full', 'other'))
# 
# 
# master_meta %>% group_by(both_class3) %>% tally() %>% ggplot(aes(x=both_class3, y=n)) + geom_col() +
#   theme_bw()+
#   theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
#   ggtitle('Co-occurrence of Metal and MDR islands',
#           subtitle = 'Metal 1st, MDR 2nd')
# 
# 
# 
# # this one for total 
# # master_meta %>% group_by(country, both_class3) %>%
# #   tally() %>% ggplot(aes(x=country, y=n, fill=both_class3)) + geom_col() +
# #   theme_bw()+scale_fill_brewer(palette = 'Set1')+
# #   theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
# #   labs(fill = "SGI4 / MDR") +
# #   ggtitle('Co-occurrence of Metal and MDR islands',
# #           subtitle = 'Metal 1st, MDR 2nd')
# # 
# # 
# ### THIS ONE!
# master_meta %>%
#   filter(tree_clade == 'south') %>% 
#   group_by(country, both_class) %>% summarise(number_genomes=n()) %>%
#   group_by(country) %>% mutate(tot_genomes_country=sum(number_genomes),
#                                percent_tot_genomes=(number_genomes/tot_genomes_country)*100) %>%
#   # filter(!(country %in% c('unknown', 'South Korea', 'Brazil', 'Mexico'))) %>%
#   ggplot(aes(x=country, y=percent_tot_genomes, fill=both_class)) + geom_col() +
#   scale_fill_brewer(palette = 'Set1')+
#   theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
#   ggtitle('Co-occurrence of SGI-4 and MDR islands, ONLY SOUTH',
#           subtitle = 'SGI-4 1st, MDR 2nd') +
#   labs(fill = "SGI4 / MDR") +
#   geom_text(aes(y=50, label=tot_genomes_country)) + 
#   ylab('Percent total genomes')
# 
# 
# 
# # coocur over time #
# # THIS ONE #
# # master_meta %>% 
# #   filter(tree_clade =='south') %>% 
# #   group_by(year, both_class) %>% summarise(number_genomes=n()) %>%
# #   group_by(year) %>% mutate(tot_genomes_country=sum(number_genomes),
# #                             percent_tot_genomes=(number_genomes/tot_genomes_country)*100) %>%
# #   filter(as.numeric(year) > 2005) %>%
# #   ggplot(aes(x=year, y=percent_tot_genomes, fill=both_class)) + geom_col() +
# #   scale_fill_brewer(palette = 'Set1')+
# #   theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
# #   ggtitle('Co-occurrence of Metal and MDR islands, ONLY SOUTH',
# #           subtitle = 'SGI-4 1st, MDR 2nd') +
# #   labs(fill = "SGI4 / MDR") +
# #   geom_text(aes(y=50, label=tot_genomes_country)) + ylab('Percent of total genomes')
# 
# 
# 
# 
# ### THIS ONE ###
# master_meta %>%
#   filter(tree_clade =='south') %>% 
#   group_by(year, both_class) %>% 
#   summarise(number_genomes=n()) %>%
#   group_by(year) %>% mutate(tot_genomes_country=sum(number_genomes),
#                             percent_tot_genomes=(number_genomes/tot_genomes_country)*100) %>%
#   filter(as.numeric(year) > 2005) %>%
#   ggplot(aes(x=year, y=percent_tot_genomes, fill=both_class)) + geom_col() +
#   scale_fill_brewer(palette = 'Set1')+
#   theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
#   ggtitle('Co-occurrence of Metal and MDR islands, ONLY SOUTH',
#           subtitle = 'Metal 1st, MDR 2nd') +
#   labs(fill = "SGI4 / MDR") +
#   geom_text(aes(y=50, label=tot_genomes_country))+ 
#   ylab('Percent of genomes')
# 
# #
# #############################
# #########################
# 
# 
# # cooccur host?
# # 
# # master_meta %>% group_by(country, both_class3) %>%
# #   tally() %>% ggplot(aes(x=country, y=n, fill=both_class3)) + geom_col() +
# #   theme_bw()+scale_fill_brewer(palette = 'Set1')+
# #   theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
# #   labs(fill = "SGI4 / MDR") +
# #   ggtitle('Co-occurrence of Metal and MDR islands',
# #           subtitle = 'Metal 1st, MDR 2nd')
# # 
# 
# ### THIS ONE!
# master_meta %>% filter(!is.na(country)) %>% 
#   group_by(country, both_class3, tree_clade) %>% summarise(number_genomes=n()) %>%
#   group_by(country) %>% mutate(tot_genomes_country=sum(number_genomes),
#                                percent_tot_genomes=(number_genomes/tot_genomes_country)*100) %>%
#   # filter(!(country %in% c('unknown', 'South Korea', 'Brazil', 'Mexico'))) %>%
#   ggplot(aes(x=country, y=percent_tot_genomes, fill=both_class3)) + geom_col() +
#   scale_fill_brewer(palette = 'Set1')+
#   theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
#   ggtitle('Co-occurrence of SGI-4 and MDR islands',
#           subtitle = 'SGI-4 1st, MDR 2nd') +
#   labs(fill = "SGI4 / MDR") +
#   geom_text(aes(y=50, label=tot_genomes_country)) + 
#   ylab('Percent total genomes') + facet_wrap(~tree_clade, ncol=1)
# 
# 
# ##
# master_meta %>% filter(!is.na(Isolation)) %>% 
#   filter(!(Isolation %in% c('human', 'unknown', 'other', 'Food', 'environmental'))) %>%
#   group_by(Isolation, both_class, tree_clade) %>%
#   summarise(number_genomes=n()) %>%
#   ungroup() %>% 
#   ggplot(aes(x=Isolation, y=number_genomes, fill=both_class)) + geom_col() +
#   scale_fill_brewer(palette = 'Set1')+
#   theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
#   ggtitle('Co-occurrence of SGI4 and MDR module by host', 'excluding human isolates') +
#   xlab('') +  labs(fill = "SGI4 / MDR") +
#   ylab('number of genomes') + facet_wrap(~tree_clade)
# 
# 
# # percent of total
# master_meta %>%
#   filter(tree_clade == 'south') %>% 
#   group_by(Isolation, both_class) %>%
#   summarise(number_genomes=n()) %>%
#   ungroup() %>% 
#   group_by(Isolation) %>%
#   mutate(tot_genomes_country=sum(number_genomes),
#          percent_tot_genomes=(number_genomes/tot_genomes_country)*100) %>%
#   ungroup() %>% 
#   # filter(!(Isolation %in% c('human', 'unknown', 'other', 'Food', 'environmental'))) %>%
#   filter(!is.na(Isolation)) %>% 
#   ggplot(aes(x=Isolation, y=percent_tot_genomes, fill=both_class)) + geom_col() +
#   scale_fill_brewer(palette = 'Set1')+
#   theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
#   ggtitle('Co-occurrence of SGI4 and MDR module by host, only SOUTH', 'excluding human isolates') +
#   geom_text(aes(y=50, label=tot_genomes_country))+
#   xlab('') +  labs(fill = "SGI4 / MDR") +
#   ylab('Percent of total genomes')
# 
# 
# 
# #
# 
# 
# 
# 
# 
# ######################
# ##################
# 
# 
# # master_meta %>% group_by(state, SGI4_class) %>% 
# #   summarise(number_genomes=n()) %>% 
# #   group_by(state) %>% mutate(tot_genomes_state=sum(number_genomes), 
# #                              percent_tot_genomes = number_genomes/tot_genomes_state *100) %>% 
# #   ggplot(aes(x=state, y=percent_tot_genomes, fill=SGI4_class)) + geom_col()
# 
# 
# 
# #
# # FULL_SGI4_STATE <- master_meta2 %>% group_by(state, SGI4_class) %>% 
# #   summarise(number_genomes=n()) %>% 
# #   group_by(state) %>% mutate(tot_genomes_state=sum(number_genomes), 
# #                              percent_tot_genomes = number_genomes/tot_genomes_state *100) %>% 
# #   filter(SGI4_class == 'full')
# # 
# # 
# # 
# # FULL_MDR_STATE <- master_meta2 %>% group_by(state, MDR_class) %>% 
# #   summarise(number_genomes=n()) %>% 
# #   group_by(state) %>% mutate(tot_genomes_state=sum(number_genomes), 
# #                              percent_tot_genomes = number_genomes/tot_genomes_state *100) %>% 
# #   filter(MDR_class == 'full')
# 
# 
# ######## by state maps ##########
# 
# # 
# # states <- map_data("state")
# # 
# # state_names <- tolower(state.name)
# # names(state.abb) <- state_names
# # 
# # names(state_names) <- state.abb
# # master_meta$region <- state_names[master_meta$state]
# # 
# # state_plots <- master_meta %>% filter(!is.na(region))
# # 
# # 
# # 
# # FULL_SGI4_STATE <- master_meta %>% group_by(region, SGI4_class) %>% 
# #   summarise(number_genomes=n()) %>% 
# #   group_by(region) %>% mutate(tot_genomes_state=sum(number_genomes), 
# #                               percent_tot_genomes = number_genomes/tot_genomes_state *100) %>% 
# #   filter(SGI4_class == 'full')
# # 
# # 
# # 
# # FULL_MDR_STATE <- master_meta %>% group_by(region, MDR_class) %>% 
# #   summarise(number_genomes=n()) %>% 
# #   group_by(region) %>% mutate(tot_genomes_state=sum(number_genomes), 
# #                               percent_tot_genomes = number_genomes/tot_genomes_state *100) %>% 
# #   filter(MDR_class == 'full')
# # 
# # 
# # 
# # FULL_BOTH_STATE <- master_meta %>% group_by(region, both_class) %>% 
# #   summarise(number_genomes=n()) %>% 
# #   group_by(region) %>% mutate(tot_genomes_state=sum(number_genomes), 
# #                               percent_tot_genomes = number_genomes/tot_genomes_state *100) %>% 
# #   filter(both_class == 'full / full')
# # 
# # 
# # states %>% left_join(FULL_SGI4_STATE) %>% ggplot() +
# #   geom_polygon(aes(x = long, y = lat, fill = number_genomes, group = group), color = "white") + 
# #   coord_fixed(1.3) + scale_fill_viridis_c() + ggtitle('Number of 4,[5],12:i:- genomes with SGI-4')
# # 
# # 
# # states %>% left_join(FULL_MDR_STATE) %>% ggplot() +
# #   geom_polygon(aes(x = long, y = lat, fill = number_genomes, group = group), color = "white") + 
# #   coord_fixed(1.3) + scale_fill_viridis_c() + ggtitle('Number of 4,[5],12:i:- genomes with full MDR module')
# # 
# # 
# # states %>% left_join(FULL_BOTH_STATE) %>% ggplot() +
# #   geom_polygon(aes(x = long, y = lat, fill = number_genomes, group = group), color = "white") + 
# #   coord_fixed(1.3) + scale_fill_viridis_c() + ggtitle('Number of 4,[5],12:i:- genomes with SGI-4 and MDR module')
# # 
# 
# #######
# # investigate leiden cluster #
# # seems like most all chicken genomes belong to cluster 1. some cow some turkey
# # cluster 2 and 3 are associated strongly with swine, cluster 2 also turkey
# # cluster 3 strongly associated with swine.
# # 
# # 
# # master_meta %>% group_by(Isolation , leid_clust, year) %>% 
# #   tally() %>% ungroup() %>% 
# #   group_by(Isolation, leid_clust) %>% 
# #   mutate(p_genomes=n/sum(n) * 100) %>% 
# #   ggplot(aes(x=Isolation, fill=year, weight=p_genomes)) + geom_bar() + 
# #   facet_wrap(~leid_clust, scales = 'free')
# # 
# # 
# # master_meta %>% group_by(Isolation , leid_clust, year) %>% 
# #   tally() %>% ungroup() %>%
# #   filter(!(year %in% c('unknown', NA))) %>% 
# #   # filter(leid_clust != 0) %>%
# #   group_by(Isolation, leid_clust) %>% 
# #   mutate(p_genomes=n/sum(n) * 100, 
# #          grp=paste(Isolation, leid_clust, sep='_'), 
# #          year_num=as.numeric(year), 
# #          leid_clust=as.factor(leid_clust)) %>% 
# #   ungroup() %>% 
# #   #filter(year>2000) %>% 
# #   ggplot(aes(x=year_num, y=n, color=leid_clust, group=grp)) + geom_line() +
# #   facet_wrap(~Isolation, scales = 'free_y') + scale_color_brewer(palette = 'Set1', direction = -1)
# # 
# # ### ABX patterns ###
# # master_meta$num_resist[is.na(master_meta$num_resist)] <- 0
# # 
# # 
# # master_meta %>% filter(!(is.na(MDR_class))) %>% 
# #   ggplot(aes(x=tree_clade, y=num_resist)) + 
# #   geom_violin() +
# #   geom_jitter(aes(color=MDR_class),alpha=.3, height = 0) +
# #   ggtitle('The number of antibiotics each genome has predicted resistance to')
# # #TODO look at AMR not related to MDR module
# # 
# # 
# # 
# # 
# # # 
# # # unique(master_meta$res_type)
# # # table(master_meta$res_type)
# # 
# # master_meta %>% filter(tree_clade == 'south') %>% 
# #   # mutate(RES=fct_lump(res_type, n=20)) %>% 
# #   group_by(res_type) %>% tally() %>% arrange(desc(n))
# # 
# # 
# # master_meta %>% filter(tree_clade == 'north') %>% 
# #   # mutate(RES=fct_lump(res_type, n=20)) %>% 
# #   group_by(res_type) %>% tally() %>% arrange(desc(n))
# # 
# # master_meta %>% filter(tree_clade == 'other') %>% 
# #   # mutate(RES=fct_lump(res_type, n=20)) %>% 
# #   group_by(res_type) %>% tally() %>% arrange(desc(n))
# 

#### accessory genome clusters ###

acc_clade_dat <- 
  master_meta %>%
  group_by(tree_clade,leiden_cluster) %>%
  tally() %>% 
  mutate(prop_clade=n/sum(n), 
         perc_clade=prop_clade *100) %>% 
  filter(tree_clade %in% c('clade_one', 'clade_two', 'clade_three'))

acc_clade_dat %>% 
  ggplot(aes(x=tree_clade, y=sqrt(n), fill=as.character(leiden_cluster))) +
  geom_col()

master_meta %>% 
  filter(tree_clade == 'clade_one' & leiden_cluster == 'B' |
           tree_clade == 'clade_three' & leiden_cluster == 'C') %>%
  dplyr::select(year, state, country, ISOLATION, isolation_source)

colnames(master_meta)


master_meta %>% filter(tree_clade == 'clade_three' & leiden_cluster == 'C') %>% pull(year)

######


save.image(file='FIGS.RData')


# interpro <-
#   read_tsv('./data/july6_bigpan_interproscan.tsv',
#            col_names = c('accno', 'MD5', 'len', 'analysis', 
#                          'sig_accno', 'sig_description', 
#                          'start', 'stop', 'score', 'status', 
#                          'date', 'interpro_anno_accno', 
#                          'interpro_anno_desc', 'GO', 'pathway')) %>% 
#   select(accno, GO) %>% 
#   filter(!is.na(GO))
# 
# ONLY_NORTH <- north_genes$accno[!(north_genes$accno %in% south_genes$accno)]
# ONLY_SOUTH <- south_genes$accno[!(south_genes$accno %in% north_genes$accno)]


### CHecking BRADS NUMBERS ###

master_meta %>% group_by(targ_presence_SGI4) %>% 
  tally() %>% 
  mutate(perc=n/sum(n))

master_meta %>% nrow()



look <- master_meta %>%
  count(country) %>% 
  mutate(perc=(n/sum(n))*100)

master_meta %>% count(tree_clade) %>%
  mutate(perc=(n/sum(n))*100)



master_meta %>%
  count(continent) %>% 
  mutate(perc=(n/sum(n))*100)


master_meta %>% count(closest_ref)
master_meta %>% 
  group_by(MLST, closest_ref) %>% 
  tally() %>% arrange(desc(n))

LOOK <- master_meta %>% 
  group_by(MLST) %>% 
  tally() %>% 
  arrange(desc(n))

master_meta %>% group_by(Source) %>% tally() %>%
  filter(!(Source %in% c('Human', 'Other'))) %>% 
  mutate(perc=(n/sum(n))*100)


master_meta %>% group_by(targ_presence_MDR, tree_clade) %>% tally()


master_meta$AMR_genotypes
master_meta$stress_genotypes


MDR_GENES <- 
  grepl('merA',master_meta$stress_genotypes) &
grepl('merC',master_meta$stress_genotypes)&
grepl('merD',master_meta$stress_genotypes)&
grepl('merE',master_meta$stress_genotypes)&
grepl('merP',master_meta$stress_genotypes)&
grepl('merR',master_meta$stress_genotypes)&
grepl('merT',master_meta$stress_genotypes)&
grepl('blaTEM-1',master_meta$AMR_genotypes)&
grepl("aph(3'')-Ib",master_meta$AMR_genotypes, fixed = T) &
grepl("aph(6)-Id",master_meta$AMR_genotypes, fixed = T) &
grepl("sul2",master_meta$AMR_genotypes, fixed = T) &
grepl("tet(B)",master_meta$AMR_genotypes, fixed = T)


sum(MDR_GENES)
master_meta$MDR_GENES <- MDR_GENES



master_meta %>% group_by(tree_clade, MDR_GENES) %>% tally()


master_meta %>% group_by(Source, MDR_GENES) %>% tally() %>% mutate(perc=(n/sum(n)*100)) %>% 
  filter(MDR_GENES ==T)

master_meta %>% group_by(targ_presence_SO4698_sopE_mTmV) %>% tally()





grepl('blaCMY-2',master_meta$AMR_genotypes) %>% sum()
grepl('blaSHV-12',master_meta$AMR_genotypes) %>% sum()
grepl('blaCTX-M',master_meta$AMR_genotypes) %>%  sum()

master_meta$blacmy2 <- grepl('blaCMY-2',master_meta$AMR_genotypes)
master_meta$blashv12 <- grepl('blaSHV-12',master_meta$AMR_genotypes)
master_meta$blactxm <- grepl('blaCTX-M',master_meta$AMR_genotypes)


master_meta %>% 
  group_by(blacmy2, tree_clade) %>%
  tally() %>% 
  mutate(perc=(n/sum(n))*100)


master_meta %>% 
  filter(tree_clade == 'clade_one') %>% 
  group_by(blacmy2, Source) %>%
  tally() %>% 
  mutate(perc=(n/sum(n))*100)

# The number (%) of Salmonella enterica serovar I 4,[5],12:i:- strain sequences based on country of origin, closest reference strain, and clade membership

master_meta %>% 
  count(country, closest_ref, tree_clade) %>% 
  mutate(percent_of_tot=(n/sum(n))*100) %>% 
  transmute(Country=country, 
            closest_reference = closest_ref, 
            phylogenetic_clade=tree_clade, 
            percent_of_tot=signif(percent_of_tot, 3)) %>% 
  arrange(desc(percent_of_tot)) %>% 
  write_tsv('Table_S2.tsv')
