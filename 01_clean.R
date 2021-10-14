library(tidyverse)

master_meta <- read_tsv('./data/JAN_UPDATE_META.tsv', guess_max = 5000) #%>%




master_meta$country <- sub('(.*):(.*)','\\1',master_meta$geo_loc_name)
# unique(master_meta$geo_loc_name)
# unique(master_meta$country)
master_meta <-
  master_meta %>%
  filter(country != 'not collected') %>% 
  mutate(continent=case_when(
    country %in% c('USA', 'Canada', 'Mexico')   ~ 'North America', 
    country %in% c('Australia')                 ~ 'Oceania', 
    country %in% c('Brazil', 'Chile',
                   'Colombia', 'Ecuador',
                   'Barbados','Guyana',
                   'Trinidad and Tobago')       ~ 'South America', 
    country %in% c('Belgium', 'Denmark', 
                   'France', 'Germany',
                   'Ireland', 'Italy',
                   'Poland', 'Spain', 
                   'Sweden', 'Switzerland',
                   'United Kingdom')           ~ 'Europe', 
    country %in% c('Cambodia', 'China', 
                   'Japan', 'Laos', 
                   'Lebanon', 'Taiwan', 
                   'Thailand', 'Viet Nam')     ~ 'Asia', 
    TRUE ~ 'MISSED ME'
  ))



table(master_meta$country)
master_meta %>% filter(continent == 'MISSED ME')
# north america
# south america
# asia
# europe
# africa
# Oceania

###
master_meta <- 
  master_meta %>% 
  mutate(global_area=case_when(
    continent %in% c('North America', 'South America')      ~ 'Americas',
    continent %in% c('Europe', 'Asia', 'Oceania', 'Africa') ~ 'Eurasia',
    TRUE                                                    ~ 'ERROR'))




###


# remove genomes with no country?
master_meta <- master_meta[!is.na(master_meta$country),]

master_meta$state <- 'unknown'
master_meta$state[grep(':', master_meta$geo_loc_name)] <- sub('(.*):(.*)','\\2', master_meta[grep(':', master_meta$geo_loc_name),]$geo_loc_name)
master_meta$state <- trimws(master_meta$state, which = 'both')
master_meta$state[grep('United Kingdom', master_meta$state)] <- NA
master_meta$state[grep('Chester', master_meta$state)] <- 'PA'
master_meta$state[grep('Lancaster', master_meta$state)] <- 'PA'
master_meta$state[grep('Knoxville', master_meta$state)] <- 'TN'
master_meta$state[grep('Nassau', master_meta$state)] <- 'NY'
master_meta$state[grep('SOUTH DAKOTA', master_meta$state)] <- 'SD'
master_meta$state[grep('MINNESOTA', master_meta$state)] <- 'MD'
master_meta$state[grep('New Mexico', master_meta$state)] <- 'NM'



# unique(master_meta$state)




# state is done #


SOURCE_WRANGLE <- 
  master_meta %>% 
  transmute(
    genome=genome, 
    host=tolower(host), 
    isolation_source=tolower(isolation_source), 
    IFSAC_category=tolower(IFSAC_category), 
    epi_type=tolower(epi_type), 
    ontological_term=ontological_term, 
    source_type=tolower(source_type), 
    collected_by=tolower(collected_by))
  



swinenames <- c('swine','pork','porcine','sow','sus','hog','pig', 'scrofa')
cownames <- c('bovine','beef','veal','cow','cattle','bos','steer', 'taurus', 'calf')
chickennames <- c('chicken', 'chick', 'gallus', 'broiler', 'egg')
turkeynames <- c('turkey', 'meleagris', 'gallopavo')
humannames <- c('human', 'homo', 'sapiens')
horsenames<- c('equine', 'equus', 'horse', 'caballus')
dognames<- c('canine', 'dog', 'canis')
catnames<- c('cat', 'felis', 'catus')

cleaner <- 
  function(column, matches_vector, return_value){
  PATTERN <- paste(matches_vector, collapse = '|')
  MATCHES <- grepl(pattern = PATTERN, x = column)
  column[MATCHES] <- return_value
  return(column)
}



# first pass through isolation_source column
SOURCE_WRANGLE <- 
  SOURCE_WRANGLE %>%
  mutate(isotmp=cleaner(column = isolation_source,
                        matches_vector = swinenames, 
                        return_value = 'swine'), 
         isotmp=cleaner(column=isotmp, 
                        matches_vector = cownames, 
                        return_value = 'cow'), 
         isotmp=cleaner(column = isotmp, 
                        matches_vector = chickennames, 
                        return_value = 'chicken'), 
         isotmp=cleaner(column = isotmp, 
                        matches_vector = turkeynames, 
                        return_value = 'turkey'),
         isotmp=cleaner(column = isotmp, 
                        matches_vector = humannames, 
                        return_value = 'human'), 
         isotmp=cleaner(column = isotmp, 
                        matches_vector = catnames, 
                        return_value = 'cat'), 
         isotmp=cleaner(column = isotmp, 
                        matches_vector = dognames, 
                        return_value = 'dog'), 
         isotmp=cleaner(column = isotmp, 
                        matches_vector = horsenames, 
                        return_value = 'horse')
    
)


finished <- SOURCE_WRANGLE %>%
  filter(isotmp %in% c('swine', 'cow', 'chicken', 'turkey', 'human', 'cat', 'dog', 'horse'))

not_done <- SOURCE_WRANGLE %>%
  filter(!(isotmp %in% c('swine', 'cow', 'chicken', 'turkey', 'human', 'cat', 'dog', 'horse')))



not_done <- 
  not_done %>% 
  mutate(isotmp=cleaner(column = host,
                        matches_vector = swinenames, 
                        return_value = 'swine'), 
         isotmp=cleaner(column=isotmp, 
                        matches_vector = cownames, 
                        return_value = 'cow'), 
         isotmp=cleaner(column = isotmp, 
                        matches_vector = chickennames, 
                        return_value = 'chicken'), 
         isotmp=cleaner(column = isotmp, 
                        matches_vector = turkeynames, 
                        return_value = 'turkey'),
         isotmp=cleaner(column = isotmp, 
                        matches_vector = humannames, 
                        return_value = 'human'), 
         isotmp=cleaner(column = isotmp, 
                        matches_vector = catnames, 
                        return_value = 'cat'), 
         isotmp=cleaner(column = isotmp, 
                        matches_vector = dognames, 
                        return_value = 'dog'), 
         isotmp=cleaner(column = isotmp, 
                        matches_vector = horsenames, 
                        return_value = 'horse')
         
  ) 

finished2 <- 
  not_done %>% 
  filter(isotmp %in% c('swine', 'cow', 'chicken', 'turkey', 'human', 'cat', 'dog', 'horse'))



not_done2 <- not_done %>%
  filter(!(isotmp %in% c('swine', 'cow', 'chicken', 'turkey', 'human', 'cat', 'dog', 'horse')))


not_done2 <- 
  not_done2 %>% 
  mutate(isotmp=ifelse(epi_type == 'clinical' & collected_by == 'cdc', 'human', isotmp), 
         isotmp=ifelse(epi_type == 'clinical' & grepl('public health', collected_by), 'human', isotmp)) 

finished3 <- not_done2 %>% 
  filter(isotmp %in% c('swine', 'cow', 'chicken', 'turkey', 'human', 'cat', 'dog', 'horse'))

not_done3 <- not_done2 %>% 
  filter(!(isotmp %in% c('swine', 'cow', 'chicken', 'turkey', 'human', 'cat', 'dog', 'horse')))


not_done3 <- 
  not_done3 %>% mutate(isotmp=
                        ifelse(epi_type == 'clinical' & grepl('stool|urine', isolation_source), 'human', isotmp), 
                       isotmp=
                         ifelse(collected_by == 'isuvdl', 'swine', isotmp))


finished4 <- not_done3 %>% 
  filter(isotmp%in% c('swine', 'cow', 'chicken', 'turkey', 'human', 'cat', 'dog', 'horse'))

not_done4 <-  not_done3 %>% 
  filter(!(isotmp%in% c('swine', 'cow', 'chicken', 'turkey', 'human', 'cat', 'dog', 'horse')))

table(not_done4$isolation_source)
table(not_done4$host)
table(not_done4$IFSAC_category)

finished5 <- not_done4 %>% mutate(isotmp='other')


cleaned_isolations <- 
  bind_rows(finished, finished2, finished3, finished4, finished5) %>%
  transmute(genome=genome,
            ISOLATION=factor(isotmp)) %>% 
  mutate(ISOLATION=fct_lump(ISOLATION, n = 6, other_level = 'other'))


cleaned_isolations %>%
  count(ISOLATION) %>% 
  arrange(desc(n))

#### year?

master_meta$year <- as.numeric(sub('(....).*','\\1',master_meta$collection_date))

tmpyear <- as.numeric(sub('(....).*','\\1',master_meta$target_creation_date))

master_meta$year <- ifelse(is.na(master_meta$year), tmpyear, master_meta$year)




####




master_meta <-
  master_meta %>% 
  left_join(cleaned_isolations) %>%
  write_tsv('./output/01_meta_JAN_UPDATE.tsv')

#################
# for dereplication #


master_meta %>%
  group_by(PDS_acc, ISOLATION) %>%
  tally() %>%
  ungroup() %>%  
  arrange(desc(n)) %>%
  filter(ISOLATION != 'human') %>% 
  filter(ISOLATION == 'swine')

# 
# master_meta %>% group_by(PDS_acc, year, country, ISOLATION) %>% tally() %>% 
#   arrange(desc(n))
# 
# leid <- read_csv('./output/leiden_R1_clust.csv')
# 
# leid %>% group_by(clusts) %>% tally()
# 
# LOOK <- leid %>% select(genome, clusts) %>% left_join(master_meta)
# 
# 
# weirdos <- LOOK %>% filter(!(clusts %in% c(0,1))) %>% pull(genome)
# 
# LOOK2 <- LOOK %>%
#   group_by(PDS_acc, clusts, ISOLATION) %>% 
#   tally() %>%
#   arrange(desc(n)) %>% 
#   filter(!is.na(ISOLATION))
# 
# 
# LOOK
# 
# test <- matr %>% select(Gene, all_of(weirdos))
# colSums(matr2)
# #### SNP TREE


