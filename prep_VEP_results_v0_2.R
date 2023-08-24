
# this is a script for Lusanda to prep his data

library(tidyverse)


# put file path in Lusanda .Rmd

VEP_results <- read_delim(file_path)

# main aim is to fix duplicates

names(VEP_results)[1] <- gsub("#", "X.", names(VEP_results)[1])

length(unique(VEP_results$X.Uploaded_variation))

VEP_results %>% 
  select(X.Uploaded_variation, Location) %>% 
  distinct() %>% 
  group_by(X.Uploaded_variation) %>% 
  count() %>% 
  filter(n != 1) %>%
  inner_join(VEP_results %>% 
               select(X.Uploaded_variation, Location)) %>% 
  print(n = Inf)

VEP_results %>% 
  select(X.Uploaded_variation, Location) %>% 
  distinct() 


VEP_results <- VEP_results %>% 
  distinct()


# split out the allele frequencies
names(VEP_results)

VEP_results_AF <- VEP_results %>% 
  select(X.Uploaded_variation, Location,
         Allele, Codons, contains("AF"))

VEP_results <- VEP_results %>% 
  select(-contains("AF")) %>% 
  distinct()

# find which other variables are the cause of duplicates

x <- map_dbl(names(VEP_results),
        ~VEP_results %>% 
          select(all_of(c("X.Uploaded_variation", .x))) %>% 
          distinct() %>% 
          group_by(X.Uploaded_variation) %>% 
          count() %>% 
          filter(n != 1) %>% 
          nrow())

x <- which(x != 0)
x <- names(VEP_results)[c(1, x)]

VEP_results %>% 
  select(all_of(x)) %>% 
  distinct() %>% 
  group_by(X.Uploaded_variation) %>% 
  count() %>% 
  filter(n != 1) %>%
  inner_join(VEP_results %>% 
               select(all_of(x))) %>% 
  print(n = Inf)

# remove the duplicate lines for these IDs with only dashes in PHENOTYPE

VEP_results %>% 
  select(X.Uploaded_variation, PHENOTYPES) %>% 
  distinct() %>% 
  print(n = Inf)

x <- VEP_results %>% 
  select(all_of(x)) %>% 
  distinct() %>% 
  group_by(X.Uploaded_variation) %>% 
  count() %>% 
  filter(n != 1) %>% 
  pull(X.Uploaded_variation)

VEP_results <- VEP_results %>% 
  filter(!(X.Uploaded_variation %in% x & PHENOTYPES == "-"))

###

x <- c("X.Uploaded_variation", "PHENOTYPES")

VEP_results %>% 
  select(all_of(x)) %>% 
  distinct() %>% 
  group_by(X.Uploaded_variation) %>% 
  count() %>% 
  filter(n != 1) %>%
  inner_join(VEP_results %>% 
               select(all_of(x))) %>% 
  print(n = Inf) %>% 
  pull(PHENOTYPES)

x <- str_extract_all(VEP_results$PHENOTYPES, ",")
x <- max(map_dbl(x, ~length(.x)))

x <- paste0("xx", 1:x)

tmp <- VEP_results %>% 
  select(PHENOTYPES) %>% 
  separate(col = PHENOTYPES, into = x, sep = ",") %>% 
  pivot_longer(everything(),names_to = "key", values_to = "value") %>% 
  mutate(across(value, ~tolower(.x))) %>% 
  filter(!is.na(value) & value != "-") %>% 
  select(-key) %>% 
  distinct() %>% 
  print(n = Inf)

# organise PHENOTYPES into tables ---------------------------------------------

# park not mentioned

pheno_nopd <- tmp %>% 
  filter(!str_detect(value, "park")) %>% 
  mutate(pheno_consolidated = "no_pd_mention")

tmp <- tmp %>% 
  filter(str_detect(value, "park")) 


# wolff-parkinson-white

pheno_wpw <- tmp %>% 
  filter(str_detect(value, "wolff")) %>% 
  mutate(pheno_consolidated = "wolff_parkinson_white")

tmp <- tmp %>% 
  filter(!str_detect(value, "wolff")) 

# parkinsonism

pheno_parkinsonism <- tmp %>% 
  filter(str_detect(value, "parkinsonism")) %>% 
  mutate(pheno_consolidated = "parkinsonism")

tmp <- tmp %>% 
  filter(!str_detect(value, "parkinsonism"))

tmp %>% 
  print(n = Inf)

# late onset parkinson disease

pheno_lateonsetpd <- tmp %>% 
  filter(str_detect(value, "parkinson_disease") &
           str_detect(value, "late") &
           str_detect(value, "onset")) %>% 
  mutate(pheno_consolidated = "pd_late_onset")
  
tmp <- tmp %>% 
  filter(!(str_detect(value, "parkinson_disease") &
           str_detect(value, "late") &
           str_detect(value, "onset")))

# young, juvenile, early, 

pheno_earlyonsetpd <- tmp %>% 
  filter(str_detect(value, "parkinson_disease") &
           (str_detect(value, "young") |
           str_detect(value, "early") |
             str_detect(value, "juvenile"))) %>% 
  mutate(pheno_consolidated = "pd_early_onset")

tmp <- tmp %>% 
  filter(!(str_detect(value, "parkinson_disease") &
           (str_detect(value, "young") |
              str_detect(value, "early") |
              str_detect(value, "juvenile")))) 

# parkinson disease, no age category

pheno_noagecategorypd <- tmp %>% 
  filter((str_detect(value, "parkinson_disease") |
           str_detect(value, "parkinson's_disease")) &
           !(str_detect(value, "young") |
              str_detect(value, "early") |
              str_detect(value, "juvenile") |
               str_detect(value, "late-onset"))) %>% 
  mutate(pheno_consolidated = "pd_no_age_category")

tmp <- tmp %>% 
  filter(!(str_detect(value, "parkinson_disease") |
            str_detect(value, "parkinson's_disease")) &
           !(str_detect(value, "young") |
               str_detect(value, "early") |
               str_detect(value, "juvenile") |
               str_detect(value, "late-onset"))) %>% 
  print(n = Inf)

# parkinson-dementia

pheno_parkdementiasynd <- tmp %>% 
  filter(str_detect(value, "parkinson-dementia")) %>%  
  mutate(pheno_consolidated = "park_dementia_synd")

tmp <- tmp %>% 
  filter(!str_detect(value, "parkinson-dementia")) %>% 
  print(n = Inf)

# parkinsonian

pheno_parkinsonian <- tmp %>% 
  filter(str_detect(value, "parkinsonian")) %>%  
  mutate(pheno_consolidated = "parkinsonian")

# other

pheno_other <- tmp %>% 
  filter(!str_detect(value, "parkinsonian")) %>%  
  mutate(pheno_consolidated = "other")

pheno_tables <- list(pheno_earlyonsetpd = pheno_earlyonsetpd, 
                     pheno_lateonsetpd = pheno_lateonsetpd,
                     pheno_noagecategorypd = pheno_noagecategorypd, 
                     pheno_parkdementiasynd = pheno_parkdementiasynd,
                     pheno_parkinsonian = pheno_parkinsonian,
                     pheno_parkinsonism = pheno_parkinsonism,
                     pheno_wpw = pheno_wpw, pheno_other = pheno_other,
                     pheno_nopd = pheno_nopd)

x <- paste0("pheno_tables_", chip, ".rds")
write_rds(pheno_tables, x) 

# Add pheno tables to data

x <- str_extract_all(VEP_results$PHENOTYPES, ",")
x <- max(map_dbl(x, ~length(.x)))


x <- paste0("xx", 1:x)

names(VEP_results)
glimpse(VEP_results)

keep <- c("X.Uploaded_variation" , "Existing_variation",  
          "Location" ,"Allele" ,"Consequence"          
          ,"IMPACT" ,"SYMBOL","Gene"                          
          ,"BIOTYPE","Amino_acids","Codons"     
          ,"STRAND" ,"SYMBOL_SOURCE"        
          ,"HGNC_ID" ,"SIFT"  ,"PolyPhen"                  
          ,"CLIN_SIG" ,"PUBMED" ,"CADD_PHRED","CADD_RAW", "pheno_consolidated")


VEP_results <- VEP_results %>% 
  separate(col = PHENOTYPES, into = x, sep = ",") %>% 
  pivot_longer(starts_with("xx"),names_to = "key", values_to = "value") %>% 
  mutate(across(value, ~tolower(.x))) %>% 
  filter(!is.na(value) & value != "-") %>% 
  select(-key, value) %>% 
  distinct() %>% 
  inner_join(pheno_tables %>% 
  bind_rows()) %>% 
  distinct() %>%  
  select(all_of(keep)) %>% 
  filter(pheno_consolidated != "no_pd_mention") %>% 
  distinct()

VEP_results_pheno <- VEP_results %>% 
  select(X.Uploaded_variation, Allele, SYMBOL, pheno_consolidated) %>% 
  distinct() 

VEP_results <- VEP_results %>% 
  select(-pheno_consolidated) %>% 
  distinct()
  

# evaluate duplicates again

x <- map_dbl(names(VEP_results),
             ~VEP_results %>% 
               select(all_of(c("X.Uploaded_variation", .x))) %>% 
               distinct() %>% 
               group_by(X.Uploaded_variation) %>% 
               count() %>% 
               filter(n != 1) %>% 
               nrow())

x <- which(x != 0)
x <- names(VEP_results)[c(1, x)]
x

VEP_results[,x] %>% 
  group_by(X.Uploaded_variation) %>% 
  count() %>% 
  filter(n != 1) %>% 
  inner_join(VEP_results[,x])

# duplicates caused by different alleles


# allele frequencies

VEP_results_AF



VEP_results_pheno

x <- paste0("rds/VEP_results_", chip, ".rds")
write_rds(VEP_results, x)
x <- paste0("rds/VEP_results_pheno_", chip, ".rds")
write_rds(VEP_results_pheno, x)
x <- paste0("rds/VEP_results_AF_", chip, ".rds")
write_rds(VEP_results_AF, x)








