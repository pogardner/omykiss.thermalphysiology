# Testing for variation in adaptive panel
library(tidyverse)
library(here)
adap_mhap <- read.csv(here("Genomics/Data/Adaptive_genotypes_mhapcode_wide.csv"))

# metadata
ctmax <- read.csv(here("Genomics/Data/Omy5_inversion_metadata_merge.csv"))

adap_mhap_pops <- ctmax %>% 
  rename(indiv.ID = NMFS_DNA_ID_use) %>% 
  filter(!indiv.ID %in% c("M145767","M145542", 
                          "M145524", "M145692")) %>% 
  dplyr::select(c(indiv.ID,WATERSHED)) %>% 
  left_join(., adap_mhap, by ="indiv.ID")

nc <- ncol(adap_mhap)
loci <- str_replace(names(adap_mhap)[seq(2, nc, by = 2)], "\\.1$", "")
# then make some long format genotypes
long_genos_adap <- adap_mhap %>% 
  gather(key = "loc", value = "Allele", -indiv.ID) %>%
  extract( # ChatGPT helped me with this code chunk. had to adapt from CKMR tutorial bc some of the gene names have other . in them
    loc,
    into = c("Locus", "gene_copy"),
    regex = "^(.*)\\.([12])$"
  ) %>%
  mutate(Allele = as.character(Allele)) %>%
  mutate(Allele = ifelse(Allele == "0", NA, Allele)) %>%
  rename(Indiv = indiv.ID)

alle_freqs_adap <- long_genos_adap %>%
  count(Locus, Allele) %>%
  group_by(Locus) %>%
  mutate(Freq = n / sum(n),
         Chrom = "Unk",
         Pos = as.integer(factor(Locus, levels = loci))) %>%
  ungroup() %>%
  dplyr::select(Chrom, Pos, Locus, Allele, Freq) %>%
  arrange(Pos, desc(Freq)) %>%
  mutate(AlleIdx = NA,
         LocIdx = NA) %>%
  filter(!is.na(Allele))

# I want to visualize genotypes. I also want to color it by location - maybe all individuals in each location are homogenous in genotypes

# Let's look at each gene group individually.
colnames(adap_mhap_pops)
adap_mhap[adap_mhap ==0] <- NA
ch28 <- adap_mhap_pops %>% 
  dplyr::select(starts_with(c("Omy_Ch28_11","indiv.ID","WATER")))
six6<- adap_mhap_pops %>% 
  dplyr::select(starts_with(c("Omy_Ch25_Six6","indiv.ID","WATER")))
greb <-adap_mhap_pops %>% 
  dplyr::select(starts_with(c("Omy_Ch28_greb1","indiv.ID","WATER")))
omy <- adap_mhap_pops %>% 
  dplyr::select(starts_with(c("Omy_Ch05","indiv.ID","WATER")))
vgll3 <- adap_mhap_pops %>% 
  dplyr::select(starts_with(c("Omy_Ch22_VGLL3","indiv.ID","WATER")))



# Ch 28
ch28_long <- ch28 %>%
  pivot_longer(
    cols = -c(indiv.ID,WATERSHED),
    names_to = "Locus",
    values_to = "Allele"
  )
ggplot(ch28_long,aes(x=Locus, y = Allele, color = indiv.ID))+ 
   geom_jitter()+
   labs(title = "Ch 28")+
  facet_wrap(~WATERSHED)+
   theme(legend.position = "none",
         axis.text.x = element_text(angle = 45, hjust = 1, size = 5))
 

# six6
six6_long <- six6 %>%
  pivot_longer(
    cols = -c(indiv.ID,WATERSHED),
    names_to = "Locus",
    values_to = "Allele"
  )
ggplot(six6_long,aes(x=Locus, y = Allele, color = indiv.ID))+ 
   geom_jitter()+
   labs(title = "six6")+
  facet_wrap(.~WATERSHED)+
   theme(legend.position = "none",
         axis.text.x = element_text(angle = 45, hjust = 1))

# vgll3
vgll3_long <- vgll3 %>%
  pivot_longer(
    cols = -c(indiv.ID,WATERSHED),
    names_to = "Locus",
    values_to = "Allele"
  )
ggplot(vgll3_long,aes(x=Locus, y = Allele, color = indiv.ID))+ 
  geom_jitter()+
  labs(title = "vgll3")+
  facet_wrap(.~WATERSHED)+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Omy05
omy_long <- omy %>%
  pivot_longer(
    cols = -c(indiv.ID,WATERSHED),
    names_to = "Locus",
    values_to = "Allele"
  )
ggplot(omy_long,aes(x=Locus, y = Allele, color = indiv.ID))+ 
  geom_jitter()+
  labs(title = "Omy05")+
  facet_wrap(.~WATERSHED)+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# greb
greb_long <- greb %>%
  pivot_longer(
    cols = -c(indiv.ID,WATERSHED),
    names_to = "Locus",
    values_to = "Allele"
  )
ggplot(greb_long,aes(x=Locus, y = Allele, color = indiv.ID))+ 
  geom_jitter()+
  labs(title = "greb1")+
  facet_wrap(.~WATERSHED)+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
