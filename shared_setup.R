# Libraries
library(tidyverse)
library(magrittr)
library(pheatmap)
library(caret)
library(phylolm)
library(treeio)
library(ggtree)
library(ape)
library(aplot)
library(ggVennDiagram)
library(ggsignif)
library(ggstance)
library(ranger)
library(cowplot)
library(exactRankTests)

# Source
source("scripts/analysis_wrapper.R")
source("phylotree/read_raxml_tree.R")
source("phylotree/tree_funs.R")
source("phylotree/phyloglmstep_fixed.R")
source("phylotree/comparative_fixed.data.R")
set.seed(1)

# define color codes
color_scale <-
  c("poplar endophyte"="#F8766D",
    "endophyte"="red",
    "saprobe"="#00BFC4",
    "canker"="purple",
    "ERMF"="gold",
    "fungal parasite"="black",
    "pathogen"="green",
    "ECM"="brown"
  )

ecology_order <-
  c("poplar endophyte",
    "endophyte",
    "ERMF",
    "ECM",
    "pathogen",
    "canker",
    "saprobe",
    "fungal parasite"
  )

# Read data
## taxa
taxa_tb <- read_csv("meta_tables/genome_relatives.csv") %>%
  select(
    spcode=portal_ID,
    strain,
    phylum,
    class,
    order,
    family,
    genus,
    species
  ) %>%
  column_to_rownames("spcode") %>%
  as_tibble(rownames = "spcode")

## genome metadata and ecology
all_genome_tb <-
  read_csv("meta_tables/genome_relatives.csv") %>%
  rename(spcode=portal_ID)

eco_tb <-
  all_genome_tb %>%
  select(spcode, ecology) %>%
  mutate(ecology=replace(ecology, ecology=="poplar endo", "poplar endophyte")) %>%
  mutate(ecology=as.factor(ecology))

contrast_tb <-
  eco_tb %>%
  left_join(select(taxa_tb, spcode, phylum, class),
            by="spcode") %>%
  select(-phylum, -class) %>%
  arrange(ecology) %>%
  mutate(ecology=fct_collapse(
    ecology,
    `poplar endophyte`=c("poplar endophyte"),
    saprobe="saprobe",
    other_level = "Other"
  )
  ) %>%
  filter(ecology!="Other") %>%
  mutate(ecology=fct_drop(ecology))

contrast_tb <-
  contrast_tb %>%
  filter(!(spcode %in% c("Amore1", "Hesve2finisherSC"))) %>% #extra saprobe
  filter(!(spcode %in% c("Gloci1"))) %>% #dulicate in one lineage
  filter(!(spcode %in% c("Ilysp1","MarPMI226"))) %>% #dulicate in one lineage
  filter(!(spcode %in% c("Hyafin1"))) #dulicate in one lineage
