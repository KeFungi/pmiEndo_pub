#' extract and make portal id column from portal_link
library(tidyverse)
library(magrittr)

intable <- read_csv("JGI_SMC_raw.csv")
intable %<>% mutate(portal_id=map_chr(portal_link, ~str_match(.x, "[^\\/]+$")[,1]))
intable %<>% select(portal_id, Genome, portal_link, everything())
write_csv(intable, "JGI_SMC.csv")
