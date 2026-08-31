# 5his code is intended to run the related package. testing on R version 4.3.3
# have to run ggplot2 for the related package to work
library(ggplot2)

# install related package. this can be done from https://www.frasierlab.ca/software
install.packages(
  "C:/Users/jgard/Downloads/related_1.0.tar.gz",
  repos = NULL,
  type = "source"
)

# load package
library(related)


# Load data --------------------------------------------


## Groupings ------------------------------------------------
# related can handle groups, which here are populations. We do expect some admixture between populations, as well.
# To do that, I put an ID at the beginning of each population and year. 
## Big Creek = BC
## Scott Creek = SC
## Carmel River = CR
## russian river = RR
## Napa River= NR

# For year i used:
## 2024: A
## 2025: B

# Ancestry -------------------------------------------------