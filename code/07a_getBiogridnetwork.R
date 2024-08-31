library(tidyverse)

# https://downloads.thebiogrid.org/BioGRID/Release-Archive/BIOGRID-4.4.235/
# just used tab 3 cause biogrid recommend to use that on their site:

# download into wd from: https://downloads.thebiogrid.org/File/BioGRID/Release-Archive/BIOGRID-4.4.235/BIOGRID-ALL-4.4.235.tab3.zip


bioGridAll <- read.delim("BIOGRID-ALL-4.4.235.tab3.txt", header = T, sep = "\t")


##human data and select 2 columns (gene symbols for each interactor)
biogridAll_human <- bioGridAll %>% filter(Organism.ID.Interactor.A == "9606" &
                                            Organism.ID.Interactor.B == "9606" &
                                            Organism.Name.Interactor.A == "Homo sapiens" &
                                            Organism.Name.Interactor.B == "Homo sapiens")%>%
  dplyr::select(Official.Symbol.Interactor.A, Official.Symbol.Interactor.B)


write.table(biogridAll_human, "raw/allNetworkData.txt", row.names = F, quote = F, sep = "\t")
