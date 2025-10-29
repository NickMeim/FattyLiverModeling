library(tidyverse)
tf_responses <- readRDS("../data/ChemPert/Transcriptional_responses.rds")
colnames(tf_responses)[1] <- "Response_ID"
metadata <- read.csv("../data/ChemPert/Information_for_transcriptional_responses.csv")

metadata_human <- metadata %>% filter(Species == "Human" & grepl("Hepatocy", Cell_Source))
tf_responses_hep <- tf_responses %>% filter(Response_ID %in% metadata_human$Response_ID)

resp_net <- NULL
for (ii in 1:nrow(tf_responses_hep)){
  target <- strsplit(tf_responses_hep$`Up-regulation`[ii], 
                     split = "; ", fixed = T)  %>% unlist() %>% unique()
  if (length(target) > 0){
    dummy <- data.frame(source = tf_responses_hep$Response_ID[ii],
                        target = target,
                        mor = 1)
  } else {
    dummy <- NULL
  }
  
  target <- strsplit(tf_responses_hep$`Down-regulation`[ii], 
                     split = "; ", fixed = T)  %>% unlist() %>% unique()
  if (length(target) > 0){
    dummy <- rbind(dummy,
                   data.frame(source = tf_responses_hep$Response_ID[ii],
                              target = target,
                              mor = -1))
  }
  
  resp_net <- rbind(resp_net, dummy)
}

resp_net <- resp_net %>% unique() %>% mutate(Response_ID = source)
resp_net <- merge(resp_net, metadata_human, by = "Response_ID")

save(resp_net, file = "chempert_response_network.RData")