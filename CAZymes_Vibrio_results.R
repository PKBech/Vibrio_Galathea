#Summaries hmm results from hmm dbcan search'

library(tidyverse)
library(dplyr)
library(stringr)




# Create an empty data frame to store the results
Vibrio_CAZymes_Summery <- data.frame(CAZyme_fam = c("AA","CBM","CE","GH","GT","PL","cohesin"))
Vibrio_CAZymes_table <- data.frame(CAZyme = c("AA1","AA10"))


# List all files with the pattern "S*.out.dm.txt"
file_list <- list.files(path = "../hmmer/hmms_Vibrios/", pattern = "*.txt", full.names = TRUE)

# Loop through the list of files
for (file in file_list) {
  
  # Read the file
  df <- read.table(file = file, header = F, sep = "\t")
  
  # Rename the columns
  colnames(df) <- c("target name","accession","tlen", "query name","accession","qlen","E-value","score",
                         "bias", "no.","of","c-Evalue","i-Evalue","score", "bias","from","to", "from", 
                         "to","from","to","acc")
  
  
  # Summarize the data
  df_sum <- df %>% 
    select(`query name`, `target name`, `E-value`) %>% 
    group_by(`query name`) %>% 
    slice(which.min(`E-value`)) %>% 
    rename(AA_ID = `query name`, CAZyme = `target name`, eval = `E-value`) %>% 
    mutate(CAZyme = gsub("\\.hmm", "", CAZyme), CAZyme_fam = str_extract(CAZyme, "[A-Za-z]+")) %>% 
    group_by(CAZyme_fam) %>% 
    dplyr::summarise(summary = n())
  
  #Create full table with all CAZymes subfamilies
  df_table <- df %>% 
    select(`query name`, `target name`, `E-value`) %>% 
    group_by(`query name`) %>% 
    slice(which.min(`E-value`)) %>% 
    rename(AA_ID = `query name`, CAZyme = `target name`, eval = `E-value`) %>% 
    mutate(CAZyme = gsub("\\.hmm", "", CAZyme), CAZyme_fam = str_extract(CAZyme, "[A-Za-z]+")) %>% 
    group_by(CAZyme) %>% 
    dplyr::summarise(Strain = n())
  
  # Extract the first five characters of the file name
  file_prefix <- str_sub(basename(file), 1, 5)
  
  # Rename the second column with the file prefix
  colnames(df_sum)[2] <- paste0(file_prefix, "_", colnames(df_sum)[2])
  colnames(df_table)[2] <- paste0(file_prefix, "_", colnames(df_table)[2])
  
  
  # Full join the summarized data to the result data frame
  Vibrio_CAZymes_Summery <- full_join(Vibrio_CAZymes_Summery, df_sum, by = "CAZyme_fam")
  Vibrio_CAZymes_table <- full_join(Vibrio_CAZymes_table, df_table, by = "CAZyme")
}

# Save the final result data frame as csv files
str(Vibrio_CAZymes_Summery)
str(Vibrio_CAZymes_table)

#Replace NAs with zeros
Vibrio_CAZymes_table[is.na(Vibrio_CAZymes_table)] <- 0
Vibrio_CAZymes_Summery[is.na(Vibrio_CAZymes_Summery)] <- 0

# Save an object to a file
saveRDS(Vibrio_CAZymes_table, file = "Vibrio_CAZymes_table.rds")
saveRDS(Vibrio_CAZymes_Summery, file = "Vibrio_CAZymes_Summery.rds")


rownames(Vibrio_CAZymes_table) <- Vibrio_CAZymes_table$CAZyme
Vibrio_CAZymes_table <- Vibrio_CAZymes_table[,-1] 
Vibrio_CAZymes_table_t <- t(Vibrio_CAZymes_table)

rownames(Vibrio_CAZymes_Summery) <- Vibrio_CAZymes_Summery$CAZyme_fam
Vibrio_CAZymes_Summery <- Vibrio_CAZymes_Summery[,-1] 
Vibrio_CAZymes_Summery_t <- t(Vibrio_CAZymes_Summery)

colnames(Vibrio_CAZymes_table_t) <- paste0("Hmm", "!", colnames(Vibrio_CAZymes_table_t))
colnames(Vibrio_CAZymes_Summery_t) <- paste0("Hmm_Big", "!", colnames(Vibrio_CAZymes_Summery_t))


#write.csv(Vibrio_CAZymes_table_t, "Vibrio_CAZymes_table.csv", row.names = TRUE, quote = FALSE)
#write.csv(Vibrio_CAZymes_Summery_t, "Vibrio_CAZymes_Summery.csv", row.names = TRUE, quote = FALSE)


