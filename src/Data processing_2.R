if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DEqMS")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("fgsea")

library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DEqMS)
library(fgsea)

#Set the working directory for the project
setwd("C:/Users/Yihuei_Yiang/Documents/TurboID/data/raw")
all_proteins <- read_excel("R-test_20230404_SA_HEK293T_TurboID_construct7_Biotin_R1.xlsx")
dim(all_proteins)

colnames(all_proteins)
dim(all_proteins)

#Add gene names from the description
#all_proteins$Gene <- apply(
#  all_proteins, 1, function(x) {
#    descr <- x["Description"]
#    #OBS: format-dependent operations
#    gn <- strsplit(descr, "GN=")[[1]][[2]]
#    strsplit(gn, " ")[[1]][[1]]
#  }
#)
#dim(all_proteins)

#Construct a function to rename the abundance
#rename_ratios <- function(df) {
#  for (c in colnames(df)) {
#    if (grepl("Abundance:", c)) {
#      new_colname <- gsub("Abundance:", "", c)
#      new_colname <- gsub(": Sample", "", new_colname)
#      new_colname <- gsub("F1", "Contruct 7 + Biotin_R1", new_colname)
#      new_colname <- gsub("F2", "Contruct 7 + Biotin_R2", new_colname)
#      new_colname <- gsub("F3", "Contruct 7 + Biotin_R3", new_colname)
#      new_colname <- gsub("F4", "Construct 7_R1", new_colname)
#      new_colname <- gsub("F5", "Construct 7_R2", new_colname)
#      new_colname <- gsub("F6", "Construct 7_R3", new_colname)
#      names(df)[names(df) == c] <- new_colname
#    }
#  }
#  df
#}


#res <- rename_ratios(all_proteins)
#all_proteins <- res[[1]]
#quan_columns <- res[[2]]
#quan_columns


#If you only want to find the quantitative columns, but not to rename them

find_ratios <- function(df) {
  #Initialize the vector of ratio columns
  ratio_columns <- character()
  for ( cname in colnames(df) ) {
    #check the columns one by one
    if ( grepl("construct7", cname) ) {
      ratio_columns <- c(ratio_columns, cname)
    }
  }
  ratio_columns
}
quan_columns <- find_ratios(all_proteins)
quan_columns

dfWide <- all_proteins %>%
  filter(!grepl("cont_",Accession)) %>%
  subset (select=c("Accession", quan_columns) ) %>%
  na.omit()

rownames(dfWide) <- dfWide$Accession
dfWide$Accession <- NULL
dfWide <- log2(dfWide)
#Look at the distribution of quan values
summary(dfWide)

#Box Plot
boxplot(
  Log2_Abund~Sample, data = gather(dfWide, Sample, Log2_Abund),
  main = "Original Log2 Ratios"
)


library(tidyr)
library(dplyr)
library(ggplot2)

#For each column, subtract the median of the column from each of it's values
dfNorm <- mapply('-', dfWide, apply(dfWide,2,median))
#Transform into a dataframe
dfNorm <- as.data.frame(dfNorm, row.names = row.names(dfWide))


# Calculate the median for the "construct7_R1", "construct7_R2", and "construct7_R3" columns
construct7_R1_median <- median(dfWide$`construct7_R1`, na.rm = TRUE)
construct7_R2_median <- median(dfWide$`construct7_R2`, na.rm = TRUE)
construct7_R3_median <- median(dfWide$`construct7_R3`, na.rm = TRUE)

# Calculate the average of medians for "construct7+Biotin_R1", "construct7+Biotin_R2", "construct7+Biotin_R3"
construct7_Biotin_medians <- c(median(dfWide$`construct7+Biotin_R1`, na.rm = TRUE),
                               median(dfWide$`construct7+Biotin_R2`, na.rm = TRUE),
                               median(dfWide$`construct7+Biotin_R3`, na.rm = TRUE))
construct7_Biotin_avg_median <- mean(construct7_Biotin_medians)

# Subtract the median of "construct7_R1", "construct7_R2", and "construct7_R3"
dfNorm$`construct7_R1` <- dfWide$`construct7_R1` - construct7_R1_median
dfNorm$`construct7_R2` <- dfWide$`construct7_R2` - construct7_R2_median
dfNorm$`construct7_R3` <- dfWide$`construct7_R3` - construct7_R3_median

# Subtract the difference between the average of medians of "construct7+Biotin_R1", "construct7+Biotin_R2", "construct7+Biotin_R3"
# and the average of medians of "construct7_R1", "construct7_R2", "construct7_R3"
avg_median_diff <- construct7_Biotin_avg_median - mean(c(construct7_R1_median, construct7_R2_median, construct7_R3_median))
dfNorm$`construct7_R1` <- dfNorm$`construct7_R1` - avg_median_diff
dfNorm$`construct7_R2` <- dfNorm$`construct7_R2` - avg_median_diff
dfNorm$`construct7_R3` <- dfNorm$`construct7_R3` - avg_median_diff

# Transform into a dataframe
dfNorm <- as.data.frame(dfNorm, row.names = row.names(dfWide))

# Draw the boxplot
boxplot(Log2_Abund ~ Sample, data = gather(dfNorm, Sample, Log2_Abund),
        main = "Normalized Log2 Ratios")



## Data filtering function
filter_valids = function(df, conditions, min_count, at_least_one = FALSE) {
  # df = data frame containing LOG2 data for filtering and organized by data type
  # conditions = a character vector dictating the grouping
  # min_count = a numeric vector of the same length as "conditions" indicating the minimum 
  #     number of valid values for each condition for retention
  # at_least_one = TRUE means to keep the row if min_count is met for at least one condition
  #     FALSE means min_count must be met across all conditions for retention
  
  log2.names = grep("construct7", names(df), value = TRUE)   # Extract LOG2 column names
  cond.names = lapply(conditions, # Group column names by conditions
                      function(x) grep(x, log2.names, value = TRUE, perl = TRUE))
  
  cond.filter = sapply(1:length(cond.names), function(i) {
    df2 = df[cond.names[[i]]]   # Extract columns of interest
    df2 = as.matrix(df2)   # Cast as matrix for the following command
    sums = rowSums(is.finite(df2)) # count the number of valid values for each condition
    sums >= min_count[i]   # Calculates whether min_count requirement is met
  })
  
  if (at_least_one) {
    df$KEEP = apply(cond.filter, 1, any)
  } else {
    df$KEEP = apply(cond.filter, 1, all)
  }
  
  return(df)  # No rows are omitted, filter rules are listed in the KEEP column
}


## Apply filtering
df.F = filter_valids(df,
                     conditions = c("Parental", "Resistant"),
                     min_count = c(2, 2),
                     at_least_one = TRUE)

