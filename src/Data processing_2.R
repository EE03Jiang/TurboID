install.packages("lintr")
install.packages("styler")

library(lintr)
library(styler)
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)


#Set the working directory for the project
all_proteins <- read_excel("data/raw/20230404_SA_HEK293T_TurboID_construct7_Biotin_R1.xlsx")
dim(all_proteins)
colnames(all_proteins)

#find the quantitative columns
find_abundance <- function(df) {
  #Initialize the vector of abundance columns
  abundance_columns <- character()
  for ( cname in colnames(df) ) {
    #check the columns one by one
    if ( grepl("construct7", cname) ) {
      abundance_columns <- c(abundance_columns, cname)
    }
  }
  abundance_columns
}
quan_columns <- find_abundance(all_proteins)
quan_columns

dfWide <- all_proteins %>%
  subset (select=c("Accession", quan_columns) ) %>%
  na.omit()

##Log2 transformation
rownames(dfWide) <- dfWide$Accession
dfWide$Accession <- NULL
dfWide <- log2(dfWide)
summary(dfWide)


#Box Plot
boxplot(
  Log2_Abund~Sample, data = gather(dfWide, Sample, Log2_Abund),
  main = "Original Log2 Ratios"
)


##Normalization
library(dplyr)
library(ggplot2)

#For each column, subtract the median of the column from each of it's values
dfNorm <- mapply('-', dfWide, apply(dfWide,2,median))
#Transform into a dataframe
dfNorm <- as.data.frame(dfNorm, row.names = row.names(dfWide))


# Calculate the median for the "construct7" columns
construct7_medians <- c(median(dfWide$`construct7_R1`, na.rm = TRUE),
                               median(dfWide$`construct7_R2`, na.rm = TRUE),
                               median(dfWide$`construct7_R3`, na.rm = TRUE))
construct7_avg_median <- mean(construct7_medians)

# Calculate the average of medians for "construct7+Biotin"
construct7_Biotin_medians <- c(median(dfWide$`construct7+Biotin_R1`, na.rm = TRUE),
                               median(dfWide$`construct7+Biotin_R2`, na.rm = TRUE),
                               median(dfWide$`construct7+Biotin_R3`, na.rm = TRUE))
construct7_Biotin_avg_median <- mean(construct7_Biotin_medians)

# Subtract the difference between the average of medians of "construct7+Biotin_R1", "construct7+Biotin_R2", "construct7+Biotin_R3"
# and the average of medians of "construct7_R1", "construct7_R2", "construct7_R3"
avg_median_diff <- construct7_Biotin_avg_median - construct7_avg_median
dfNorm$`construct7_R1` <- dfNorm$`construct7_R1` - avg_median_diff
dfNorm$`construct7_R2` <- dfNorm$`construct7_R2` - avg_median_diff
dfNorm$`construct7_R3` <- dfNorm$`construct7_R3` - avg_median_diff

# Transform into a dataframe
dfNorm <- as.data.frame(dfNorm, row.names = row.names(dfWide))

# Draw the boxplot
boxplot(Log2_Abund ~ Sample, data = gather(dfNorm, Sample, Log2_Abund),
        main = "Normalized Log2 Ratios")


## Imputation
