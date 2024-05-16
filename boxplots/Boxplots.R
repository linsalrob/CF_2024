##################################################################
#                                                                #
# Boxplot code                                                   #
#                                                                #
# (c) 2023-2024 Nick Falk                                        #
#                                                                #
# These codes were contributed by Nick                           #
#                                                                #
##################################################################


####STEP 1 Data Wrangling
#set the working directory where the taxa and metadata files are
setwd('')

# Load necessary libraries
library(vegan)  # For OTU data (Example: OTU table)
library(dplyr)  # For data manipulation
library(stats)  # For Kruskal-Wallis test
library(ggplot2)

# Import taxa table
genus_otu <- read.csv(file = "genus_revised_minion_202304.csv", header = TRUE, row.names = 1)

# Visualize the taxa table
#head(genus_otu)

#Normalization/Transformation choice (square-root, log, relative abundance)
#square root transform all the values in the otu table
sqrt_genus <- sqrt(genus_otu)

#log transform all the values in the otu table
#sqrt_genus <- log1p(genus_otu)

#uses decostand within vegan to transfrom the OTU table to relative abundance. "total" is the method for relative abudnance, with "2" meaning to normalize by column
#sqrt_genus_otu <- decostand(genus_otu, "total", 2)

#re-name otu table if no square root is done, but still want to proceed with the code that uses a dataframe called "sqrt_genus_otu"
sqrt_genus_otu <- sqrt_genus

#head(sqrt_genus_otu)

# Transpose the taxa table
genus_otu2 <- t(sqrt_genus_otu)

# Visualize the transposed and transformed taxa table
#head(genus_otu2)

# Import metadata table
metadata <- read.csv(file = "FINAL_METADATA_MINION.csv", header = TRUE, row.names = 1)

# Add the "X" character to the beginning of each rowname in the metadata table to match the taxa table; may or not be applicable
rownames(metadata) <- paste("X", rownames(metadata), sep = "")

# Visualize the revised metadata table
#head(metadata)

# Merge the taxa table and metadata table by the shared rownames
merged1 <- merge(metadata, genus_otu2, by='row.names')

# Visualize merged table
#head(merged1)

# Visualize the structure of the merged table and ensure the categories you require to test are factors and not integers (if 0 and 1's are used, for example)
#str(merged1)

#Convert the columns of interest to a factor (i.e. the 1/0 for culture positive/negative)
merged1$CMS_Pseudomonas.aeruginosa <- factor(merged1$CMS_Pseudomonas.aeruginosa)
merged1$CMS_Stenophotomonas.maltophilia  <- factor(merged1$CMS_Stenophotomonas.maltophilia)
merged1$NTM <- factor(merged1$NTM)
merged1$CMS_Mycobacteroides.abscessus <- factor(merged1$CMS_Mycobacteroides.abscessus)
merged1$CMS_Mycobacterium.intracellulare <- factor(merged1$CMS_Mycobacterium.intracellulare)
merged1$CMS_Staphylococcus..aureus <- factor(merged1$CMS_Staphylococcus..aureus)
merged1$CMS_Achromobacter.xylosoxidans <- factor(merged1$CMS_Achromobacter.xylosoxidans)
merged1$CMS_Burkholderia.cepacia <- factor(merged1$CMS_Burkholderia.cepacia)
merged1$CMS_Haemophilus.influenzae <- factor(merged1$CMS_Haemophilus.influenzae)

####STEP 2 Create Boxplots
#Create an object that will tell the boxplot code to plot only the "0" and "1" levels for the culture categories, and not the "NAs"
specific_levels <- c("0", "1")

create_boxplot <- function(data, culture_taxa, taxa) {
  ggplot(subset(data, get(culture_taxa) %in% specific_levels), aes_string(x = culture_taxa, y = taxa, fill = culture_taxa)) +
    geom_boxplot(outlier.shape = NA, color = "black") +  # Set outline color to black
    geom_jitter(width = 0.25, color = "black") +  # Set color to black for data points
    scale_fill_manual(values = c("0" = "#F8766D", "1" = "#00BFC4")) +  # Set fill colors for factor levels
    theme_minimal() +  # Set background color to white
    ggtitle("") +
    xlab("") +
    ylab("Abundance in Patient") +
    scale_x_discrete(labels = c("0" = "Negative Culture", "1" = "Positive Culture")) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),  # Center the plot title and set font size
      axis.title.x = element_text(size = 14),  # Set font size for x-axis title
      axis.title.y = element_text(size = 14, vjust = 2.7),   # Set font size for y-axis title
      axis.text.x = element_text(size = 12),    # Set font size for x-axis tick labels
      axis.text.y = element_text(size = 12)     # Set font size for y-axis tick labels
    ) +
    scale_y_continuous(labels = scales::number_format(scale = 1)) +
    guides(fill = guide_legend(title = NULL))  # Remove legend title
}

# Example usage:
create_boxplot(merged1, "CMS_Pseudomonas.aeruginosa", "g__Pseudomonas")

####Step 3 Additionally Stats Tests
#perform Kruskal-Wallis Test
kruskal.test(g__Pseudomonas ~ CMS_Pseudomonas.aeruginosa, data = merged1)
# calculate median
g__CMS_Pseudomonas.aeruginosa_median <- aggregate(g__Pseudomonas ~ CMS_Pseudomonas.aeruginosa, merged1, median)
print(g__CMS_Pseudomonas.aeruginosa_median)
#Filter the result to get the median for the specified group; this becomes the threshold value for hits in the "0" grouping
threshold_CMS_Pseudomonas.aeruginosa <- g__CMS_Pseudomonas.aeruginosa_median[g__CMS_Pseudomonas.aeruginosa_median$CMS_Pseudomonas.aeruginosa == 1, "g__Pseudomonas"]
print(threshold_CMS_Pseudomonas.aeruginosa)
#aggregate
aggregate(g__Pseudomonas ~ CMS_Pseudomonas.aeruginosa, merged1, median)
# Find row names where the 'Value' column is above the threshold
above_threshold_rownames_CMS_Pseudomonas.aeruginosa <- merged1 %>%
  filter(CMS_Pseudomonas.aeruginosa == "0" & g__Pseudomonas > threshold_CMS_Pseudomonas.aeruginosa) %>%
  select(Row.names)
# Print the row names
print(above_threshold_rownames_CMS_Pseudomonas.aeruginosa)
write.csv(above_threshold_rownames_CMS_Pseudomonas.aeruginosa, file = "CMS_Pseudomonas.aeruginosa")
#Run pairwise test on groups within a column
pairwise.wilcox.test(merged1$g__Pseudomonas, merged1$CMS_Pseudomonas.aeruginosa, p.adjust.method = "BH")

