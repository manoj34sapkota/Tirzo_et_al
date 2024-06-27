setwd("/Users/biuser/Desktop/Parastou/Manuscript")


df<-read.csv("data.csv")
# Load required libraries
library(dplyr)
library(tidyr)
library(survival)
library(survminer)

# Assuming df is the name of your data frame
# Reformat the data to handle multiple rows for the same case number

reformatted_data <- df %>%
  group_by(Case) %>%
  summarize(
    Age = first(Age),
    Disorder = first(Disorder),
    Sex = first(Sex),
    Results = first(Results),
    Gene..MPN. = first(Gene..MPN.),
    Polymorphism..MPN. = first(Polymorphism..MPN.),
    VAF..MPN. = first(VAF..MPN.),
    Gene = list(Gene),
    Polymorphism = list(Polymorphism),
    VAF = list(VAF),
    Final.Dx=first(Final.Dx)
  ) %>%
  unnest(c(Gene, Polymorphism, VAF))

# Convert VAF column to numeric
reformatted_data$VAF..MPN. <- as.numeric(sub("%", "", reformatted_data$VAF..MPN.))
reformatted_data$VAF <- as.numeric(sub("%", "", reformatted_data$VAF))

# Display the reformatted data
head(reformatted_data)



######### analysis


# Assuming reformatted_data is the name of your reformatted data frame
# Make sure the 'Results' column is a factor for survival analysis
reformatted_data$Results <- as.factor(reformatted_data$Results)
reformatted_data$Age<- as.numeric(reformatted_data$Age)
reformatted_data$Sex<-as.factor(reformatted_data$Sex)
reformatted_data$Disorder<-as.factor(reformatted_data$Disorder)
reformatted_data$Gene..MPN.<-as.factor(reformatted_data$Gene..MPN.)
reformatted_data$Gene<-as.factor(reformatted_data$Gene)
reformatted_data$Final.Dx<-as.factor(reformatted_data$Final.Dx)

# Descriptive Statistics
summary(reformatted_data)

# Additional Visualizations (customize based on your needs)
# Create a boxplot with only unique cases
unique_plot_data <- reformatted_data %>%
  distinct(Case, .keep_all = TRUE) %>%
  mutate(Results = factor(Results, levels = c("NEGATIVE", "POSITIVE")))  # Ensure correct order of factor levels

# Count the number of cases for each result
case_counts <- count(unique_plot_data, Results)

# Create the boxplot
ggplot(unique_plot_data, aes(x = Results, y = Age, fill = Results)) +
  geom_boxplot() +
  geom_text(data = case_counts, aes(x = Results, y = 100, label = paste("n =", n)),
            position = position_dodge(0.8), vjust = -0.5, inherit.aes = FALSE) +
  labs(title = "Boxplot of Age by Results (Unique Cases)",
       x = "Results",
       y = "Age") +
  theme_minimal()

# Assuming reformatted_data is your data frame
result_groups <- split(unique_plot_data$Age, unique_plot_data$Results)

# Perform t-test
t_test_result <- t.test(result_groups$POSITIVE, result_groups$NEGATIVE)

# Display the results
print(t_test_result)

# Example: Barplot for Gender distribution
# Filter out rows with missing values in 'Sex' or 'Results' and select only "M" and "F"

# Create a bar plot
plot <- ggplot(unique_plot_data, aes(x = Results, fill = Sex)) +
  geom_bar(position = "dodge", color = "black") +
  labs(title = "Bar Diagram of Result vs Sex",
       x = "Results",
       y = "Count") +
  theme_minimal()

# Add text labels with the count above each bar
plot +
  geom_text(stat = "count", aes(label = ..count..),
            position = position_dodge(width = 0.9), vjust = -0.5)

table(unique_plot_data$Results, unique_plot_data$Sex)

sex.contingency.table<-table(unique_plot_data$Results, unique_plot_data$Sex)
# Perform the chi-square test
chi_squared_test <- chisq.test(sex.contingency.table)

# Print the test results
print(chi_squared_test)

##################
# Assuming reformatted_data is the name of your reformatted data frame

# Create a contingency table of gene frequencies by Results
gene_table <- table(reformatted_data$Gene, reformatted_data$Results)

# Perform Fisher's exact test with simulated p-value
fisher_test_result <- fisher.test(gene_table, simulate.p.value = TRUE)

# Display the observed frequencies
observed_frequencies <- gene_table
print(observed_frequencies)

# Display the results of Fisher's exact test
print(fisher_test_result)


#In simpler terms, the small p-value (0.0004998) suggests that the observed distribution of genes is unlikely to occur by chance alone, and there is evidence to conclude that the presence of certain genes is associated with the diagnosis.
#Therefore, you may infer that the frequency of certain genes is significantly different between the POSITIVE and NEGATIVE cases. Keep in mind that the actual interpretation may depend on the context and domain knowledge related to your specific dataset.

###################
# Assuming reformatted_data is the name of your reformatted data frame

# Create an empty data frame to store results
gene_association_results <- data.frame(Gene = character(0), P_Value = numeric(0))

# List of unique genes in the dataset
unique_genes <- unique(reformatted_data$Gene)

# Loop through each gene and perform Fisher's exact test
for (gene in unique_genes) {
  # Create a contingency table for the current gene
  gene_table <- table(reformatted_data$Results, reformatted_data$Gene == gene)
  
  # Perform Fisher's exact test with simulated p-value
  fisher_test_result <- fisher.test(gene_table, simulate.p.value = TRUE)
  
  # Extract the p-value
  p_value <- fisher_test_result$p.value
  
  # Append results to the data frame
  gene_association_results <- rbind(gene_association_results, data.frame(Gene = gene, P_Value = p_value))
}

# Display the results
print(gene_association_results)




########## removing empty cells or cases with no additional gene identified
# Filter out cases where the gene information is empty
filtered_data <- reformatted_data[reformatted_data$Gene != "", ]

# Create a contingency table of gene frequencies by Results
gene_table <- table(filtered_data$Gene, filtered_data$Results)

# Perform Fisher's exact test with simulated p-value
fisher_test_result <- fisher.test(gene_table, simulate.p.value = TRUE)

# Display the observed frequencies
observed_frequencies <- gene_table
print(observed_frequencies)

# Display the results of Fisher's exact test
print(fisher_test_result)

# Assuming reformatted_data is the name of your reformatted data frame

# Filter out cases where the gene information is empty
filtered_data <- reformatted_data[reformatted_data$Gene != "", ]

# Create an empty data frame to store results
gene_association_results <- data.frame(Gene = character(0), P_Value = numeric(0))

# List of unique genes in the filtered dataset
unique_genes <- unique(filtered_data$Gene)

# Loop through each gene and perform Fisher's exact test
for (gene in unique_genes) {
  # Create a contingency table for the current gene
  gene_table <- table(filtered_data$Results, filtered_data$Gene == gene)
  
  # Perform Fisher's exact test with simulated p-value
  fisher_test_result <- fisher.test(gene_table, simulate.p.value = TRUE)
  
  # Extract the p-value
  p_value <- fisher_test_result$p.value
  
  # Append results to the data frame
  gene_association_results <- rbind(gene_association_results, data.frame(Gene = gene, P_Value = p_value))
}

# Display the results
print(gene_association_results)



#################
# Assuming gene_association_results is the name of your results data frame
library(ggplot2)

# Set a significance threshold (e.g., 0.05)
significance_threshold <- 0.05

# Create a volcano plot
ggplot(gene_association_results, aes(x = -log10(P_Value), y = Gene, color = factor(P_Value < significance_threshold))) +
  geom_point(size = 3) +
  scale_color_manual(values = c("grey", "red"), guide = "none") +
  labs(title = "Volcano Plot of Gene Association with Diagnosis",
       x = "-log10(P-Value)",
       y = "Gene") +
  theme_minimal()


# Load necessary libraries
library(ggplot2)

# Create a stacked bar plot
ggplot(reformatted_data, aes(x = Gene, fill = Results)) +
  geom_bar(position = "stack") +
  labs(title = "Gene vs. Positive/Negative Cases",
       x = "Gene",
       y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("NEGATIVE" = "red", "POSITIVE" = "green"))


# Create a new dataset with unique cases
unique_cases <- reformatted_data %>%
  distinct(Case, .keep_all = TRUE)

# Encode the presence or absence of additional gene information as a binary variable
unique_cases$Gene_Presence <- ifelse(unique_cases$Gene != "", "Present", "Absent")

# Perform a chi-squared test
contingency_table <- table(unique_cases$Gene_Presence, unique_cases$Results)
chi_squared_test <- chisq.test(contingency_table)

# Print the p-value
print(chi_squared_test)


# Create a bar plot
barplot(contingency_table, beside = TRUE, col = c("lightblue", "lightcoral"),
        main = "Association between Gene Presence and Results",
        xlab = "Gene Presence", ylab = "Frequency")

# Add legend
legend("topright", legend = rownames(contingency_table), fill = c("lightblue", "lightcoral"))

# Add text labels with frequencies
text(x = barplot(as.matrix(contingency_table), beside = TRUE, plot = FALSE),
     y = contingency_table - 20, label = as.vector(contingency_table))

# Add p-value text
text(3, max(contingency_table) -20, paste("p-value =", format(chi_squared_test$p.value, digits = 3)),
     cex = 0.8, col = "darkgreen")



#disorder vs Gene

# Assuming reformatted_data is your data frame
contingency_table_disorder_gene <- table(reformatted_data$Disorder, reformatted_data$Gene)

# Display the contingency table
print(contingency_table_disorder_gene)

write.csv(contingency_table_disorder_gene,"disordervsgeen.csv")


###### focussing only in positive cases

# Assuming reformatted_data is your data frame
library(dplyr)

# Filter rows with positive results
positive_data <- reformatted_data %>% filter(Results == "POSITIVE")

# Count the occurrences of each additional gene in positive cases
gene_counts <- positive_data %>%
  select(Gene) %>%
  table() %>%
  as.data.frame()

# Rename columns
colnames(gene_counts) <- c("Gene", "Count")

# Sort the table by count in descending order
gene_counts <- gene_counts[order(-gene_counts$Count), ]

# Display the summary table
print(gene_counts)


# Assuming gene_counts is your data frame
library(ggplot2)

# Filter out rows with blank gene names
gene_counts <- gene_counts[gene_counts$Gene != "", ]

# Sort the data frame by Count in descending order
gene_counts <- gene_counts[order(-gene_counts$Count), ]

# Plot the bar diagram
ggplot(gene_counts, aes(x = reorder(Gene, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "Gene Counts in Positive Cases",
       x = "Gene",
       y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank())



###### focussing only in negative cases

# Assuming reformatted_data is your data frame
library(dplyr)

# Filter rows with positive results
negative_data <- reformatted_data %>% filter(Results == "NEGATIVE")

# Count the occurrences of each additional gene in positive cases
gene_counts <- negative_data %>%
  select(Gene) %>%
  table() %>%
  as.data.frame()

# Rename columns
colnames(gene_counts) <- c("Gene", "Count")

# Sort the table by count in descending order
gene_counts <- gene_counts[order(-gene_counts$Count), ]

# Display the summary table
print(gene_counts)


# Assuming gene_counts is your data frame
library(ggplot2)

# Filter out rows with blank gene names
gene_counts <- gene_counts[gene_counts$Gene != "", ]

# Sort the data frame by Count in descending order
gene_counts <- gene_counts[order(-gene_counts$Count), ]

# Plot the bar diagram
ggplot(gene_counts, aes(x = reorder(Gene, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "Gene Counts in Negative Cases",
       x = "Gene",
       y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank())




# Assuming reformatted_data is your data frame
library(dplyr)

# Filter positive cases
positive_cases <- reformatted_data[reformatted_data$Results == "POSITIVE", ]
positive_disorder<-positive_cases
positive_disorder$Age<-NULL
positive_disorder[3:10]<-NULL


# Remove duplicate cases and get count summary of disorders
unique_positive_disorders <- positive_disorder %>%
  distinct(Case, .keep_all = TRUE) %>%
  group_by(Disorder) %>%
  summarise(Count = n())

# Print the result
print(unique_positive_disorders)


# Assuming reformatted_data is your data frame
library(dplyr)

# Filter positive cases
negative_cases <- reformatted_data[reformatted_data$Results == "NEGATIVE", ]
negative_disorder<-negative_cases
negative_disorder$Age<-NULL
negative_disorder[3:10]<-NULL


# Remove duplicate cases and get count summary of disorders
unique_negative_disorders <- negative_disorder %>%
  distinct(Case, .keep_all = TRUE) %>%
  group_by(Disorder) %>%
  summarise(Count = n())

# Print the result
print(unique_negative_disorders)


library(ggplot2)

# Filter out disorders with zero count
filtered_positive_disorders <- unique_positive_disorders %>%
  filter(Count > 0, Disorder !="")

# Arrange disorders in descending order of count
filtered_positive_disorders <- filtered_positive_disorders[order(-filtered_positive_disorders$Count), ]

# Plot the bar diagram
ggplot(filtered_positive_disorders, aes(x = reorder(Disorder, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  labs(title = "Count of Disorders in Positive Cases",
       x = "Disorder",
       y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(clip = 'off') +
  geom_text(aes(label = Count), vjust = -0.5)



library(ggplot2)

# Filter out disorders with zero count and empty string
filtered_negative_disorders <- unique_negative_disorders %>%
  filter(Count > 0, Disorder != "")

# Arrange disorders in descending order of count
filtered_negative_disorders <- filtered_negative_disorders[order(-filtered_negative_disorders$Count), ]

# Plot the bar diagram
ggplot(filtered_negative_disorders, aes(x = reorder(Disorder, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "salmon", color = "black") +
  labs(title = "Count of Disorders in Negative Cases",
       x = "Disorder",
       y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(clip = 'off') +
  geom_text(aes(label = Count), vjust = -0.5)


library(ggplot2)

# Filter out disorders with zero count and empty string
filtered_negative_disorders <- unique_negative_disorders %>%
  filter(Count > 0, Disorder != "")

# Arrange disorders in descending order of count
filtered_negative_disorders <- filtered_negative_disorders[order(-filtered_negative_disorders$Count), ]

# Plot the pie chart
ggplot(filtered_negative_disorders, aes(x = "", y = Count, fill = reorder(Disorder, -Count))) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y") +
  labs(title = "Pie Chart of Disorders in Negative Cases",
       fill = "Disorder") +
  theme_minimal()


library(ggplot2)

# Filter out disorders with zero count and empty string
filtered_positive_disorders <- unique_positive_disorders %>%
  filter(Count > 0, Disorder != "")

# Plot the pie chart
ggplot(filtered_positive_disorders, aes(x = "", y = Count, fill = reorder(Disorder, -Count))) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y") +
  labs(title = "Count of Disorders in Positive Cases",
       fill = "Disorder",
       y = "") +
  theme_minimal()



# Assuming reformatted_data is your data frame
# Create a subset with relevant columns
gene_disorder_subset <- reformatted_data %>%
  select(Case, Gene, Disorder, Results) %>%
  filter(Results == "POSITIVE")  # Considering only positive cases

# Create a contingency table
contingency_table <- table(gene_disorder_subset$Gene, gene_disorder_subset$Disorder)

# Assuming contingency_table has been created

# Remove genes with all zero count
contingency_table <- contingency_table[rowSums(contingency_table) > 0, ]

# Remove disorders with all zero count
contingency_table <- contingency_table[, colSums(contingency_table) > 0]

# Remove the gene with no name ("")
contingency_table <- contingency_table[rownames(contingency_table) != "", ]

# Remove the first column with no column name
contingency_table <- contingency_table[, colnames(contingency_table) != ""]

# Display the updated contingency table
print(contingency_table)


# Plot a heatmap of the updated contingency table with rotated y-axis labels
heatmap(contingency_table, 
        col = colorRampPalette(c("white", "blue"))(20),
        #main = "Contingency Table Heatmap",
        ylab = "Gene",
        las = 2,  # Rotate y-axis labels
        cex.axis = 0.7  # Adjust the size of axis labels for better readability
)

# Add a legend
legend("topright", legend = c("Low Count", "High Count"), fill = colorRampPalette(c("white", "blue"))(2), title = "Count")

  # Perform chi-squared test on the updated contingency table
  chi_squared_test <- chisq.test(contingency_table)
  
  # Print the results
  print(chi_squared_test)
  
  # Perform Fisher's exact test
  fisher_test_result <- fisher.test(contingency_table)
  
  # Print the results
  print(fisher_test_result)
  
  
  # List to store chi-squared test results for each gene
  chi_squared_results_genes <- list()
  
  # Perform chi-squared test for each gene
  for (gene in rownames(contingency_table)) {
    gene_contingency_table <- contingency_table[gene, , drop = FALSE]
    
    # Check if the gene has any non-zero counts
    if (sum(gene_contingency_table) > 0) {
      chi_squared_results_genes[[gene]] <- chisq.test(gene_contingency_table)
    } else {
      warning(paste("Insufficient data for chi-squared test for gene", gene))
    }
  }
  
  # Extract p-values from the results
  p_values_genes <- sapply(chi_squared_results_genes, function(x) {
    if (!is.null(x)) {
      return(x$p.value)
    } else {
      return(NA)
    }
  })
  
  # Create a data frame with gene names and corresponding p-values
  gene_p_values_df <- data.frame(Gene = names(p_values_genes), P_Value = p_values_genes)
  
  # Print the results
  print(gene_p_values_df)
  

  # Create a subset with relevant columns
  gene_disorder_subset <- reformatted_data %>%
    select(Case, Gene, Disorder, Results) %>%
    filter(Results == "NEGATIVE")  # Considering only positive cases
  
  # Create a contingency table
  contingency_table <- table(gene_disorder_subset$Gene, gene_disorder_subset$Disorder)
  
  # Assuming contingency_table has been created
  
  # Remove genes with all zero count
  contingency_table <- contingency_table[rowSums(contingency_table) > 0, ]
  
  # Remove disorders with all zero count
  contingency_table <- contingency_table[, colSums(contingency_table) > 0]
  
  # Remove the gene with no name ("")
  contingency_table <- contingency_table[rownames(contingency_table) != "", ]
  
  # Remove the first column with no column name
  contingency_table <- contingency_table[, colnames(contingency_table) != ""]
  
  # Display the updated contingency table
  print(contingency_table)
  # Plot a heatmap of the updated contingency table with rotated y-axis labels and a continuous color scale
  # Display the updated contingency table
  print(contingency_table)
  
  # Plot a heatmap of the updated contingency table with rotated y-axis labels
  heatmap(contingency_table, 
          col = colorRampPalette(c("white", "blue"))(20),
          main = "Contingency Table Heatmap",
          ylab = "Gene",
          las = 2,  # Rotate y-axis labels
          cex.axis = 0.7  # Adjust the size of axis labels for better readability
  )
  
  # Add a legend
  legend("topright", legend = c("Low Count", "High Count"), fill = colorRampPalette(c("white", "blue"))(2), title = "Count")
  
  # Perform Fisher's exact test
  fisher_test_result <- fisher.test(contingency_table)
  
  # Print the results
  print(fisher_test_result)
  
  
  # List to store chi-squared test results for each gene
  chi_squared_results_genes <- list()
  
  # Perform chi-squared test for each gene
  for (gene in rownames(contingency_table)) {
    gene_contingency_table <- contingency_table[gene, , drop = FALSE]
    
    # Check if the gene has any non-zero counts
    if (sum(gene_contingency_table) > 0) {
      chi_squared_results_genes[[gene]] <- fisher.test(gene_contingency_table)
    } else {
      warning(paste("Insufficient data for chi-squared test for gene", gene))
    }
  }
  
  # Extract p-values from the results
  p_values_genes <- sapply(chi_squared_results_genes, function(x) {
    if (!is.null(x)) {
      return(x$p.value)
    } else {
      return(NA)
    }
  })
  
  # Create a data frame with gene names and corresponding p-values
  gene_p_values_df <- data.frame(Gene = names(p_values_genes), P_Value = p_values_genes)
  
  # Print the results
  print(gene_p_values_df)
  

######   comutation plot positive cases
  # Specify the order of gene names
  gene_order <- c("ABL1", "NF1", "KRAS", "RUNX1", "GNAS", "PPM1D", "TP53", "U2AF1", "SF3B1", "SRSF2", "IDH2", "ASXL1", "DNMT3A", "TET2")
  
  # Specify the order of cases
  case_order <- c(
    17, 19, 40, 43, 51, 66, 85, 87, 114, 192, 228, 242, 252, 267, 323, 429, 472, 511, 563, 570,
    579, 636, 668, 670, 716, 719, 734, 798, 855, 869, 935, 987, 995, 1045, 1051, 1071, 1113, 1124,
    1131, 1135, 1149, 1206, 75, 210, 738, 945, 292, 527, 768, 914, 1144, 142, 173, 227, 271, 1053, 1209
  )
  
  call_order <- c(1, 2, 3, 4, 5, 6, 8, 9, 10, 13, 16, 17, 18, 19, 22, 23, 24, 25, 27, 28, 29, 30, 31, 32, 33, 34, 35, 38, 39, 40, 42, 44, 45, 46, 47, 49, 50, 51, 52, 53, 55, 56, 7, 14, 36, 43, 21, 26, 37, 41, 54, 11, 12, 15, 20, 48, 57)
  
  # Filter positive cases and remove rows with empty gene names
  positive_cases <- reformatted_data[reformatted_data$Results == "POSITIVE" & reformatted_data$Gene != "", ]
  
  # Extract unique genes and cases
  unique_genes <- unique(positive_cases$Gene)
  unique_genes <- unique_genes[unique_genes %in% gene_order]  # Filter genes based on specified order
  unique_cases <- unique(positive_cases$Case)
  
  # Create a matrix to store gene presence (1 for present, 0 for absent)
  gene_presence_matrix <- matrix(0, nrow = length(unique_genes), ncol = length(unique_cases))
  
  # Match gene names and cases to matrix indices and set presence to 1
  for (i in 1:nrow(positive_cases)) {
    gene_index <- match(positive_cases$Gene[i], unique_genes)
    case_index <- match(positive_cases$Case[i], unique_cases)
    
    if (!is.na(gene_index) && !is.na(case_index)) {
      gene_presence_matrix[gene_index, case_index] <- 1
    }
  }
  
  # Convert the matrix to a data frame for plotting
  comutation_data_long <- as.data.frame(gene_presence_matrix)
  colnames(comutation_data_long) <- unique_cases
  rownames(comutation_data_long) <- unique_genes
  comutation_data_long$Gene <- rownames(comutation_data_long)
  
  # Reshape the data for plotting
  comutation_data_long <- tidyr::gather(comutation_data_long, key = "Case", value = "Gene_Presence", -Gene)
  
  # Add a column for call number
  comutation_data_long$Call_Number <- 0 + cumsum(c(TRUE, diff(as.numeric(as.factor(comutation_data_long$Case))) != 0))
  
  # Plot the comutation matrix for positive cases with specified gene and case order
  library(ggplot2)
  ggplot(comutation_data_long, aes(x = factor(Call_Number, levels=call_order), y = factor(Gene, levels = gene_order), fill = factor(Gene_Presence))) +
    geom_tile(color = "black", size = 0.1, alpha=0.5) +  # Add borders
    scale_fill_manual(values = c("white", "green")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=6)) +
    labs(x = "Call Number", y = "Gene", fill = "Gene Presence") +
    scale_x_discrete(labels = case_order)
  
  
  
  ######   comutation plot negative cases
  # Filter negative cases and remove rows with empty gene names
  negative_cases <- reformatted_data[reformatted_data$Results == "NEGATIVE" & reformatted_data$Gene != "", ]
  
  # Reverse order for the genes
  gene_order <- c("DNMT3A", "TET2", "ASXL1", "IDH2", "SRSF2", "SF3B1", "U2AF1", "ZRSR2", "DDX41", "PRPF8", "KRAS", "NRAS", "CBL", "SH2B3", "NF1", "SETBP1", "RUNX1", "BCOR", "STAG2", "CEBPA", "ETV6", "SMC1A", "TP53", "PPM1D", "RB1", "MYD88", "CXCR4")
  
  # Reverse the order
  gene_order <- rev(gene_order)
  
  # Extract unique genes and cases
  unique_genes <- unique(negative_cases$Gene)
  unique_genes <- unique_genes[unique_genes %in% gene_order]  # Filter genes based on specified order
  unique_cases <- unique(negative_cases$Case)
  
  # Create a matrix to store gene presence (1 for present, 0 for absent)
  gene_presence_matrix <- matrix(0, nrow = length(unique_genes), ncol = length(unique_cases))
  
  # Match gene names and cases to matrix indices and set presence to 1
  for (i in 1:nrow(negative_cases)) {
    gene_index <- match(negative_cases$Gene[i], unique_genes)
    case_index <- match(negative_cases$Case[i], unique_cases)
    
    if (!is.na(gene_index) && !is.na(case_index)) {
      gene_presence_matrix[gene_index, case_index] <- 1
    }
  }
  
  # Convert the matrix to a data frame for plotting
  comutation_data_long <- as.data.frame(gene_presence_matrix)
  colnames(comutation_data_long) <- unique_cases
  rownames(comutation_data_long) <- unique_genes
  comutation_data_long$Gene <- rownames(comutation_data_long)
  
  # Reshape the data for plotting
  comutation_data_long <- tidyr::gather(comutation_data_long, key = "Case", value = "Gene_Presence", -Gene)
  
  # Add a column for call number
  comutation_data_long$Call_Number <- 0 + cumsum(c(TRUE, diff(as.numeric(as.factor(comutation_data_long$Case))) != 0))
  
  # Plot the comutation matrix for negative cases with specified gene order
  library(ggplot2)
  ggplot(comutation_data_long, aes(x = as.factor(Call_Number), y = factor(Gene, levels = gene_order), fill = factor(Gene_Presence))) +
    geom_tile(color = "black", size = 0.05, alpha=0.7) +  # Add borders
    scale_fill_manual(values = c("white", "red")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=3)) +
    labs(x = "Case Number", y = "Gene", fill = "Gene Presence") +
    scale_x_discrete(labels = unique_cases)
  
  
  
  
  ##PCA
  # Install and load the required packages if not already installed
  library(FactoMineR)
  library(factoextra)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  
  # Assuming your data is already in the "reformatted_data" dataframe
  pca_data <- subreformatted_data[, c("Case", "Gene")]
  
  # Add a column for presence, set to 1
  pca_data$presence <- 1
  
  # Combine rows with the same case by summing up the presence values
  pca_data_combined <- pca_data %>%
    group_by(Case, Gene) %>%
    summarise(presence = sum(presence))
  
  # Spread the data to wide format
  pca_data_wide <- pca_data_combined %>%
    spread(Gene, presence, fill = 0)
  
  
  # Assuming your data is named pca_data
  pca_data_wide <- pca_data %>%
    mutate(presence = 1) %>%
    spread(Gene, presence, fill = 0) %>%
    group_by(Case) %>%
    summarise_all(max) %>%
    column_to_rownames(var = "Case")
  
  
  # Assuming your data is named pca_data_wide
  pca_data_wide <- pca_data_wide %>%
    mutate_at(vars(-Case), ~ifelse(. != 0, 1, 0))
  
  # Print the resulting dataframe
  print(pca_data_wide)
  
  # Assuming your data is named pca_data_wide
  # Add the Final.Dx column to pca_data_wide
  pca_data_wide$Final.Dx <- reformatted_data$Final.Dx[match(pca_data_wide$Case, reformatted_data$Case)]
  
  # Create a data frame with unique Case and Final.Dx combinations
  color_df <- unique(reformatted_data[, c("Case", "Final.Dx")])
  
  # Perform PCA
  pca_result <- prcomp(pca_data_wide[, -c(1, ncol(pca_data_wide))], scale. = TRUE)
  
  # Create a data frame with PCA components and Final.Dx
  pca_df <- data.frame(Case = pca_data_wide$Case, Final.Dx = pca_data_wide$Final.Dx, pca_result$x)
  
  # Merge with color_df to get the color information
  pca_df <- merge(pca_df, color_df, by = "Case")
  
  # Plot PCA with ggplot2
  ggplot(pca_df, aes(x = PC1, y = PC2, color = Final.Dx.x)) +
    geom_point() +
    labs(title = "PCA of Cases Colored by Final.Dx",
         x = "Principal Component 1",
         y = "Principal Component 2") +
    theme_minimal()
  
  # Assuming your data is already in the "reformatted_data" dataframe
  count_data <- subreformatted_data[, c("Case", "Gene", "Final.Dx")]
  
  # Assuming your data is named "count_data"
  # Create a new column indicating gene presence (1 for gene, 0 for no gene)
  count_data$Gene_Presence <- ifelse(count_data$Gene == "", 0, 1)
  
  # Remove duplicates based on Case
  unique_count_data <- unique(count_data[, c("Case", "Gene_Presence", "Final.Dx")])
  
  # Generate a frequency table
  frequency_table <- table(unique_count_data$Gene_Presence, unique_count_data$Final.Dx)
  
  # Print the frequency table
  print(frequency_table)
  
  # Perform Fisher's exact test with simulated p-value
  fisher_result <- fisher.test(frequency_table, simulate.p.value = TRUE)
  
  # Extract simulated p-value
  p_value_fisher <- fisher_result$p.value
  
  # Print the result
  print(p_value_fisher)
  
  
  # Assuming your data is already in the "reformatted_data" dataframe
  raw_count_data <- subreformatted_data[, c("Case", "Gene", "Final.Dx")]
  count_data<-subset(raw_count_data,raw_count_data$Final.Dx!="")
  
  # Assuming your data is named "count_data"
  # Create a new column indicating gene presence (1 for gene, 0 for no gene)
  count_data$Gene_Presence <- ifelse(count_data$Gene == "", 0, 1)
  
  # Remove duplicates based on Case
  unique_count_data <- unique(count_data[, c("Case", "Gene_Presence", "Final.Dx")])
  
  # Generate a frequency table
  frequency_table <- table(unique_count_data$Gene_Presence, unique_count_data$Final.Dx)
  
  
  
  # Assuming your data is named "gene_presence_data"
  # Transpose the data to have genes as rows and factors of Final.Dx as columns
  gene_presence_data<-read.csv("finaldx.csv", row.names = 1)
  
  # Assuming that gene_presence_data is your data frame
  
  # Transpose the data frame
  gene_presence_data_t <- as.data.frame(t(gene_presence_data))
  gene_presence_data_t$Gene <- row.names(gene_presence_data_t)
  
  # Install and load necessary libraries
  library(ggplot2)
  library(reshape2)  # for melt function
  
  # Melt the data frame for ggplot
  melted_data <- melt(gene_presence_data_t, id.vars = "Gene", variable.name = "Final.Dx", value.name = "Count")
  
  # Create a heatmap using ggplot2
  ggplot(melted_data, aes(x = Gene, y = Final.Dx, fill = Count)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "white", high = "red") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Gene Presence Heatmap",
         x = "Gene",
         y = "Final.Dx")
  
  
  