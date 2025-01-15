### This code replicates all the analysis and figures done in the paper entitled: 
# "Comparing methods for estimating geographic ranges in freshwater fishes"
# by Valencia-Rodríguez et al. Submitted to Freshwater Biology

#-----------------------------------------------------------------------------#
####             3.1 Comparison of range sizes between methods             ####
#-----------------------------------------------------------------------------#

## This analysis only considers ranges restricted to water bodies for each method

# Load libraries
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# Set your working directory
setwd("~/yourpath")

## Load the "ranges_size.csv" table with the range size estimates for each method. 
# This file is available in section 2.6 or in the tables folder. 
data <- read_csv("ranges_size.csv")

### Prepare the data for analysis ###

# Select only the columns with methods to evaluate range sizes restricted to water bodies
ranges <- data %>% select(Convex_network,Static_network,Dynamic_network,SDM,Expert_network) %>% 
  rename(Convex=Convex_network, Static=Static_network, Dynamic=Dynamic_network, Expert=Expert_network)

## Kolmogorov-Smirnov (KS) test to compare range sizes estimated by each pair of methods

# Function to perform KS test between pairs of methods
ks_test <- function(col1, col2, data) {
  filtered_data <- data %>% select(all_of(col1), all_of(col2)) %>% drop_na()
  test_result <- ks.test(filtered_data[[col1]], filtered_data[[col2]])
  return(c(col1, col2, test_result$statistic, test_result$p.value, nrow(filtered_data)))
}

# Create a list of all possible method pairs
method_pairs <- combn(names(ranges), 2, simplify = FALSE)
# KS test for each pair of methods
ks_results <- lapply(method_pairs, function(pair) {
  ks_test(pair[1], pair[2], ranges)
})

# Merge the results into an array
ks_results_df <- do.call(rbind, ks_results)
colnames(ks_results_df) <- c("Method1", "Method2", "KS_Statistic", "p_KS", "Sample_Size")

# Convert results to a data frame and ensure correct data types
ks_results_df <- as.data.frame(ks_results_df)
ks_results_df$KS_Statistic <- as.numeric(as.character(ks_results_df$KS_Statistic))
ks_results_df$p_KS <- as.numeric(as.character(ks_results_df$p_KS))
ks_results_df$Sample_Size <- as.numeric(as.character(ks_results_df$Sample_Size))
# Print the results
print(ks_results_df)

## Note: The "ks_results_df" object will be used after constructing FIGURE 3 to generate Tables S3 ##

### FIGURE 3 ###

# Convex hull
hull <- ggplot(ranges, aes(Convex/100000)) +
  geom_histogram(fill="grey80", color="black", bins=20) +
  labs(x = expression("Species range size"), y="Number of species") +
  scale_y_continuous(limits=c(0, 45), breaks=seq(0, 45, by=5), expand=c(0,0)) +
  scale_x_continuous(breaks=seq(0, 12, by=2), expand=c(0,0))+
  theme_classic() +
  theme(text = element_text(size=18),
        axis.line=element_line(color="black"),
        axis.text=element_text(color="black", size=16),
        axis.title=element_text(color="black"),
        plot.title=element_text(hjust=0.5, size=16))+
  ggtitle("Convex hull")

# Static alpha
static <- ggplot(ranges, aes(Static/100000)) +
  geom_histogram(fill="grey80", color="black", bins=20) +
  labs(x = expression("Species range size"), y="") +
  scale_y_continuous(limits=c(0, 45), breaks=seq(0, 45, by=5), expand=c(0,0)) +
  scale_x_continuous(breaks=seq(0, 12, by=2), expand=c(0,0)) +
  theme_classic() +
  theme(text=element_text(size=18),
        axis.line=element_line(color="black"),
        axis.text=element_text(color="black", size=16),
        axis.title=element_text(color="black"),
        plot.title=element_text(hjust=0.5, size=16))+
  ggtitle("Static alpha")

# Dynamic alpha
dynamic <- ggplot(ranges, aes(Dynamic/100000)) +
  geom_histogram(fill="grey80", color="black", bins=20) +
  labs(x = expression("Species range size"), y="") +
  scale_y_continuous(limits=c(0, 45), breaks=seq(0, 45, by=5), expand=c(0,0)) +
  scale_x_continuous(breaks=seq(0, 12, by=2), expand=c(0,0)) +
  theme_classic() +
  theme(text=element_text(size=18),
        axis.line=element_line(color="black"),
        axis.text=element_text(color="black", size=16),
        axis.title=element_text(color="black"),
        plot.title=element_text(hjust=0.5, size=16))+
  ggtitle("Dynamic alpha")

# Expert maps
expert <- ggplot(ranges, aes(Expert/100000)) +
  geom_histogram(fill="grey80", color="black", bins=20) +
  labs(x = expression("Species range size"), y="") +
  scale_y_continuous(limits=c(0, 45), breaks=seq(0, 45, by=5), expand=c(0,0)) +
  scale_x_continuous(breaks=seq(0, 12, by=1), expand=c(0,0)) +
  theme_classic() +
  theme(text = element_text(size = 18),
        axis.line=element_line(color="black"),
        axis.text=element_text(color="black", size=16),
        axis.title=element_text(color="black"),
        plot.title=element_text(hjust=0.5, size=16))+
  ggtitle("Expert map")

# Species Distribution Model (SDM)
sdm <- ggplot(ranges, aes(SDM/100000)) +
  geom_histogram(fill="grey80", color="black", bins=20) +
  labs(x = expression("Species range size"), y="") +
  scale_y_continuous(limits=c(0, 45), breaks=seq(0, 45, by=5), expand=c(0,0)) +
  scale_x_continuous(breaks=seq(0, 12, by=2), expand=c(0,0)) +
  theme_classic() +
  theme(text=element_text(size = 18),
        axis.line=element_line(color="black"),
        axis.text=element_text(color="black", size=16),
        axis.title=element_text(color="black"),
        plot.title=element_text(hjust=0.5, size=16))+
  ggtitle("Species Distribution Model")

# Combine images into a single figure and export
group.fig.3 <- wrap_elements(hull) | wrap_elements(static) | wrap_elements(dynamic) | 
  wrap_elements(expert) | wrap_elements(sdm)

Fig.3 <- group.fig.3 + theme(plot.margin = margin(1, 1, 1, 1)) # Adjust margins

# Save figure
ggsave(Fig.3, filename="Fig_3.tiff", width=18, height=7,
       units="cm", scale=2, dpi=400)


### Table S3 ### 

## Prepare the data for analysis
# Convert the "ranges" object to long format
data_long <- pivot_longer(ranges,cols=c(Convex,Static,Dynamic,SDM,Expert),
                          names_to="Method", values_to="RangeSize") %>%
  drop_na() # Remove rows with NA

# Kruskal-Wallis test to test differences between range sizes
kruskal_test <- kruskal.test(RangeSize ~ Method, data = data_long)

# Function to perform Wilcoxon test between all possible pairs of methods
wilcox_test <- function(col1, col2, ranges) {
  filtered_data <- ranges %>% select(all_of(col1), all_of(col2)) %>% drop_na()
  test_result <- wilcox.test(filtered_data[[col1]], filtered_data[[col2]])
  return(c(col1, col2, test_result$statistic, test_result$p.value))
}

# Create a list of all method pairs
method_pairs <- combn(names(ranges), 2, simplify = FALSE)

# Wilcoxon test for each pair of methods
wilcox_results <- lapply(method_pairs, function(pair) {
  wilcox_test(pair[1], pair[2], ranges)
})

# Combine results into a matrix
# Using do.call to combine Wilcoxon test results into a matrix and then convert to data.frame
wilcox_results_df <- do.call(rbind, wilcox_results)
colnames(wilcox_results_df) <- c("Method1", "Method2", "W", "p_Wilcoxon")

# Ensure numerical columns are of the correct type
wilcox_results_df <- as.data.frame(wilcox_results_df)
wilcox_results_df$W <- as.numeric(wilcox_results_df$W)
wilcox_results_df$p_Wilcoxon <- as.numeric(wilcox_results_df$p_Wilcoxon)

# Merge tables using Method1 and Method2 as key columns
## Note: The "ks_results_df" object was constructed in the first part of this section
merged_df <- merge(ks_results_df, wilcox_results_df, by = c("Method1", "Method2"))

# Add column with compared methods
merged_df$Comparative_Method <- paste(merged_df$Method1, "-", merged_df$Method2)

# Rearrange columns and rename for the final format of Table S3
final_df <- merged_df[,c("Comparative_Method","W","p_Wilcoxon",
                         "KS_Statistic","p_KS", "Sample_Size")]
colnames(final_df) <- c("Comparative Method","W","p(W)","D","p(D)","Sample size")

# Format p-values
final_df$`p(W)` <- ifelse(final_df$`p(W)` < 0.05, paste0(formatC(final_df$`p(W)`, format="f", digits=3), "*"), formatC(final_df$`p(W)`, format="f", digits=2))
final_df$`p(D)` <- ifelse(final_df$`p(D)` < 0.05, paste0(formatC(final_df$`p(D)`, format="f", digits=3), "*"), formatC(final_df$`p(D)`, format="f", digits=2))

# Display the result
print(final_df)

# Save the results as Table S3
write.csv(final_df, "Table_S3.csv", row.names = FALSE)

#-----------------------------------------------------------------------------#
#### 3.2 Unrestricted vs. restricted range size within each polygon method ####
#-----------------------------------------------------------------------------#

## This analysis only considers ranges derived from each polygon method

# Load libraries
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

### Load and prepare the data ###
# Load the table with estimated range sizes for each method
data <- read_csv("ranges_size.csv") # This is the same table from the previous analysis

### Transform the table ###
# Reshape the data from a wide format to a long format, combining method and shape information
# into separate columns and ensuring that each row represents a single observation of range size.
df <- data %>%
  pivot_longer(cols= -c(Species, 'Number of occurrence records', 'Alpha value'), 
               names_to=c("Method_Shape"),values_to="Value") %>%
  separate(Method_Shape, into=c("Method", "Shape"), sep="_", fill="right") %>%
  mutate(Method=recode(Method, Convex="Convex hull", Dynamic="Dynamic alpha",
                       Static="Static alpha", Expert="Expert map"),
    Shape = ifelse(is.na(Shape), "Polygon", recode(Shape, network="Water bodies", netw="Water bodies"))) %>% 
  select(Method, Shape, Value) %>%
  drop_na(Value)

# Kruskal-Wallis test to detect differences within methods
kruskal <- kruskal.test(df$Value ~ df$Shape, data = df)

# Filter the data to remove NAs and exclude SDM, which only has restricted ranges in water bodies
df.wilcoxon <- df %>% filter(!is.na(Shape) & Method != "SDM")

# List of methods
methods <- unique(df.wilcoxon$Method)
# List to store Wilcoxon test results within each method
results_list <- list()

for (method in methods) {
  df_method <- df.wilcoxon %>% filter(Method == method)
  if (length(unique(df_method$Shape)) == 2) {
    test <- wilcox.test(Value ~ Shape, data = df_method, paired = FALSE)
    results_list[[method]] <- list(p.value = test$p.value, W = test$statistic)
  } else {
    results_list[[method]] <- list(p.value = NA, W = NA)
  }
}

# Create dataframe with the results
results_df <- do.call(rbind, lapply(names(results_list), function(method) {
  data.frame(Method=method, P_Value=results_list[[method]]$p.value, W=results_list[[method]]$W)}))


### FIGURE 4 ###

orden_nom <- c("Convex hull", "Static alpha", "Dynamic alpha", "Expert map", "SDM")
df$Method <- factor(df$Method, levels = orden_nom)
df$Shape <- factor(df$Shape, levels = c("Polygon", "Water bodies"))

# Filter out SDM rows
df_filtered <- df %>% filter(Method != "SDM")

Fig.4 <- 
  ggplot(df_filtered, aes(x = Method, y = Value / 100000, fill = Shape)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#E69F00", "#999999"))+
  labs(x = "", y = expression("Range size (" * km^2 * " x " *10^5* ")")) +
  scale_y_continuous(limits = c(0, 100), breaks=seq(0, 100, by=10), expand=c(0,0))+
  theme_classic() +
  theme(text = element_text(size = 18, color = "black"),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black", size = 16))

# Save figure
ggsave(Fig.4, filename = "Fig_4.tiff", width = 18, height = 7,
       units ="cm", scale = 2, dpi = 400)

#-----------------------------------------------------------------------------#
####      3.3 Differences and similarities among species' range sizes      ####
#-----------------------------------------------------------------------------#

## This analysis only considers ranges restricted to water bodies for each method

# Load libraries
library(dplyr)
library(ggplot2)
library(gridExtra)

## Apply standard regression models to compare differences in range sizes estimated by each method

## Note: We call the data table adjusted in section 3.1 (object "ranges"),
# which contains the range sizes of each method adjusted to water bodies

# Define the methods to compare
methods <- c("Convex", "Static", "Dynamic", "SDM", "Expert")

# Function to perform linear regressions between pairs of methods
lm_test <- function(col1, col2, data) {
  formula <- as.formula(paste(col1, "~", col2))
  lm_result <- lm(formula, data = data)
  summary_result <- summary(lm_result)
  
  intercept <- summary_result$coefficients[1, 1]
  intercept_t <- summary_result$coefficients[1, 3]
  intercept_p <- summary_result$coefficients[1, 4]
  
  slope <- summary_result$coefficients[2, 1]
  slope_t <- summary_result$coefficients[2, 3]
  slope_p <- summary_result$coefficients[2, 4]
  
  r_squared <- summary_result$r.squared
  
  return(c(col1, col2, intercept, intercept_t, intercept_p, slope, slope_t, slope_p, r_squared))
}

# Create a list of all possible pairs of methods
method_pairs <- combn(methods, 2, simplify = FALSE)

# Perform linear regressions for each pair of methods
lm_results <- lapply(method_pairs, function(pair) {
  lm_test(pair[1], pair[2], ranges)
})

# Combining the rows of the list into a matrix and assign column names to the resulting matrix
lm_results_df <- do.call(rbind, lm_results)
colnames(lm_results_df) <- c("Comparative Method", "Method2", "Intercept", "Intercept_t", 
                             "Intercept_p", "Slope", "Slope_t", "Slope_p", "R-squared")

# Convert the matrix to a data frame
lm_results_df <- as.data.frame(lm_results_df)

# Verify and convert data types in the resulting data frame to ensure 
# numerical values are treated correctly
lm_results_df$Intercept <- as.numeric(as.character(lm_results_df$Intercept))
lm_results_df$Intercept_t <- as.numeric(as.character(lm_results_df$Intercept_t))
lm_results_df$Intercept_p <- as.numeric(as.character(lm_results_df$Intercept_p))
lm_results_df$Slope <- as.numeric(as.character(lm_results_df$Slope))
lm_results_df$Slope_t <- as.numeric(as.character(lm_results_df$Slope_t))
lm_results_df$Slope_p <- as.numeric(as.character(lm_results_df$Slope_p))
lm_results_df$`R-squared` <- as.numeric(as.character(lm_results_df$`R-squared`))

# Concatenate the names of the comparative methods
lm_results_df$`Comparative Method` <- paste(lm_results_df$`Comparative Method`, lm_results_df$Method2, sep = " - ")
lm_results_df <- lm_results_df %>% select(-Method2)
# Print the results
print(lm_results_df)

# Save the results table (Table S4)
write.csv(lm_results_df, "Table_S4.csv", row.names = FALSE)


### FIGURE 5 ###

# Define a general theme for all the plots (Figure 5)
mytheme <- theme(text=element_text(family='Calibri'),
                 axis.text.y=element_text(colour="black", size=11, face="plain"), 
                 axis.text.x=element_text(colour="black", face="plain", size=11),
                 legend.text=element_text(size=11, face="plain", colour="black"), 
                 axis.title.x=element_text(size =13, colour ="black"),
                 axis.title.y=element_text(size =13, colour ="black"),
                 panel.background =element_blank(),
                 panel.border=element_rect(colour= "black", fill=NA, size=0.8),
                 legend.key=element_blank(),
                 axis.line=element_line(size=0.1, color='black'))

# Convex hull ~ Static 
hull.static <- 
  ggplot(ranges, aes(x=ranges$Static/100000, y=ranges$Convex/100000)) +
  geom_point(color="gray60", size=1.5) +
  geom_smooth(method="lm", formula=y~x, se=FALSE, linetype="dashed", color="red", linewidth=1) +
  geom_abline(slope=1, intercept=0, color="grey", linewidth=0.2) +
  labs(x=expression("Static alpha"), y="Convex hull") +
  scale_x_continuous(limits=c(0, 12), breaks=seq(0, 12, by=3), expand=c(0,0)) +
  scale_y_continuous(limits=c(0, 12), breaks=seq(0, 12, by=3), expand=c(0,0)) +
  mytheme+
  annotate("text", x=3.6, y=11, label=expression(paste(italic(R)^2,"=","0.93")), size=3, hjust=1)

# Convex hull ~ Dynamic 
hull.dynamic <- 
  ggplot(ranges,aes(x=ranges$Dynamic/100000, y=ranges$Convex/100000))+
  geom_point(color="gray60", size=1.5) +
  geom_smooth(method="lm", formula=y~x, se=FALSE, linetype="dashed", color="red", linewidth=1)+
  geom_abline(slope=1, intercept=0, color="grey", linewidth=0.2) +
  labs(x=expression("Dynamic alpha"), y="Convex hull") +
  scale_x_continuous(limits=c(0, 12), breaks=seq(0, 12, by=3), expand=c(0,0)) +
  scale_y_continuous(limits=c(0, 12), breaks=seq(0, 12, by=3), expand=c(0,0)) +
  mytheme +
  annotate("text", x=3.6, y=11, label=expression(paste(italic(R)^2,"=","0.95")), size=3, hjust=1)

# Convex hull ~ Expert
hull.expert <-  
  ggplot(ranges,aes(x=ranges$Expert/100000, y=ranges$Convex/100000))+
  geom_point(color="gray60", size=1.5) +
  geom_smooth(method="lm", formula=y~x, se=FALSE, linetype="dashed", color="red", linewidth=1)+
  geom_abline(slope=1, intercept=0, color="grey", linewidth=0.2) +
  labs(x=expression("Expert map"), y="Convex hull") +
  scale_x_continuous(limits=c(0, 12), breaks=seq(0, 12, by=3), expand=c(0,0)) +
  scale_y_continuous(limits=c(0, 12), breaks=seq(0, 12, by=3), expand=c(0,0)) +
  mytheme +
  annotate("text", x=3.6, y=11, label=expression(paste(italic(R)^2,"=","0.88")), size=3, hjust=1)

# Static ~ Dynamic
static.dynamic <-
  ggplot(ranges,aes(x=ranges$Dynamic/100000, y=ranges$Static/100000))+
  geom_point(color="gray60", size=1.5) +
  geom_smooth(method="lm", formula=y~x, se=FALSE, linetype="dashed", color="red", linewidth=1)+
  geom_abline(slope=1, intercept=0, color="grey", linewidth=0.2) +
  labs(x=expression("Dynamic alpha"), y="Static alpha") +
  scale_x_continuous(limits=c(0, 12), breaks=seq(0, 12, by=3), expand=c(0,0)) +
  scale_y_continuous(limits=c(0, 12), breaks=seq(0, 12, by=3), expand=c(0,0)) +
  mytheme +
  annotate("text", x=3.6, y=11, label=expression(paste(italic(R)^2,"=","0.91")), size=3, hjust=1)

# Static ~ Expert
static.expert <-
  ggplot(ranges,aes(x=ranges$Expert/100000, y=ranges$Static/100000))+
  geom_point(color="gray60", size=1.5) +
  geom_smooth(method="lm", formula=y~x, se=FALSE, linetype="dashed", color="red", linewidth=1)+
  geom_abline(slope=1, intercept=0, color="grey", linewidth=0.2) +
  labs(x=expression("Expert map"), y="Static alpha") +
  scale_x_continuous(limits=c(0, 12), breaks=seq(0, 12, by=3), expand=c(0,0)) +
  scale_y_continuous(limits=c(0, 12), breaks=seq(0, 12, by=3), expand=c(0,0)) +
  mytheme +
  annotate("text", x=3.6, y=11, label=expression(paste(italic(R)^2,"=","0.85")), size=3, hjust=1)

# Dynamic ~ Expert
dynamic.expert <- 
  ggplot(ranges,aes(x=ranges$Expert/100000, y=ranges$Dynamic/100000))+
  geom_point(color="gray60", size=1.5) +
  geom_smooth(method="lm", formula=y~x, se=FALSE, linetype="dashed", color="red", linewidth=1)+
  geom_abline(slope=1, intercept=0, color="grey", linewidth=0.2) +
  labs(x=expression("Expert map"), y ="Dynamic alpha") +
  scale_x_continuous(limits=c(0, 12), breaks=seq(0, 12, by=3), expand=c(0,0)) +
  scale_y_continuous(limits=c(0, 12), breaks=seq(0, 12, by=3), expand=c(0,0)) +
  mytheme +
  annotate("text", x=3.6, y=11, label=expression(paste(italic(R)^2,"=","0.81")), size=3, hjust=1)

# Convex hull ~ SDM 
hull.sdm <-
  ggplot(ranges,aes(x=ranges$SDM/100000, y=ranges$Convex/100000))+
  geom_point(color="gray60", size=1.5) +
  geom_smooth(method="lm", formula=y~x, se=FALSE, linetype="dashed", color="red", linewidth=1)+
  geom_abline(slope=1, intercept=0, color="grey", linewidth=0.2) +
  labs(x=expression("SDMs"), y="Convex hull") +
  scale_x_continuous(limits=c(0, 12), breaks=seq(0, 12, by=3), expand=c(0,0)) +
  scale_y_continuous(limits=c(0, 12), breaks=seq(0, 12, by=3), expand=c(0,0)) +
  mytheme +
  annotate("text", x=3.6, y=11, label=expression(paste(italic(R)^2,"=","0.76")), size=3, hjust=1)

# Static ~ SDM
static.sdm <-
  ggplot(ranges,aes(x=ranges$SDM/100000, y=ranges$Static/100000))+
  geom_point(color="gray60", size=1.5) +
  geom_smooth(method="lm",formula=y~x,se=FALSE,linetype="dashed", color="red",linewidth=1)+
  geom_abline(slope=1, intercept=0, color="grey", linewidth=0.2) +
  labs(x=expression("SDMs"), y="Static alpha") +
  scale_x_continuous(limits=c(0, 12), breaks=seq(0, 12, by=3), expand=c(0,0)) +
  scale_y_continuous(limits=c(0, 12), breaks=seq(0, 12, by=3), expand=c(0,0)) +
  mytheme +
  annotate("text", x=3.6, y=11, label=expression(paste(italic(R)^2,"=","0.73")), size=3, hjust=1)

# Dynamic ~ SDM
dynamic.sdm <-
  ggplot(ranges,aes(x=ranges$SDM/100000, y=ranges$Dynamic/100000))+
  geom_point(color = "gray60", size=1.5) +
  geom_smooth(method="lm",formula=y~x,se=FALSE,linetype="dashed", color="red",linewidth=1)+
  geom_abline(slope=1, intercept=0, color="grey", linewidth=0.2) +
  labs(x=expression("SDMs"), y ="Dynamic alpha") +
  scale_x_continuous(limits=c(0, 12), breaks=seq(0, 12, by=3), expand=c(0,0)) +
  scale_y_continuous(limits=c(0, 12), breaks=seq(0, 12, by=3), expand=c(0,0)) +
  mytheme +
  annotate("text", x=3.6, y=11, label=expression(paste(italic(R)^2,"=","0.68")), size=3, hjust=1)

# SDM ~ Expert
sdm.expert <-
  ggplot(ranges,aes(x=ranges$Expert/100000, y=ranges$SDM/100000))+
  geom_point(color="gray60", size=1.5) +
  geom_smooth(method="lm",formula=y~x,se=FALSE,linetype="dashed", color="red",linewidth=1)+
  geom_abline(slope=1, intercept=0, color="grey", linewidth=0.2) +
  labs(x=expression("Expert map"), y="SDMs") +
  scale_x_continuous(limits=c(0, 12), breaks=seq(0, 12, by=3), expand=c(0,0)) +
  scale_y_continuous(limits=c(0, 12), breaks=seq(0, 12, by=3), expand=c(0,0)) +
  mytheme +
  annotate("text", x=3.6, y=11, label=expression(paste(italic(R)^2,"=","0.72")), size=3, hjust=1)

# Combine all the plots into a grid
Fig.5 <- grid.arrange(hull.static, hull.dynamic, hull.expert, hull.sdm, static.dynamic,
                      static.expert, static.sdm, dynamic.expert, dynamic.sdm, sdm.expert,
                      ncol=3, nrow=4, layout_matrix=rbind(c(1, 2, 3), c(4, 5, 6),
                                                         c(7, 8, 9), c(NA, 10, NA)))

# Save figure
ggsave(Fig.5, filename="Fig_5.tiff", width=18, height=18,
       units="cm", scale=1, dpi=400)
