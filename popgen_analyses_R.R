# Popgen analyses of Furcellaria lumbricalis using R

# ----- Identify clones and divide into MLLs, following procedure by Black, Rippe and Matz (https://doi.org/10.1111/eva.70126) ------
# Their Scripts can be found here: https://github.com/z0on/FLKeys-Coral-Seascape-Gen

# Setup environment
HOME <- "C:/Users/xsevej/OneDrive - University of Gothenburg/Genetics/Steffis pipeline/Analyses 2025" # directory
DATA <- paste0(HOME, "/data/") # data folder
RESULT <- paste0(HOME, "/results/") #result folder
SCRATCH <- paste0(HOME, "/scratch/") #result folder
setwd(HOME)

library(viridis)
library(sf)
library(spatstat)
library(pheatmap)
library(ggplot2)
library(ggdendro)

# Load meta data
meta_df <- read.csv("DATA/furcellaria_meta.csv")

# genetic distances
IBS=as.matrix(read.table("DATA/myresult2.ibsMat"))
dim(IBS)
samples <- meta_df 
samples <- samples$indv
# samples=meta_df$indv

# adding sample names to everything
dimnames(IBS)=list(samples,samples)

# ------ hierarchical clustering tree to look for clones, wrong species collections
hc=hclust(as.dist(IBS),"ave")
plot(hc,cex=0.5) # clustering of samples by IBS (great to detect clones or closely related individuals)
# there are a few "low-hanging" groups of closely related individuals
# Note: to make sure that two samples are clones, their distance must be compared to the distance
# between genotyping replicates (the same sample genotyped twice)
abline(h=0.141,col="red") # this seems like a reasonable "low-hang" threshold for calling related groups
abline(h=0.149,col="red") # this seems like a reasonable conservative threshold for calling related groups
pheatmap(1-IBS)

# Dendrogram with branches going to 0 distance
# Convert to dendrogram and then to dendro_data for ggplot
dend <- as.dendrogram(hc)
dend_data <- dendro_data(dend, type = "rectangle")

# Extract labels and match to population
labels_df <- label(dend_data)
labels_df$population <- meta_df$population[match(labels_df$label, samples)]

# Assign colors to populations
pop_levels <- unique(labels_df$population)
pop_colors <- setNames(viridis(length(pop_levels)), pop_levels)

pop_colors <- c(
  "A-FL-A" = "#ffbb78",  # Light Blue
  "A-FL-M" = "#e6550d",  # Dark Blue
  "K-FL-A" = "#aec7e8",  # Light Orange
  "K-FL-M" = "#084594",  # Dark Orange
  "N-FL-A" = "#98df8a",  # Light Green
  "N-FL-M" = "#006d2c",  # Dark Green
  "R-FL-A" = "#c5b0d5",   # Light Purple
  "R-FL-M" = "#6a3d9a")   # Dark Purple

# Merge colors into segment data
seg_df <- segment(dend_data)

# Add a column to indicate if segment ends at a leaf
seg_df$is_leaf <- seg_df$yend == 0

# For leaf segments, assign population color; otherwise gray
seg_df$color <- ifelse(seg_df$is_leaf,
                       pop_colors[labels_df$population[seg_df$xend]],
                       "gray50")

# Create the ggplot object for horizontal dendrogram
ggplot() +
  geom_segment(data = seg_df,
               aes(x = y, y = x, xend = yend, yend = xend, color = color),
               lineend = "round") +
  scale_color_manual(values = pop_colors, name = "Population") +  # legend title
  scale_y_continuous(breaks = 1:length(labels(dend)),
                     labels = labels_df$label) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 6,
                                   color = pop_colors[labels_df$population])) +
  labs(x = "IBS distance", y = "", title = "") +
  geom_vline(xintercept = 0.141, col = "black") +
  geom_vline(xintercept = 0.149, col="black", linetype = "dashed")
  
# --- retaining a single representative of each highly-related group ---
cuts=cutree(hc,h=0.141)
goods=c();i=1
for (i in unique(cuts)) {
  goods[i]=names(cuts[cuts==i])[1]
}
length(goods)  # how many samples are left?

# For cutoff=0.141
# Since we get some MLLs with several samples from each form we will re-add some
# samples so that those MLLs are represented for both attached and drifting (analogous
# to using the clonecorrect function within site*type strata from the "poppr" package).
# For cutoff=0.141 these MLLs are #22, #24, #31 and #34. These will in addition
# to the attached samples previously retained also each be represented by a drifting sample.
testGoods <- append(goods, "R-FL-M-18", 29)
testGoods <- append(testGoods, "R-FL-M-3", 30)
testGoods <- append(testGoods, "A-FL-M-10", 36)
testGoods <- append(testGoods, "A-FL-M-15", 38)
goods <- testGoods

# For cutoff=0.149, the MLLs are #15, #5, #17, #18, #19. 
# These will have drifting samples added.
testGoods <- append(goods, "R-FL-A-11", 16)
testGoods <- append(testGoods, "R-FL-M-13", 18)
testGoods <- append(testGoods, "R-FL-M-14", 19)
testGoods <- append(testGoods, "R-FL-M-16", 20)
testGoods <- append(testGoods, "A-FL-M-10", 23)
testGoods <- append(testGoods, "A-FL-M-15", 24)
goods <- testGoods

# subsetting all data for only the retained samples
IBS.p=IBS[goods,goods]

# Subset meta file
meta_df.p <- meta_df[meta_df$indv %in% goods, ]

# How many MLLs per population?
counts <- table(meta_df.p$population)
print(counts)
sum(counts)

# Before re-adding drifting samples to ensure they are represented in MLLs:
# A-FL-A   A-FL-M  K-FL-A  K-FL-M   N-FL-A  N-FL-M  R-FL-A   R-FL-M   <- if cutoff = 0.141
#   4        1       15      2        18      3       9       4       total = 56

# After re-adding drifting samples to ensure they are represented in MLLs:
# A-FL-A   A-FL-M  K-FL-A  K-FL-M   N-FL-A  N-FL-M  R-FL-A   R-FL-M   <- if cutoff = 0.141
#   4        3       15      2        18      3       9       6       total = 60

# Before re-adding drifting samples to ensure they are represented in MLLs:
# A-FL-A   A-FL-M  K-FL-A  K-FL-M   N-FL-A  N-FL-M  R-FL-A   R-FL-M   <- if cutoff = 0.149
#   2        0       15      2        15      1       2       1       total = 38

# After re-adding drifting samples to ensure they are represented in MLLs:
# A-FL-A   A-FL-M  K-FL-A  K-FL-M   N-FL-A  N-FL-M  R-FL-A   R-FL-M   <- if cutoff = 0.149
#   2        2       15      2        15      1       3       4       total = 44

# Visualize the new ibs
hc=hclust(as.dist(IBS.p),"ave")
plot(hc,cex=0.5)
pheatmap(1-IBS.p)

# Save new ibsMat and meta_df. Do this process for both cutoffs, chaining cuts
write.table(IBS.p, paste0(DATA, "myresult2_noclones_cutoff0141.ibsMat"), row.names=T)
write.csv(meta_df.p, paste0(DATA, "furcellaria_meta_noclones_cutoff0141.csv"))


# ----- Produce PCoA on new IBS mat -----
library(vegan)

HOME <- "C:/Users/xsevej/OneDrive - University of Gothenburg/Genetics/Steffis pipeline/Analyses 2025" # directory
DATA <- paste0(HOME, "/data/") # data folder
RESULT <- paste0(HOME, "/results/") #result folder
SCRATCH <- paste0(HOME, "/scratch/") #result folder
setwd(HOME)
meta_df <- read.csv("DATA/furcellaria_meta_noclones_cutoff0141.csv")
bams <- meta_df$indv # list of bam files
goods=c(1:length(bams))

# reading table of pairs of replicates (tab-delimited) - skip if there are no clones
# clonepairs=read.table("clonepairs.tab",sep="\t")
# repsa= clonepairs[,1]
# repsb= clonepairs[,2]
# removing "b" replicates
# goods=which(!(bams %in% repsb))

# loading individual to population correspondences
i2p <- meta_df[,c("indv", "population")] # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list or vcf file.
row.names(i2p)=i2p[,1]
i2p=i2p[goods,]
population=i2p[,2]

# settign up colors for plotting
palette(rainbow(length(unique(population))))
colors=as.numeric(as.factor(population))
colpops=as.numeric(as.factor(sort(unique(population))))

# clustering / PCoA based on identity by state (IBS) based on single read resampling
# (for low and/or uneven coverage)

ma = as.matrix(read.table("DATA/myresult2_noclones_cutoff0141.ibsMat"))
hc=hclust(as.dist(ma),"ave")
plot(hc,cex=0.5)  # this shows how similar clones are

ma=ma[goods,goods]
dimnames(ma)=list(bams[goods],bams[goods])
hc=hclust(as.dist(ma),"ave")
plot(hc,cex=0.7) # without clones

# performing PCoA and CAP
conds=data.frame(cbind(population))
pp0=capscale(ma~1)
pp=capscale(ma~population,conds)

# significance of by-site divergence
adonis2(ma~population,conds)

# eigenvectors
plot(pp0$CA$eig, main="cutoff = 0.141")

axes2plot=c(1,2)  
library(adegenet) # for transp()
cmd=pp0
plot(cmd,choices=axes2plot,display="sites",type="n") # choices - axes to display
points(cmd,choices=axes2plot,pch=19,col=transp(colors,alpha=0.7)) # add points
# ordihull(cmd,choices= axes2plot,groups= conds$grp,draw="polygon",col=1+as.numeric(unique(as.factor(conds$grp))),label=T)
ordispider(cmd,choices= axes2plot,groups=conds$population,col="grey80")
ordiellipse(cmd,choices= axes2plot,groups= conds$population,draw="polygon",col=colpops,label=T)

# unscaled, to identify outliers
plot(cmd$CA$u[,axes2plot],pch=19,col=colors)
ordispider(cmd$CA$u[,axes2plot],groups=conds$population,col="grey80")
ordiellipse(cmd$CA$u[,axes2plot],groups= conds$population,draw="polygon",col=colpops,label=T)
identify(cmd$CA$u[,axes2plot],labels=colnames(ma),n=3,cex=0.7)

# Make nicer looking graph!
# Extract site scores (coordinates)
site_scores <- as.data.frame(cmd$CA$u)
site_scores$Sample <- rownames(site_scores)

# Merge with metadata
plot_df <- merge(site_scores, meta_df, by.x = "Sample", by.y = "indv")

# Get values for MDS 1 and 2 for graph
cmd$CA$eig[c(1:2)] # For labels in plot below

library(ggplot2)

# Change names for publication
plot_df <- plot_df %>%
  mutate(population = case_when(
    population == "A-FL-A" ~ "B1-A",
    population == "A-FL-M" ~ "B1-D",
    population == "K-FL-A" ~ "S1-A",
    population == "K-FL-M" ~ "S1-D",
    population == "N-FL-A" ~ "S2-A",
    population == "N-FL-M" ~ "S2-D",
    population == "R-FL-A" ~ "S3-A",
    population == "R-FL-M" ~ "S3-D",
    TRUE ~ population
  ))
plot_df$population=factor(plot_df$population)

# Define custom colors for each population (4 sites × 2 types)
custom_colors <- c(
  "B1-A" = "#ffbb78",  # Light Blue
  "B1-D" = "#e6550d",  # Dark Blue
  "S1-A" = "#aec7e8",  # Light Orange
  "S1-D" = "#084594",  # Dark Orange
  "S2-A" = "#98df8a",  # Light Green
  "S2-D" = "#006d2c",  # Dark Green
  "S3-A" = "#c5b0d5",   # Light Purple
  "S3-D" = "#6a3d9a"   # Dark Purple
)

# Define custom shapes
custom_shapes <- c("A" = 17, "M" = 16)

# Create plot. Added - for MDS2 for cutoff=0.149 to match direction of cutoff=0.141
PC0141 <- ggplot(plot_df, aes(x = MDS1, y = -MDS2, color = population, shape = type)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(type = "t", show.legend = FALSE) +
  scale_color_manual(values = custom_colors) +
  scale_shape_manual(values = custom_shapes) +
  theme_classic() +
  labs(
    title = "",
    x = "MDS1 (5.2%)",
    y = "MDS2 (4.5%)",
    color = "Population",
    shape = "Type")+
  theme(legend.title=element_blank())
  #  plot.title = element_text(hjust = 0.5, size = 16),
  #  axis.title = element_text(size = 14),
  #  legend.title = element_text(size = 12),
  #  legend.text = element_text(size = 10)
  #)

grid.arrange(PC0141, PC0149, ncol=2) # 569px x 275 px, or 675 x 250 if you include the legends from R.

# ----- Test isolation by distance (continue from PCoA script above) -----
# IBD using IBS matrix and pairwise distances between sites
# Use IBS distance matrix created above and the site info from meta_df
gen.mat = ma
site <- meta_df$site

sites <- unique(site)
n <- length(sites)
gen.dist.site <- matrix(NA, n, n, dimnames = list(sites, sites))

for (i in 1:n) {
  for (j in 1:n) {
    ind.i <- which(site == sites[i])
    ind.j <- which(site == sites[j])
    gen.dist.site[i, j] <- mean(gen.mat[ind.i, ind.j])
  }
}

# Convert to dist object
gen.dist.site <- as.dist(gen.dist.site)


# --- Distance between sites using marmap ---

# Packages
library(marmap)
library(dplyr)

# Define sites. Use decimal degrees: longitude (x) and latitude (y)
sites <- tibble::tribble(
  ~site,        ~lon,      ~lat,
  "A",     17.645251034376972, 58.812896695829494,  
  "K",     11.613455694001821, 58.12659256331583,    
  "N",     11.683171694864555, 57.88575259143868,  
  "R",     11.7044508997673, 57.82031128209617)

# Least-cost (over water) distances using marmap
# Choose a bathymetry extent that covers all sites with a margin
# Adjust xlim/ylim to region. Resolution 'res' can be 1 or 2 minutes (smaller is finer but slower).
x_margin <- 3
y_margin <- 3
xlim <- range(sites$lon) + c(-x_margin, x_margin)
ylim <- range(sites$lat) + c(-y_margin, y_margin)

# Fetch NOAA ETOPO1 bathymetry
bathy <- getNOAA.bathy(lon1 = xlim[1], lon2 = xlim[2],
                       lat1 = ylim[1], lat2 = ylim[2],
                       resolution = 0.5)  # increase to 1 for finer grid

# Build the transition matrix: cells with depth >= min.depth are passable; land (>= 0) is not.
# For near-shore work, set min.depth slightly negative to ensure shallow water remains passable.
# Note: In marmap, depths are negative (water), and land is 0 or positive.
tm <- trans.mat(bathy, min.depth = -0.5)  # allow paths in very shallow water

# Compute pairwise least-cost distances (km) over water.
# lc.dist expects a data.frame with columns x (lon), y (lat)
coords <- data.frame(x = sites$lon, y = sites$lat)

# Compute pairwise least-cost distances (km) over water
coords <- as.matrix(coords)
rownames(coords) <- c("A", "K", "N", "R")
geo.dist.marmap <- lc.dist(tm, coords, res = "dist", meters = T, round = 0)
print(geo.dist.marmap) # Distances make sense

# (Optional) Quick plot to inspect bathymetry and sites
# plot.bathy is useful to visually confirm your area and that sites are on water
plot(bathy, image = TRUE, land = TRUE, axes = TRUE)
points(sites$lon, sites$lat, pch = 21, bg = "gold", cex = 1.2)
text(sites$lon, sites$lat, labels = sites$site, pos = 3, cex = 0.8)

# Input the distances in same order as gen.dist.site
# Create distance matrix (empty 4x4)
sites <- unique(site)
geo.mat <- matrix(0,
                  nrow = length(sites),
                  ncol = length(sites),
                  dimnames = list(sites, sites))

# Fill matrix with distance between sites in meters
geo.mat["R", "N"] <- 7807
geo.mat["N", "R"] <- geo.mat["R", "N"]

geo.mat["R", "K"] <- 45032
geo.mat["K", "R"] <- geo.mat["R", "K"]

geo.mat["R", "A"] <- 828224
geo.mat["A", "R"] <- geo.mat["R", "A"]

geo.mat["N", "K"] <- 38367
geo.mat["K", "N"] <- geo.mat["N", "K"]

geo.mat["N", "A"] <- 834610
geo.mat["A", "N"] <- geo.mat["N", "A"]

geo.mat["K", "A"] <- 868948
geo.mat["A", "K"] <- geo.mat["K", "A"]

print(geo.mat)

# Convert to a 'dist' object for Mantel
geo.dist.site <- as.dist(geo.mat)

library(vegan)
mantel_res <- mantel(gen.dist.site, geo.dist.site, method = "pearson", permutations = 9999)
mantel_res # p=0.042 for cutoff 0141, p=0.551 for cutoff 0149

# Prepare data for plotting
ibd_df <- data.frame(
  gen = as.vector(gen.dist.site),
  geo = as.vector(geo.dist.site))

# Fit a linear model (genetic ~ geographic)
ibd_lm <- lm(gen ~ geo, data = ibd_df)
summary(ibd_lm)

# Make the IBD plot
ggplot(ibd_df, aes(x = geo, y = gen)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "Geographic distance (m)",
    y = "Genetic distance",
    title = ""
  ) +
  theme_classic(base_size = 14)

# Subset to investigate Attached (A) and Drifting (M) separately. IBS distance matrix can be subset without re-running ANGSD
ibs_by_type <- lapply(split(meta_df, meta_df$type), function(df){
  idx <- which(meta_df$type %in% df$type)
  gen.mat[idx, idx]
})
ibs.A <- ibs_by_type$A

ibs_by_type_2 <- lapply(split(meta_df, meta_df$type), function(df){
  idx <- which(meta_df$type %in% df$type)
  gen.mat[idx, idx]
})
ibs.M <- ibs_by_type_2$M

# ALso split site info into A and M
meta_df.A <- subset(meta_df, meta_df$type=="A")
meta_df.M <- subset(meta_df, meta_df$type=="M")
sites.A <- meta_df.A$site
sites.M <- meta_df.M$site
n <- length(sites)

# Create 2 distance objects
sites <- unique(site)
gen.dist.A <- matrix(NA, n, n, dimnames = list(sites, sites))
gen.dist.M <- matrix(NA, n, n, dimnames = list(sites, sites))

for (i in 1:n) {
  for (j in 1:n) {
    ind.i <- which(sites.A == sites[i])
    ind.j <- which(sites.A == sites[j])
    gen.dist.A[i, j] <- mean(ibs.A[ind.i, ind.j])
  }
}

for (i in 1:n) {
  for (j in 1:n) {
    ind.i <- which(sites.M == sites[i])
    ind.j <- which(sites.M == sites[j])
    gen.dist.M[i, j] <- mean(ibs.M[ind.i, ind.j])
  }
}

# Convert to dist object
gen.dist.A <- as.dist(gen.dist.A)
gen.dist.M <- as.dist(gen.dist.M)

mantel_res.A <- mantel(gen.dist.A, geo.dist.site, method = "pearson", permutations = 999)
mantel_res.M <- mantel(gen.dist.M, geo.dist.site, method = "pearson", permutations = 999)
mantel_res.A # p=0.042 for cutoff 0141, p=0.042 for cutoff 0149.
mantel_res.M # p=0.208 for cutoff 0141, p=0.69 for cutoff 0149.


# Plot A and M in 1 graph

# Prepare data for plotting
dfA <- data.frame(
  geo = as.vector(geo.dist.site),
  gen = as.vector(gen.dist.A),
  type = "A"
)

dfM <- data.frame(
  geo = as.vector(geo.dist.site),
  gen = as.vector(gen.dist.M),
  type = "D"
)

ibd_df_type <- bind_rows(dfA, dfM)
ibd_df_type$type <- factor(ibd_df_type$type, levels=c('D', 'A')) # Sort order of factor
ibd_df_type$geo <- ibd_df_type$geo/1000 # Convert to KM

# Plot IBD for both types
IBD0149 <- ggplot(ibd_df_type, aes(x = geo, y = gen, color = type, shape = type)) +
  geom_smooth(method = "lm", se = TRUE) +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "A" = "#FF9797",
    "D" = "#CA0000" 
  )) +
  scale_fill_manual(values = c(
    "A" = "#FF9797",
    "D" = "#CA0000"
  )) +
  labs(
    x = "Geographic distance (km)",
    y = "Genetic distance",
    title = ""
  ) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = "top")

grid.arrange(IBD0141, IBD0149, ncol=2)# Save as 569 x 300


# ----- Examine log-likelihoods and silhouette width against K to evaluate appropriate K for admixture -----
# Load required libraries
library(ggplot2)
library(dplyr)

# Set working directory to the location of log-likelihood files
setwd(DATA)

# List all log-likelihood files (assuming they are named like 'K2_run1.log', 'K2_run2.log', etc.)
log_files <- list.files(pattern = "myresult4_noclones_cutoff0141_ngsadmix_K\\d+.log")

# Initialize a data frame to store results
log_data <- data.frame(K = integer(), LogLikelihood = numeric())

# Extract K and log-likelihood from each file
for (file in log_files) {
  
  # Read the file and extract the last lines containing the k & log-likelihood
  lines <- readLines(file)
  # Extract K from filename
  k_line <- grep("nPop=", lines, value = TRUE)
  k_value <- as.numeric(sub(".*nPop=([0-9\\.]+).*", "\\1", k_line))
  # Extract log-likelihood
  like_line <- grep("best like=", lines, value = TRUE)
  log_likelihood <- as.numeric(sub(".*best like=([-0-9\\.]+).*", "\\1", like_line))
  
  # Append to data frame
  log_data <- rbind(log_data, data.frame(K = k_value, LogLikelihood = log_likelihood))
}

# Calculate mean and standard deviation for each K
summary_data <- log_data %>%
  group_by(K) %>%
  summarise(MeanLogLikelihood = mean(LogLikelihood),
            SDLogLikelihood = sd(LogLikelihood))

# Plot mean log-likelihood vs. K with error bars
library(scales)

LL0149 <- ggplot(summary_data, aes(x = K, y = MeanLogLikelihood)) +
  geom_line(color = "blue") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = MeanLogLikelihood - SDLogLikelihood,
                    ymax = MeanLogLikelihood + SDLogLikelihood),
                width = 0.2) +
  labs(title = "",
       x = "Number of Ancestral Populations (K)",
       y = "Mean Log-Likelihood") +
  theme_classic()+
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale()))

# We don't really get a plateu for any value of K... but K=3 is the closest.


# -- Measure sillhouette width per K --
# Load required libraries
library(cluster)
library(dplyr)

# Set working directory to where your .qopt files are stored
setwd(DATA)

# List all Q-matrix files (assuming they are named like 'K2_run1.qopt', etc.)
qopt_files <- list.files(pattern = "myresult4_noclones_cutoff0141_ngsadmix_K\\d+.qopt")

# Initialize a data frame to store results
silhouette_data <- data.frame(K = integer(), SilhouetteWidth = numeric())

# Loop through each file
for (file in qopt_files) {
  # Extract K from filename
  k_value <- as.numeric(match(file, qopt_files)+1)
  
  # Read Q-matrix
  qmatrix <- as.matrix(read.table(file))
  
  # Perform clustering (K-means)
  km <- kmeans(qmatrix, centers = k_value, nstart = 10)
  
  # Calculate silhouette width
  sil <- silhouette(km$cluster, dist(qmatrix))
  avg_sil_width <- mean(sil[, "sil_width"])
  
  # Store result
  silhouette_data <- rbind(silhouette_data, data.frame(K = k_value, SilhouetteWidth = avg_sil_width))
}

# Summarize and plot
summary_data <- silhouette_data %>%
  group_by(K) %>%
  summarise(MeanSilhouette = mean(SilhouetteWidth))

# Plot
SW0149 <- ggplot(summary_data, aes(x = K, y = MeanSilhouette)) +
  geom_line(color = "darkgreen") +
  geom_point(size = 3) +
  labs(title = "",
       x = "Number of Ancestral Populations (K)",
       y = "Mean Silhouette Width") +
  theme_classic()

# We get the widest width of silhouette for K=3 for cutoff=0.141, while width is greatest for K=2 for cutoff=0.149, but log-likelihood better for K=3.
library(gridExtra)
grid.arrange(LL0141, LL0149, SW0141, SW0149, ncol=2)
ggsave("AdmixK.svg", grid.arrange(LL0141, LL0149, SW0141, SW0149, ncol=2), units="px", width=569, height=400, dpi=96)

# ----- Plot admixture from NGSadmixture -----
HOME <- "C:/Users/xsevej/OneDrive - University of Gothenburg/Genetics/Steffis pipeline/Analyses 2025" # directory
DATA <- paste0(HOME, "/data/") # data folder
RESULT <- paste0(HOME, "/results/") #result folder
SCRATCH <- paste0(HOME, "/scratch/") #result folder
setwd(HOME)

# loading individual to population correspondences
meta_df <- read.csv("DATA/furcellaria_meta_noclones_cutoff0141.csv")
# meta_df <- read.csv("DATA/furcellaria_meta_noclones_cutoff0149_correct.csv")
i2p <- meta_df[-c(1,2)] # Tab-delimited table containing individual assignments to populations; must be in the same order as samples in the bam list or vcf file.
row.names(i2p)=i2p$indv
# i2p <- i2p[i2p$site == "R", ] # If you need to subset the meta file.

source('SCRATCH/plot_admixture_v5_function.R')

# assembling the input table
inName="myresult4_noclones_cutoff0141_ngsadmix_K3.qopt" # name of the input file to plot, output of ngsAdmix or ADMIXTURE run
tbl=read.table("DATA/myresult4_noclones_cutoff0141_ngsadmix_K3.qopt") # name of the input file to plot, output of ngsAdmix or ADMIXTURE run
pops=i2p # Tab-delimited table containing individual assignments to populations; must be in the same order as samples in the bam list or vcf file.
npops=as.numeric(sub(".+(\\d+)\\..+","\\1",inName))

names(i2p)=c("ind", "bam", "pop", "site", "type")
tbl=cbind(tbl,i2p)
row.names(tbl)=tbl$ind
head(tbl)

# Change pop names for publication
tbl <- tbl %>%
  mutate(pop = case_when(
    pop == "A-FL-A" ~ "B1-A",
    pop == "A-FL-M" ~ "B1-D",
    pop == "K-FL-A" ~ "S1-A",
    pop == "K-FL-M" ~ "S1-D",
    pop == "N-FL-A" ~ "S2-A",
    pop == "N-FL-M" ~ "S2-D",
    pop == "R-FL-A" ~ "S3-A",
    pop == "R-FL-M" ~ "S3-D",
    TRUE ~ pop
  ))
tbl$pop=factor(tbl$pop)

# putting populaitons in desired order (edit pop names as needed or skip to plot them alphabetically)
tbl$pop=factor(tbl$pop,levels=c("S1-A","S1-D","S2-A","S2-D","S3-A","S3-D","B1-A","B1-D")) # Along coast

ords=plotAdmixture(data=tbl, npops=npops, grouping.method="distance", vshift=-0.08, hshift = -2.5, cex=0.92, srt = 90) # Save as 569 px x 300 px.


# --- Make ADMIXTURE graph for each site with K=2 ---
# loading individual to population correspondences
meta_df <- read.csv("DATA/furcellaria_meta_noclones_cutoff0141.csv")
i2p <- meta_df[-c(1,2)] # Tab-delimited table containing individual assignments to populations; must be in the same order as samples in the bam list or vcf file.
row.names(i2p)=i2p$indv
names(i2p)=c("ind", "bam", "pop", "site", "type")

source('SCRATCH/plot_admixture_v5_function.R')

par(mfrow=c(2,2))

for(i in c("K", "N", "R", "A")){
  
  tbl <- read.table(paste0("DATA/",i,"_cutoff0141_ngsadmix_K2.qopt")) # name of the input file to plot, output of ngsAdmix
  i2p_2 <- i2p[i2p$site == i, ]
  tbl=cbind(tbl,i2p_2)
  row.names(tbl)=tbl$ind
  tbl <- tbl %>%
    mutate(pop = case_when(
    pop == "A-FL-A" ~ "B1-A",
    pop == "A-FL-M" ~ "B1-D",
    pop == "K-FL-A" ~ "S1-A",
    pop == "K-FL-M" ~ "S1-D",
    pop == "N-FL-A" ~ "S2-A",
    pop == "N-FL-M" ~ "S2-D",
    pop == "R-FL-A" ~ "S3-A",
    pop == "R-FL-M" ~ "S3-D",
    TRUE ~ pop
    ))
  tbl$pop=factor(tbl$pop)
  pops=i2p_2 # Tab-delimited table containing individual assignments to populations; must be in the same order as samples in the bam list or vcf file.
  plotAdmixture(data=tbl, npops=2, grouping.method="distance", vshift=-0.05, hshift=-2.5, cex=0.92)

}
par(mfrow=c(1,1)) # Save as 569 x 430 px


# ----- Bayescan ------
# assembling the input table
HOME <- "C:/Users/xsevej/OneDrive - University of Gothenburg/Genetics/Steffis pipeline/Analyses 2025" # directory
DATA <- paste0(HOME, "/data/") # data folder
RESULT <- paste0(HOME, "/results/") #result folder
SCRATCH <- paste0(HOME, "/scratch/") #result folder
setwd(HOME)
meta_df <- read.csv("DATA/furcellaria_meta_noclones.csv")
i2p <- meta_df[-c(1,2)] # Tab-delimited table containing individual assignments to populations; must be in the same order as samples in the bam list or vcf file.
row.names(i2p)=i2p$indv

source('SCRATCH/plot_R.r')

#install.packages("boa")
library(boa)
dat=read.table("DATA/best.bayescan.all_fst.txt",header=T)
head(dat)
table(dat[,"qval"]<0.01)
outs=which(dat[,"qval"]<0.01)
plot_bayescan(dat,FDR=0.1,add_text=F,size=0.5,highlight=outs)

mydata=read.table("DATA/best.bayescan.all.sel",colClasses="numeric")
head(mydata)
plot(density(mydata[,parameter]), xlab=parameter, main=paste(parameter,"posterior distribution"),xlim=c(0,0.25),ylim=c(0,500))
boa.hpd(mydata[,parameter],0.05)
lines(density(mydata[,"Fst2"]))
lines(density(mydata[,"Fst3"]))
lines(density(mydata[,"Fst4"]))
boa.hpd(mydata[,parameter],0.05)


# From RPubs (https://rpubs.com/lbenestan/outlier) -> make a nicer graph
bayescan <- read.table("DATA/best.bayescan.cutoff0141_fst.txt",header=T)
#bayescan <- read.table("DATA/best.bayescan.cytoff_fst.txt", header=T)
attach(bayescan)
colnames(bayescan)=c("PROB","LOG_PO","Q_VALUE","ALPHA","FST") 

# Change the value of the Q_VALUE column: 0 == 0.0001.
bayescan$Q_VALUE <- as.numeric(bayescan$Q_VALUE) 
bayescan[bayescan$Q_VALUE<=0.0001,"Q_VALUE"]=0.0001 

# Round the values.
bayescan$LOG_PO <- (round(bayescan$LOG_PO, 4)) 
bayescan$Q_VALUE <- (round(bayescan$Q_VALUE, 4)) 
bayescan$ALPHA <- (round(bayescan$ALPHA, 4)) 
bayescan$FST <- (round(bayescan$FST, 6))

# Add a column for the type of selection grouping based on a Q-VALUE < 0.01.
bayescan$SELECTION <- ifelse(bayescan$ALPHA>=0&bayescan$Q_VALUE<=0.01,"diversifying",ifelse(bayescan$ALPHA>=0&bayescan$Q_VALUE>0.01,"neutral","balancing")) 
bayescan$SELECTION<- factor(bayescan$SELECTION)
levels(bayescan$SELECTION) 

# Add actual position of sites
#vcf <- fread("DATA/myresult4_noclones_cutoff0149.vcf", skip = "#CHROM", select = c(1,2), col.names = c("CHROM", "POS"))
vcf <- fread("DATA/myresult4_noclones_cutoff0141.vcf", skip = "#CHROM", select = c(1,2), col.names = c("CHROM", "POS"))
bayescan <- cbind(bayescan, vcf)

# Save the results of the SNPs potentially under positive (divergent) and balancing selection (qvalue < 0.05).
positive <- bayescan[bayescan$SELECTION=="diversifying",] 
neutral <- bayescan[bayescan$SELECTION=="neutral",] 
balancing <- bayescan[bayescan$SELECTION=="balancing",]

# Check the number of SNPs belonging to each category.
xtabs(data=bayescan, ~SELECTION)

# Transformation Log of the Q value in order to create the ggplot graph.
range(bayescan$Q_VALUE)
bayescan$LOG10_Q <- -log10(bayescan$Q_VALUE)

# Use ggplot to create a nice graph
ggplot(bayescan,aes(x=LOG10_Q, y=FST)) +
  geom_point(aes(fill=SELECTION), pch=21, size=3)+ 
  scale_fill_manual(name="",values=c("orange","red","white"))+ 
  labs(x="Log(q-value)")+ 
  labs(y="Fst")+   
  theme_classic()+
  theme(text = element_text(size=15))

write.table(positive[,c(7,8)], "diversifying_bayescan_cutoff0141.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Manhattan plot
# Load required libraries
library(ggplot2)
library(dplyr)

# Read BayeScan output
bayescan <- read.table("DATA/best.bayescan.all_fst.txt", header = TRUE)

# Add -log10(Q) and significance
bayescan[bayescan$qval<=0.0001,"qval"]=0.0001
bayescan$log10q <- -log10(bayescan$qval)
bayescan$significant <- bayescan$qval < 0.01

# Add locus index for plotting
bayescan$locus <- 1:nrow(bayescan)

# Manhattan plot
manAll <- ggplot(bayescan, aes(x = locus, y = fst, color = significant)) +
  geom_point(size = 1, alpha = 0.6) +
  scale_color_manual(values = c("gray25", "red")) +
  theme_classic() +
  labs(title = "",
       x = "Locus Index",
       y = "Fst",
       color = "Significant (Q < 0.01)") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "none")


# ----- Selection scan PCangsd -----

# Load necessary libraries
library(RcppCNPy)
library(ggplot2)
library(dplyr)
library(qvalue)
library(data.table)

HOME <- "C:/Users/xsevej/OneDrive - University of Gothenburg/Genetics/Steffis pipeline/Analyses 2025" # directory
DATA <- paste0(HOME, "/data/") # data folder
RESULT <- paste0(HOME, "/results/") #result folder
SCRATCH <- paste0(HOME, "/scratch/") #result folder
setwd(HOME)

# --- For all samples before clone correction ---
# Load selection statistics from PCAngsd.
setwd(DATA)
# selection_stats <- npyLoad('pcangsd_cutoff0141_sel.selection.npy') # Previous PCAngsd version saved as .npy
selection_stats <- fread('pcangsd_all_sel.selection')
snp_metadata <- fread('myresult_5.mafs.gz')
sites <- read.csv('pcangsd_all_sel.sites', header=FALSE)
cov_matrix <- read.table("pcangsd_all_sel.cov")
setwd(HOME)
selection_stats <- as.data.frame(selection_stats)
meta_df <- read.csv("DATA/furcellaria_meta_noclones.csv")

# Inspect your PCA to identify PCs which drive the differentiation of interest. This analysis identifies SNPs driving that differentiation.
# Perform PCA
pca_result <- prcomp(cov_matrix, scale. = FALSE)

# Create a data frame for plotting
pca_df <- data.frame(PC1 = pca_result$x[,1],
                     PC2 = pca_result$x[,2])

# Add population info
pca_df$Population <- meta_df$population  # Make sure the order matches!

# Plot with color by population
ggplot(pca_df, aes(x = PC1, y = PC2, color = Population)) +
  geom_point(alpha = 1, cex = 3) +
  labs(title = "PCA from PCAngsd Covariance Matrix",
       x = "PC1",
       y = "PC2") +
  theme_minimal() +
  scale_color_brewer(palette = "Set2")  # You can change the palette if needed

# Convert chi-square statistics to p-values. Here you set which PC you want to identify outlier SNPs for. V1 for PC1, V2 for PC2, and so on.
p_values <- pchisq(selection_stats$V1, df=1, lower.tail=FALSE)

# Apply FDR correction
qobj <- qvalue(p = p_values)
q_values <- qobj$qvalues

# Identify outliers (e.g., q-value < 0.05)
outliers <- which(q_values < 0.01)
length(outliers)

# Add chromosome and position for identifying SNPs
snp_metadata <- as.data.frame(snp_metadata)
snp_metadata_sub <- snp_metadata[sites == 1,]

# Create a data frame for plotting
SNP = 1:length(p_values)
snp_data <- data.frame(
  PValue = p_values,
  QValue = q_values,
  Outlier = SNP %in% outliers,
  chromo = snp_metadata_sub$chromo,
  position = snp_metadata_sub$position
)

outliers <- subset(snp_data, snp_data$Outlier==TRUE)

# Manhattan-style plot
PCAll <- ggplot(snp_data, aes(x = SNP, y = -log10(PValue), color = Outlier)) +
  geom_point(size = 1, alpha = 0.6) +
  scale_color_manual(values = c("gray25", "red")) +
  labs(title = "", x = "Locus index", y = "-log(p-value)") +
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "none")

write.table(outliers[,c(4,5)], "diversifying_PCAngsd_all.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# --- For MLLs of each clonal cutoff ---
# Load selection statistics from PCAngsd
setwd(DATA)
# selection_stats <- npyLoad('pcangsd_cutoff0141_sel.selection.npy') # Previous PCAngsd version saved as .npy
selection_stats <- fread('pcangsd_cutoff0141_sel.selection')
snp_metadata <- fread('myresult4_noclones_cutoff0141.mafs.gz')
sites <- read.csv('pcangsd_cutoff0141_sel.sites', header=FALSE)
cov_matrix <- read.table("pcangsd_cutoff0141_sel.cov")
setwd(HOME)
selection_stats <- as.data.frame(selection_stats)
meta_df <- read.csv("DATA/furcellaria_meta_noclones_cutoff0141.csv")

# Inspect your PCA to identify PCs which drive the differentiation of interest. This analysis identifies SNPs driving that differentiation.
# Perform PCA
pca_result <- prcomp(cov_matrix, scale. = FALSE)

# Create a data frame for plotting
pca_df <- data.frame(PC1 = pca_result$x[,1],
                     PC2 = pca_result$x[,2])

# Add population info
pca_df$Population <- meta_df$population  # Make sure the order matches!

# Plot with color by population
ggplot(pca_df, aes(x = PC1, y = PC2, color = Population)) +
  geom_point(alpha = 1, cex = 3) +
  labs(title = "PCA from PCAngsd Covariance Matrix",
       x = "PC1",
       y = "PC2") +
  theme_minimal() +
  scale_color_brewer(palette = "Set2")  # You can change the palette if needed

# Convert chi-square statistics to p-values. Here you set which PC you want to identify outlier SNPs for. V1 for PC1, V2 for PC2, and so on.
p_values <- pchisq(selection_stats$V1, df=1, lower.tail=FALSE)

# Apply FDR correction
qobj <- qvalue(p = p_values)
q_values <- qobj$qvalues

# Identify outliers (e.g., q-value < 0.05)
outliers <- which(q_values < 0.01)
length(outliers)

# Add chromosome and position for identifying SNPs
snp_metadata <- as.data.frame(snp_metadata)
snp_metadata_sub <- snp_metadata[sites == 1,]

# Create a data frame for plotting
SNP = 1:length(p_values)
snp_data <- data.frame(
  PValue = p_values,
  QValue = q_values,
  Outlier = SNP %in% outliers,
  chromo = snp_metadata_sub$chromo,
  position = snp_metadata_sub$position
)

outliers <- subset(snp_data, snp_data$Outlier==TRUE)

# Manhattan-style plot
PC0141 <- ggplot(snp_data, aes(x = SNP, y = -log10(PValue), color = Outlier)) +
  geom_point(size = 1, alpha = 0.6) +
  scale_color_manual(values = c("gray25", "red")) +
  labs(title = "", x = "Locus index", y = "-log(p-value)") +
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "none")

write.table(outliers[,c(4,5)], "diversifying_PCAngsd_cutoff0149.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)


# --- Selective sweep of MLLs at each site ---
# Load selection statistics from PCAngsd
setwd(DATA)
# selection_stats <- npyLoad('pcangsd_K_cutoff0141_sel.selection.npy')
selection_stats <- fread('pcangsd_K_cutoff0141_sel.selection')
snp_metadata <- fread('K_cutoff0141.mafs.gz')
sites <- read.csv('pcangsd_K_cutoff0141_sel.sites', header=FALSE)
cov_matrix <- read.table("pcangsd_K_cutoff0141_sel.cov")
setwd(HOME)
selection_stats <- as.data.frame(selection_stats)
meta_df <- read.csv("DATA/furcellaria_meta_noclones_cutoff0141.csv")
meta_df_site <- subset(meta_df, meta_df$site=="K")

# Inspect your PCA to identify PCs which drive the differentiation of interest. This analysis identifies SNPs driving that differentiation.
# Perform PCA
pca_result <- prcomp(cov_matrix, scale. = FALSE)

# Create a data frame for plotting
pca_df <- data.frame(PC1 = pca_result$x[,1],
                     PC2 = pca_result$x[,2])

# Add population info
pca_df$type <- meta_df_site$type  # Make sure the order matches!

# Plot with color by population
ggplot(pca_df, aes(x = PC1, y = PC2, color = type)) +
  geom_point(alpha = 1, cex = 3) +
  labs(title = "PCA from PCAngsd Covariance Matrix R",
       x = "PC1",
       y = "PC2") +
  theme_minimal() +
  scale_color_brewer(palette = "Set2")

# Load required packages
library(RcppCNPy)
library(data.table)
library(qvalue)
library(ggplot2)

HOME <- "C:/Users/xsevej/OneDrive - University of Gothenburg/Genetics/Steffis pipeline/Analyses 2025" # directory
DATA <- paste0(HOME, "/data") # data folder

# List of sites
sites_to_analyze <- c("A", "K", "N", "R")

# Function to process one site
process_site <- function(site_name) {
  # Load data for this site
  selection_stats <- fread(file.path(DATA, paste0("pcangsd_", site_name, "_cutoff0141_sel.selection")))
  snp_metadata <- fread(file.path(DATA, paste0(site_name, "_cutoff0141.mafs.gz")))
  sites <- read.csv(file.path(DATA, paste0("pcangsd_", site_name, "_cutoff0141_sel.sites")), header = FALSE)
  cov_matrix <- read.table(file.path(DATA, paste0("pcangsd_", site_name, "_cutoff0141_sel.cov")))
  # Convert to data frame
  selection_stats <- as.data.frame(selection_stats)
  # Meta info
  meta_df <- read.csv(file.path(DATA, "furcellaria_meta_noclones_cutoff0141.csv"))
  meta_df_site <- subset(meta_df, meta_df$site == site_name)
  # PCA
  pca_result <- prcomp(cov_matrix, scale. = FALSE)
  # Chi-square to p-values
  p_values <- pchisq(selection_stats$V1, df = 1, lower.tail = FALSE)
  # FDR correction
  qobj <- qvalue(p = p_values)
  q_values <- qobj$qvalues
  # Outliers
  outliers <- which(q_values < 0.01)
  # SNP metadata subset
  snp_metadata_sub <- as.data.frame(snp_metadata)[sites == 1,]
  # Store outlier SNP info
  outlier_snps <- snp_metadata_sub[outliers, c("chromo", "position")]
  return(outlier_snps)
}

# Run for all sites
outlier_list <- lapply(sites_to_analyze, process_site)
names(outlier_list) <- sites_to_analyze

# Compare overlaps between sites
# Example: find common SNPs across all sites
common_snps <- Reduce(function(x, y) merge(x, y, by = c("chromo", "position")), outlier_list)

# Optional: pairwise overlaps between sites
pairwise_overlaps <- combn(names(outlier_list), 2, function(pair) {
  merge(outlier_list[[pair[1]]], outlier_list[[pair[2]]], by = c("chromo", "position"))
}, simplify = FALSE)

write.table(outlier_list$K[,c(1,2)], "diversifying_PCAngsd_K_cutoff0141.txt", row.names = FALSE, col.names = FALSE)


# ----- Identifying outlier SNPs identified by both Bayescan and PCAngsd -----
HOME <- "C:/Users/xsevej/OneDrive - University of Gothenburg/Genetics/Steffis pipeline/Analyses 2025" # directory
DATA <- paste0(HOME, "/data/") # data folder
RESULT <- paste0(HOME, "/results/") #result folder
SCRATCH <- paste0(HOME, "/scratch/") #result folder
setwd(HOME)

# Input the lists saved above from the separate analysis
A <- read.table("DATA/diversifying_PCAngsd_cutoff0141.txt")
B <- read.table("DATA/diversifying_bayescan_cutoff0141.txt")

# Check if any SNPs match. V1 is chromosome, V2 is position.
in_both <- merge(A, B, by = c("V1", "V2"))

print(in_both)

# ----- Diversity metrics, He, Ho, Fis, genotypic richness, Paretos Beta -----
# Inbreeding coefficient
HOME <- "C:/Users/xsevej/OneDrive - University of Gothenburg/Genetics/Steffis pipeline/Analyses 2025" # directory
DATA <- paste0(HOME, "/data/") # data folder
RESULT <- paste0(HOME, "/results/") #result folder
SCRATCH <- paste0(HOME, "/scratch/") #result folder
setwd(HOME)

tbl <- read.table("DATA/pcangsd_cutoff0141.inbreed.samples")
pops <- read.csv("DATA/furcellaria_meta_noclones_cutoff0141.csv")
pops <- pops[,-c(1,2)]
names(pops)=c("ind", "bam", "pop", "site", "type")
tbl <- cbind(tbl,pops)

boxplot(V1~pop, data=tbl, ylab="Inbreeding F")

summary_Fis <- tbl %>%
  group_by(pop) %>%
  summarise(
    mean = mean(V1, na.rm = TRUE),
    se = sd(V1, na.rm = TRUE) / sqrt(n()))

# Heterozygozity (code from heterozygosity_beagle.R at https://github.com/z0on/2bRAD_denovo)

### if your input is ANGSD beagle genotype likelihoods (-doGlf 2) -> ...beagle.gz
fin = commandArgs(T)
glf <- read.table("DATA/myresult4_noclones_cutoff0141.beagle.gz", header=TRUE)[,-c(1:3)]
glf <- array(c(t(as.matrix(glf))), c(3, ncol(glf)/3, nrow(glf)))

meta_df <- read.csv("DATA/furcellaria_meta_noclones_cutoff0141.csv")
i2p <- meta_df[-c(1,2)] # Tab-delimited table containing individual assignments to populations; must be in the same order as samples in the bam list or vcf file.
row.names(i2p)=i2p$indv
names(i2p)=c("ind", "bam", "pop", "site", "type")

# ### if your input is TEXT beagle genotype probs (-doGeno 8)
# glf <- read.table("myglf.geno.gz", header=FALSE)[,-c(1:2)]
# glf <- array(c(t(as.matrix(glf))), c(3, ncol(glf)/3, nrow(glf)))

# glf[,2,1:10] # CHECKME: geno likes across sfs bins (rows) and first 10 sites (cols) for sample 2

### in all cases, the SAF for a given site is just the genotype likelihoods ...
### from Nielsen et al 2012 PLoS One, we have SAF[0] = GL[0/0]/choose(2,0), SAF[1] = 2*GL[0/1]/choose(2,1), SAF[2] = GL[1/1]/choose(2,0)
### ... and of course everything cancels except the GL[./.]

### so: we can just use the typical EM algorithm from realSFS: 

EMstep <- function (sfs, GL) rowMeans(prop.table(sfs * GL, 2))

SFS <- matrix(1/3,3,dim(glf)[2])

maxiter <- 200
tol <- 1e-8

for(sample in 1:dim(glf)[2]){
  for (iter in 1:maxiter)
  {
    upd <- EMstep (SFS[,sample], glf[,sample,])
    if (sqrt(sum((upd - SFS[,sample])^2)) < tol)
      break;
    SFS[,sample] <- upd
  }
  message(sample," ",round(SFS[2,sample],4))
  if (iter == maxiter) warning("increase maximum number of iterations")
}
print(c(fin,round(summary(SFS[2,]),3)),quote=F)
# save(SFS,file=paste(fin,"myresult3_zygosity.RData",sep="")) # SFS[2,] is estimated heterozygosity per individual
tSFS <- t(SFS) # Transpose to match SampleAndSite.txt

tSFS <- cbind(tSFS[,2], i2p)
colnames(tSFS)=c("heterozygosity", "ind", "bam", "pop", "site", "type")

boxplot(heterozygosity~pop, data=tSFS, ylab="Heterozygosity H")

tSFS %>%
  ggplot(aes(x=pop, y=heterozygosity, fill=pop))+
  geom_boxplot(show.legend=F)+
  theme_minimal()+
  theme(text = element_text(size=15))+
  labs(x='', y='Heterozygosity')+
  geom_jitter(color="black", size=1.5, alpha=0.9)

summary_H <- tSFS %>%
  group_by(pop) %>%
  summarise(
    mean = mean(heterozygosity, na.rm = TRUE),
    se = sd(heterozygosity, na.rm = TRUE) / sqrt(n()))


# --- Calculate Expected Heterozygosity from Ho and Fis ---
He <- cbind(tSFS, tbl$V1)
He$He <- He$heterozygosity / (1 - He$`tbl$V1`)

summary_He <- He %>%
  group_by(pop) %>%
  summarise(
    mean = mean(He, na.rm = TRUE),
    se = sd(He, na.rm = TRUE) / sqrt(n()))


# --- Genotypic richness (R) ---
# Number of MLGs (G) divided by number of individuals (N). R=(G-1)/(N-1).
# Import meta data for cutoff
setwd(HOME)
meta_df <- read.csv("DATA/furcellaria_meta.csv")

# genetic distances
IBS=as.matrix(read.table("DATA/myresult2.ibsMat"))
dim(IBS)
samples <- meta_df$indv

# adding sample names to everything
dimnames(IBS)=list(samples,samples)

# hierarchical clustering tree to look for clones, wrong species collections
hc=hclust(as.dist(IBS),"ave")

# retaining a single representative of each highly-related group
cuts=cutree(hc,h=0.141)
goods=c();i=1
for (i in unique(cuts)) {
  goods[i]=names(cuts[cuts==i])[1]
}
length(goods)  # how many samples are left?

# For cutoff=0.141
# Since we get some MLLs with several samples from each form we will re-add some
# drifting samples so that those MLLs are represented for both attached and drifting.
# For cutoff=0.141 these MLLs are #22, #24, #31 and #34. These will in addition
# to the attached samples retained also each be represented by a drifting sample.
testGoods <- append(goods, "R-FL-M-18", 29)
testGoods <- append(testGoods, "R-FL-M-3", 30)
testGoods <- append(testGoods, "A-FL-M-10", 36)
testGoods <- append(testGoods, "A-FL-M-15", 38)
goods <- testGoods
length(goods)  # how many samples are left?

# For cutoff=0.149, the MLLs are #15, #5, #17, #18, #19. 
# These will have drifting samples added.
testGoods <- append(goods, "R-FL-A-11", 16)
testGoods <- append(testGoods, "R-FL-M-13", 18)
testGoods <- append(testGoods, "R-FL-M-14", 19)
testGoods <- append(testGoods, "R-FL-M-16", 20)
testGoods <- append(testGoods, "A-FL-M-10", 23)
testGoods <- append(testGoods, "A-FL-M-15", 24)
goods <- testGoods
length(goods)  # how many samples are left?

# Subset meta file
meta_df.p <- meta_df[meta_df$indv %in% goods, ]

# How many MLLs per population?
counts <- table(meta_df.p$population)
print(counts)

# Import meta_df without technical replicates for correct number of samples per pop
meta_df.nc <- read.csv("DATA/furcellaria_meta_noclones.csv")

# Calculate R for each population
pops <- unique(meta_df$population)
R_pop <- data.frame(pops=unique(meta_df$population))

# Initialize an empty data frame
summary_R <- data.frame(population = character(), R = numeric(), stringsAsFactors = FALSE)

# Loop through each population
for (i in unique(meta_df$population)) {
  N <- nrow(subset(meta_df.nc, population == i))  # Total individuals
  G <- nrow(subset(meta_df.p, population == i))   # Unique genotypes
  R <- (G-1) / (N-1)                              # Genotypic richness formula
  
  # Append to summary table
  summary_R <- rbind(summary_R, data.frame(Population = i, GenotypicRichness = R))
}


# --- Stats on MLG count table, clonal richness, Paretos Beta etc. ---
# Load required libraries
library(poppr)
library(poweRlaw)
library(dplyr)

# Get meta data
meta_df <- read.csv("DATA/furcellaria_meta.csv")

# Get correct names for publication
meta_df <- meta_df %>%
  mutate(population = case_when(
    population == "A-FL-A" ~ "B1-A",
    population == "A-FL-M" ~ "B1-D",
    population == "K-FL-A" ~ "S1-A",
    population == "K-FL-M" ~ "S1-D",
    population == "N-FL-A" ~ "S2-A",
    population == "N-FL-M" ~ "S2-D",
    population == "R-FL-A" ~ "S3-A",
    population == "R-FL-M" ~ "S3-D",
    TRUE ~ population
  ))

# Define your populations
#populations <- unique(meta_df$population)

# Add MLL assignment to meta_df using cutree results from above
meta_df$MLG <- cuts[meta_df$indv]

# When doing the clone correction technical replicates are included, but we dont
# want these included when calculating stats -> remove so only 1 row per thalli sample in meta_df
meta_df <- meta_df[-c(30,37,45,56,72,84,93,103,109,115,116,119,127,129,132,136,138,140,150),]

# Create count matrix: rows = populations, columns = MLGs, values = number of individuals per MLG
mlg_counts <- table(meta_df$population, meta_df$MLG)
mlg_counts <- as.matrix(mlg_counts)

# --- Look at distribution of MLLs ---
# Determine the maximum count across all populations for consistent y-axis scaling
max_count <- max(mlg_counts)

# Set up plotting area
par(mfrow=c(4,2))
# Loop through each population
for (i in unique(meta_df$population)) {
  # Extract counts for the current population
  counts <- mlg_counts[i, ]
  # Filter out MLGs with zero count
  counts <- counts[counts > 0]
  # Sort counts in decreasing order
  counts <- sort(counts, decreasing = TRUE)
  # Create barplot with consistent y-axis limit
  barplot(counts, main = paste(i), xlab="MLL", ylim=c(0, max_count))
}
par(mfrow=c(1,1)) # Save as 569 px x 800 px

# Function for number of individuals
nInd <- function(x){
  x <- drop(as.matrix(x))
  if (length(dim(x)) > 1){ # if it's a matrix
    res <- rowSums(x)
  } else {                 # if it's a vector
    res <- sum(x)
  }
  return(res)
}

# Function for number of MLGs
nMLG <- function(x){
  x <- drop(as.matrix(x))
  if (length(dim(x)) > 1){ # if it's a matrix
    res <- rowSums(x>0)
  } else {                 # if it's a vector
    res <- sum(x>0)
  }
  return(res)
}

# Function for calculating clonal richness
nCR <- function(x){
  x <- drop(as.matrix(x))
  if (length(dim(x)) > 1){ # if it's a matrix
    res <- (rowSums(x > 0)-1)/(rowSums(x)-1)
  } else {                 # if it's a vector
    res <- (sum(x > 0)-1)/(sum(x)-1)
  }
  return(res)
}

# Function for calculating genotypic eveness D*
D <- function(x) {
  D_star_values <- apply(x, 1, function(counts) {
    # Remove zero counts
    counts <- counts[counts > 0]
    # Total individuals in population
    N <- sum(counts)
    # Proportions
    p <- counts / N
    # Simpson's complement (D*)
    D_star <- 1 - sum(p^2)
    return(D_star)
  })
  return(D_star_values)
}

# Define Pareto beta function
power_law_beta <- function(x){
  xpow <- displ(x[x > 0])  # Create distribution object
  xpow$setPars(estimate_pars(xpow))  # Estimate parameters
  xdat <- plot(xpow, draw = FALSE)  # Extract data
  xlm <- lm(log(y) ~ log(x), data = xdat)  # Fit log-log model
  return(-coef(xlm)[2])  # Return negative slope
}

# Wrapper for diversity_stats
Beta <- function(x){
  x <- drop(as.matrix(x))
  if (length(dim(x)) > 1){
    res <- apply(x, 1, power_law_beta)
  } else {
    res <- power_law_beta(x)
  }
  return(res)
}

div_stat <- diversity_stats(mlg_counts, nInd = nInd, nMLG = nMLG, nCR = nCR, B = Beta, D = D) #
div_stat

# --- Nucleotide diversity and Tajima's D ---
# Load libraries
library(ggplot2)
library(dplyr)
library(readr)

# Import data. Before this, I combined the data for each population (from each .thetas.idx.pestPG file) into 1 table.
# Taking the sum of tP (per-site pi values across chromosomes) divided by the sum of nSites = genome-wide pi per population.
df <- read_csv("DATA/pestPG_cutoff0141.csv")
boxplot(tP~Pop, data=df)

# Summarize: mean and standard error per population
summary_pi <- df %>%
  group_by(Pop) %>%
  summarise(
    mean_pi = mean(tP, na.rm = TRUE),
    se_pi = sd(tP, na.rm = TRUE) / sqrt(n())
  )

# We don't want mean pi per chromosome for each population, we instead want pi per genome-site pooled across all chromosomes per population.
summary_pi_nSite <- df %>%
  group_by(Pop) %>%
  summarise(
    mean_pi = sum(tP, na.rm = TRUE) / sum(nSites)
  )

# --- Output table ---
output <- cbind(div_stat, summary_H$mean, summary_Fis$mean, summary_pi$mean_pi, summary_pi_nSite$mean_pi, summary_He$mean)
colnames(output) <- c("H", "G", "lambda", "E.5", "nInd", "nMLG", "R", "B", "D*", "Ho", "Fis", "pi_chrom", "pi_site", "He")
output
write.csv(output, paste0(DATA, "MLLstats_cutoff_0141.csv"))


# ----- Compare MLG stats between attached and drifting with one-way anova -----
library(performance)
library(car)
library(ggplot2)
library(ggpubr)

HOME <- "C:/Users/xsevej/OneDrive - University of Gothenburg/Genetics/Steffis pipeline/Analyses 2025" # directory
DATA <- paste0(HOME, "/data/") # data folder
RESULT <- paste0(HOME, "/results/") #result folder
SCRATCH <- paste0(HOME, "/scratch/") #result folder
setwd(HOME)

df <- read.table("DATA/MLLstats_cutoff_0141.csv", header=TRUE, sep=",")
colnames(df) <- c("X", "H", "G", "lambda", "E.5", "nInd", "MLL", "R", "B", "D", "Ho", "Fis", "pi_chrom", "pi", "He")
df$type <- as.factor(c("A","D","A","D","A","D","A","D"))
df$site <- as.factor(c("B1","B1","S1","S1","S2","S2","S3","S3"))

# Data exploration
par(mfrow=c(2,4))
for (i in c("MLL", "R","B", "Ho", "Fis", "pi", "D", "He")){
  boxplot(df[,i]~df$type, main=paste(i))
}
par(mfrow=c(1,1))

# Define test function
run_lm_test <- function(df, variable_name) {
  # Create formula dynamically
  formula <- as.formula(paste(variable_name, "~ type"))
  
  # Fit linear model
  model <- lm(formula, data = df)
  
  # Assumption checks
  levene_p <- leveneTest(formula, data = df)$`Pr(>F)`[1]
  shapiro_p <- shapiro.test(resid(model))$p.value
  x <- resid(model)
  y <- pnorm(summary(x), mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))
  ks_p <- ks.test(x, y)$p.value
  # Clean up
  rm(x, y)
  
  # ANOVA result
  anova_result <- anova(model)
  
  # Summary data frame
  assumption_df <- data.frame(
    Variable = variable_name,
    LeveneTest_p = levene_p,
    ShapiroTest_p = shapiro_p,
    KS_Test_p = ks_p
  )
  return(list(Assumptions = assumption_df, ANOVA = anova_result))
}

# Run the test
run_lm_test(df=df, variable_name="MLL")


# Compare Ho and Fis within each site. get tSFS object from where Ho and Fis are calculated

# Ho
# subset data frame to use
site_df <- subset(tSFS, tSFS$site=="R")
run_lm_test(df=site_df, variable_name="heterozygosity")

# Fis
# subset data frame to use
site_df <- subset(tbl, tbl$site=="R")
run_lm_test(df=site_df, variable_name="V1")


# plot results
custom_colors <- c(
  "B1-A" = "#ffbb78",  # Light Blue
  "B1-D" = "#e6550d",  # Dark Blue
  "S1-A" = "#aec7e8",  # Light Orange
  "S1-D" = "#084594",  # Dark Orange
  "S2-A" = "#98df8a",  # Light Green
  "S2-D" = "#006d2c",  # Dark Green
  "S3-A" = "#c5b0d5",   # Light Purple
  "S3-D" = "#6a3d9a"   # Dark Purple
)

shapes <- c(17, 16, 17, 16, 17, 16, 17, 16)

plot_boxpoints_by_type <- function(df, variable_name, site_col = "site", type_col = "type") {
  ggplot(df, aes(x = .data[[type_col]], y = .data[[variable_name]])) +
    geom_boxplot(outlier.shape = NA, alpha = 0.2) +
    geom_jitter(width = 0.1, size = 3, alpha = 0.9, color = custom_colors, shape = shapes) +
    theme_classic() +
    labs(x = "", y = variable_name, color = "site")
}

# Call plot function
plot_boxpoints_by_type(df, "R")

# Run for loop for all variables. Put them in correct order for publication
plot_list <- list()
for (i in c("MLL", "R", "D", "Ho")) {
  plot_list[[i]] <- plot_boxpoints_by_type(df, i)
}

# To add some more space for R since we want to get in significance notations above boxes
plot_list$R <- ggplot(df, aes(x = .data[["type"]], y = .data[["R"]])) +
  geom_boxplot(outlier.shape = NA, alpha = 0.2) +
  geom_jitter(width = 0.1, size = 3, alpha = 0.9, color = custom_colors, shape = shapes) +
  theme_classic() +
  labs(x = "", y = "R", color = "Site")+
  scale_y_continuous(limits=c(-0.1, 1.3))+
  scale_colour_manual(values="black")

# Customize D*
plot_list$D <- ggplot(df, aes(x = .data[["type"]], y = .data[["D"]])) +
  geom_boxplot(outlier.shape = NA, alpha = 0.2) +
  geom_jitter(width = 0.1, size = 3, alpha = 0.9, color = custom_colors, shape = shapes) +
  theme_classic() +
  labs(x = "", y = "D", color = "Site")+
  scale_y_continuous(limits=c(-0.1, 1))+
  scale_colour_manual(values="black")

# Customize MLL
plot_list$MLL <- ggplot(df, aes(x = .data[["type"]], y = .data[["MLL"]])) +
  geom_boxplot(outlier.shape = NA, alpha = 0.2) +
  geom_jitter(width = 0.1, size = 3, alpha = 0.9, color = custom_colors, shape = shapes) +
  theme_classic() +
  labs(x = "", y = "MLL", color = "Site")+
  scale_y_continuous(limits=c(0, 20))+
  scale_colour_manual(values="black")

# Customize Ho
plot_list$Ho <- ggplot(df, aes(x = .data[["type"]], y = .data[["Ho"]])) +
  geom_boxplot(outlier.shape = NA, alpha = 0.2) +
  geom_jitter(width = 0.1, size = 3, alpha = 0.9, color = custom_colors, shape = shapes) +
  theme_classic() +
  labs(x = "", y = "Ho", color = "Site")+
  scale_y_continuous(limits=c(0, 0.4))+
  scale_colour_manual(values="black")

ggarrange(plotlist = plot_list, nrow = 2, ncol = 2) # Save as 569 x 350 px


# ----- Compare Fis between before and after clone correction -----
# Inbreeding coefficient
library(ggplot2)
library(ggpattern)

HOME <- "C:/Users/xsevej/OneDrive - University of Gothenburg/Genetics/Steffis pipeline/Analyses 2025" # directory
DATA <- paste0(HOME, "/data/") # data folder
RESULT <- paste0(HOME, "/results/") #result folder
SCRATCH <- paste0(HOME, "/scratch/") #result folder
setwd(HOME)

# Cutoff 0141
tbl <- read.table("DATA/pcangsd_cutoff0141.inbreed.samples")
pops <- read.csv("DATA/furcellaria_meta_noclones_cutoff0141.csv")
pops <- pops[,-c(1,2)]
colnames(pops)=c("ind", "bam", "pop", "site", "type")
tbl <- cbind(tbl,pops)

summary_Fis_cutoff0141 <- tbl %>%
  group_by(pop) %>%
  summarise(
    mean = mean(V1, na.rm = TRUE),
    se = sd(V1, na.rm = TRUE) / sqrt(n()))

MLLs0141 <- data.frame(MLL = rep("NA", 8))
rownames(MLLs0141) <- summary_Fis_cutoff0141$pop
for (i in unique(tbl$pop)){
  MLLs0141[i,] <- nrow(subset(pops, pops$pop == i))
}

# Cutoff 0149
tbl <- read.table("DATA/pcangsd_cutoff0149.inbreed.samples")
pops <- read.csv("DATA/furcellaria_meta_noclones_cutoff0149.csv")
pops <- pops[,-c(1,2)]
colnames(pops)=c("ind", "bam", "pop", "site", "type")
tbl <- cbind(tbl,pops)

summary_Fis_cutoff0149 <- tbl %>%
  group_by(pop) %>%
  summarise(
    mean = mean(V1, na.rm = TRUE),
    se = sd(V1, na.rm = TRUE) / sqrt(n()))

MLLs0149 <- data.frame(MLL = rep("NA", 8))
rownames(MLLs0149) <- summary_Fis_cutoff0149$pop
for (i in unique(tbl$pop)){
  MLLs0149[i,] <- nrow(subset(pops, pops$pop == i))
}

# All data before clone correction
tbl <- read.table("DATA/myresult2_noclones_pcangsd.inbreed.samples")
pops <- read.csv("DATA/furcellaria_meta.csv")
pops <- pops[,-c(1,2)]
colnames(pops)=c("bam", "pop", "site", "type")

# Subset meta file to remove technical replicates
bams_noclones <- read.table("DATA/bams_noclones")
pops <- pops[pops$bam %in% bams_noclones$V1, ]

tbl <- cbind(tbl,pops)

summary_Fis_alldata <- tbl %>%
  group_by(pop) %>%
  summarise(
    mean = mean(V1, na.rm = TRUE),
    se = sd(V1, na.rm = TRUE) / sqrt(n()))

MLLsalldata <- data.frame(MLL = rep("NA", 8))
rownames(MLLsalldata) <- summary_Fis_alldata$pop
for (i in unique(tbl$pop)){
  MLLsalldata[i,] <- nrow(subset(pops, pops$pop == i))
}


# Create a data frame for plotting
Fis_data <- data.frame(
  pop = c(summary_Fis_cutoff0141$pop, summary_Fis_cutoff0149$pop, summary_Fis_alldata$pop),
  mean = c(summary_Fis_cutoff0141$mean, summary_Fis_cutoff0149$mean, summary_Fis_alldata$mean),
  se = c(summary_Fis_cutoff0141$se, summary_Fis_cutoff0149$se, summary_Fis_alldata$se),
  cutoff = c(rep("cutoff=0141", 8), rep("cutoff=0149", 8), rep("All samples", 8)),
  MLL = c(MLLs0141$MLL, MLLs0149$MLL, MLLsalldata$MLL))

# Colors for plot
custom_colors <- c(
  "B1-A" = "#ffbb78",  # Light Blue
  "B1-D" = "#e6550d",  # Dark Blue
  "S1-A" = "#aec7e8",  # Light Orange
  "S1-D" = "#084594",  # Dark Orange
  "S2-A" = "#98df8a",  # Light Green
  "S2-D" = "#006d2c",  # Dark Green
  "S3-A" = "#c5b0d5",   # Light Purple
  "S3-D" = "#6a3d9a")   # Dark Purple

# Fix pop names for plotting
Fis_data <- Fis_data %>%
  mutate(pop = case_when(
    pop == "A-FL-A" ~ "B1-A",
    pop == "A-FL-M" ~ "B1-D",
    pop == "K-FL-A" ~ "S1-A",
    pop == "K-FL-M" ~ "S1-D",
    pop == "N-FL-A" ~ "S2-A",
    pop == "N-FL-M" ~ "S2-D",
    pop == "R-FL-A" ~ "S3-A",
    pop == "R-FL-M" ~ "S3-D",
    TRUE ~ pop))

# Order pops for plotting
Fis_data$pop=factor(Fis_data$pop,levels=c("S1-A","S1-D","S2-A","S2-D","S3-A","S3-D","B1-A","B1-D")) # Along coast

# Plot comparison, save as 569px x 250px
ggplot(Fis_data, aes(x=pop, y=mean, fill=cutoff))+
  geom_bar(stat='identity', position='dodge', colour='black')+
  geom_errorbar(aes(ymin=mean-se, ymax=mean, width=0.2), position=position_dodge(0.9))+
  labs(x='', y='Fis')+
  scale_fill_manual(name='Clone correction', values=c('cutoff=0141'='#8DD2C7', 'cutoff=0149'='#FB7F71', 'All samples'='#BEB9D9'))+
  theme_classic()+
  theme(legend.title=element_blank())


# Define patterns for treatments
pattern_map <- c(
  "cutoff=0141" = "stripe",
  "cutoff=0149" = "circle",
  "All samples" = "none")

# Plot with colors by population and patterns by cutoff
# save as 569px x 250 px
ggplot(Fis_data, aes(x = pop, y = mean,
                     fill = pop, pattern = cutoff)) +
  geom_bar_pattern(stat = "identity",
                   position = "dodge",
                   colour = "black",
                   pattern_fill = "black",  # pattern color
                   pattern_angle = 45,
                   pattern_density = 0.2,
                   pattern_spacing = 0.04) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean),
                position = position_dodge(0.9), width = 0.2) +
  geom_text(data = Fis_data, 
            aes(x = pop, y = mean-se, label = MLL), 
            size = 3, 
            color = "black", 
            hjust = 0.43, 
            vjust = 2, 
            fontface = "bold", 
            position=position_dodge(0.9))+
  labs(x = "", y = "Fis") +
  scale_fill_manual(values = custom_colors) +
  scale_pattern_manual(values = pattern_map) +
  theme_classic() +
  theme(legend.title = element_blank())+
  scale_y_continuous(limits=c(-0.75, 0))


# ----- Map with sampling site and metrics -------
# Load required libraries
library(ggplot2)
library(ggspatial)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggrepel)
library(dplyr)

HOME <- "C:/Users/xsevej/OneDrive - University of Gothenburg/Genetics/Steffis pipeline/Analyses 2025" # directory
DATA <- paste0(HOME, "/data/") # data folder
RESULT <- paste0(HOME, "/results/") #result folder
SCRATCH <- paste0(HOME, "/scratch/") #result folder
setwd(HOME)

# Import table with metrics for each pop
df <- read.csv("DATA/MLLstats_cutoff_0141.csv")
# import meta data
meta_df <- read.csv("DATA/map_meta.csv")
df <- cbind(df, meta_df)

# Get Sweden map
sweden <- ne_countries(scale = "large", returnclass = "sf")
sweden <- subset(sweden, admin %in% c("Sweden", "Norway", "Denmark"))

myColors <- c(
  "B1-A" = "#ffbb78",
  "B1-D" = "#e6550d",
  "S1-A" = "#aec7e8",
  "S1-D" = "#084594",
  "S2-A" = "#98df8a",
  "S2-D" = "#006d2c",
  "S3-A" = "#c5b0d5",
  "S3-D" = "#6a3d9a")

# Plot of whole southern half of Sweden
swemap <- ggplot(data = sweden) +
  geom_sf(fill = "gray90") +
  geom_point(data = df, aes(x = lon, y = lat, size = 5), shape = 18, alpha = 1, color = myColors) +
  #geom_text_repel(data = df, aes(x = lon, y = lat, label = type), size = 3, max.overlaps = Inf) +
  #scale_size_continuous(range = c(3, 10)) +
  coord_sf(xlim = c(10.5, 19), ylim = c(55, 64), expand = FALSE) +
  labs(title = "",
       subtitle = "",
       x = "Longitude", y = "Latitude") +
  theme_classic()+
  annotation_scale()+
  theme(legend.position = "none")

# --- Higher resolution zoomed in on sites ---
# Read in the high-res version of the european map, downloaded from: https://data.dtu.dk/articles/dataset/Shapefile_of_European_countries/23686383?file=41565822
sweden_highres <- st_read("DATA/Europe_merged.shp")

# --- West coast sites ---
# Import table with metrics for each pop
df <- read.csv("DATA/MLLstats_cutoff_0141.csv")
# import meta data
meta_df <- read.csv("DATA/map_meta.csv")
df <- cbind(df, meta_df)

# Define bounding box for the west coast region (adjust as needed)
west_coast_bbox <- st_bbox(c(xmin = 11.4, xmax = 11.9, ymin = 57.78, ymax = 58.16), crs = st_crs(4326))

# Offset types slightly to avoid overlap if you're plotting metrics.
# df <- df %>%
#   mutate(
#     lon = ifelse(type == "A", lon - 0.01, lon + 0),
#     lat = ifelse(type == "A", lat - 0.05, lat + 0.05))

# Convert coordinates to sf object
sites_sf <- st_as_sf(df, coords = c("lon", "lat"), crs = 4326)

# Crop Sweden map and sites to this bounding box
sweden_zoomed <- st_crop(sweden_highres, west_coast_bbox)
sites_zoomed <- st_crop(sites_sf, west_coast_bbox)

myColors2 <- c("#aec7e8", "#084594", "#c5b0d5", "#6a3d9a", "#98df8a", "#006d2c")
shape2 <- c(23, 23, 23, 23)

# Plot the high-resolution map with ggplot2
zoomedwest <- ggplot() +
  geom_sf(data = sweden_zoomed, fill = "lightgrey", color = "grey40") + # High-resolution Sweden map (cropped)
  geom_sf(data = sites_zoomed, aes(fill = site, size = 5), color = myColors2, shape = 23, alpha = 1) + # Sites layer (cropped)
  theme_classic()+
  annotation_scale()+
  theme(legend.position = "none")

zoomedwest

# --- Map of east coast site ---
# Import table with metrics for each pop
df <- read.csv("DATA/MLLstats_cutoff_0141.csv")
# import meta data
meta_df <- read.csv("DATA/map_meta.csv")
df <- cbind(df, meta_df)

# Define bounding box for the east coast region (adjust as needed)
east_coast_bbox <- st_bbox(c(xmin = 17.45, xmax = 17.75, ymin = 58.75, ymax = 58.87), crs = st_crs(4326))

# Offset types slightly to avoid overlap.
# df <- df %>%
#  mutate(
#    lon = ifelse(type == "A", lon - 0.01, lon + 0),
#    lat = ifelse(type == "A", lat - 0.05, lat + 0.05))

# Crop Sweden map and sites to this bounding box
sweden_zoomed_e <- st_crop(sweden_highres, east_coast_bbox)
sites_zoomed_e <- st_crop(sites_sf, east_coast_bbox)

myColors3 <- c("#ffbb78")

# "B1-D" = "#e6550d",

# Plot the high-resolution map with ggplot2
zoomedeast <- ggplot() +
  geom_sf(data = sweden_zoomed_e, fill = "lightgrey", color = "grey40") + # High-resolution Sweden map (cropped)
  geom_sf(data = sites_zoomed_e, aes(fill = site, size = 5), color = myColors3, shape = 23, alpha = 1) + # Sites layer (cropped)
  theme_classic()+
  annotation_scale()+
  theme(legend.position = "none")

zoomedeast

# Combined maps for final figure
library(cowplot)

right_panel <- plot_grid(
  zoomedwest, zoomedeast,
  ncol = 1,
  align = "v")

final_plot <- plot_grid(
  swemap, right_panel,
  ncol = 2,
  rel_widths = c(1, 1))

final_plot # Save as 569 x 673 px for correct dimensions.

