library(ape)
library(phangorn)

# Alternative paths (examples)
# path = '/Users/anisus/Library/Mobile Documents/com~apple~CloudDocs/В работе/Deep-Sea Pacific/Rissoidae-Bering/DNA-Rissoidae/R/'
# path = '/Users/anisus/Yandex.Disk.localized/Rissoidae/28S_rissoidae/'

# Main working directory
path = '/Users/anisus3/Library/CloudStorage/OneDrive-Personal/Моллюски Казахстана/Arion vulgaris/Филогения/'

# Read DNA alignment
All_genes = read.dna(paste(path, 'COI_all.fa', sep = ''), format = 'fasta')

# Convert to phyDat format
All_genes_Phy = as.phyDat(All_genes)

# Compute genetic distances under JC69 model
distance = dist.dna(All_genes, model = "JC69")

# Build a Neighbor-Joining tree
tree = NJ(distance)

#### TESTING SUBSTITUTION MODELS ####
model = modelTest(All_genes_Phy, tree, multicore = TRUE)

# Display models sorted by AICc
model[order(model$AICc), ]

#############

# Number of bootstrap replicates
boot_n = 100000

# Fit Maximum Likelihood tree
fit = pml(tree, All_genes_Phy)

# Optimize tree under TPM3u model, with invariant sites and gamma rate heterogeneity
fit = optim.pml(fit, model = "TPM3u", optInv = TRUE, optGamma = TRUE)

# Log-likelihood of the optimized tree
logLik(fit)

# Bootstrap analysis
bs <- bootstrap.pml(
  fit,
  bs = boot_n,
  optNni = TRUE,
  multicore = TRUE,
  control = pml.control(trace = 0)
)

# Create output file path using date, time, and boot_n
current_time <- format(Sys.time(), "%y-%m-%d_%H-%M")
file_name <- paste0("ml-trees/", current_time, "_", boot_n, ".pdf")
file_path <- file.path(path, file_name)
file_path <- file.path(sub("/$", "", path), file_name)

# Open the PDF device for saving the tree
pdf(file_path, width = 8, height = 6)

# Plot the ML tree with bootstrap support values
bs_tree = plotBS(
  midpoint(fit$tree),
  bs,
  p = 0.5,
  type = "p",
  cex = 0.3
)

# Add a scale bar
add.scale.bar(cex = 0.3, font = 2, col = "red")

# Close the PDF device
dev.off()

