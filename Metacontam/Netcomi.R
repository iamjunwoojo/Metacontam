# Load necessary library
library(NetCoMi)

# Command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]
sparseMethod_threshold <- as.numeric(args[3])

# Load and preprocess the data
df <- t(read.csv(input_file, sep="\t", row.names = 1))

# Construct the network
net_pears <- netConstruct(df,
                          measure = "pearson",
                          normMethod = "clr",
                          zeroMethod = "multRepl",
                          zeroPar = list(z.delete = FALSE),
                          sparsMethod = "threshold",
                          dataType = "counts",
                          thresh = sparseMethod_threshold,
                          verbose = 3)

# Save the edge information to the output file
write.table(net_pears$edgelist1, file = output_file, sep="\t")
