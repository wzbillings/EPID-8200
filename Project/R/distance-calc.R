###
# Calculate pairwise strain distance matrices
# Zane Billings
# 2024-05-03
# Take the multiple sequence alignments and calculate various distance
# matrices to use for phylogenetics.
###

box::use(
	readr,
	here,
	phangorn,
	Racmacs
)

source(here::here("R", "utils.R"))

# Data loading ====
# First we need to load in the aligned sequence data. For now we'll
# only look at the protein sequences.
# TODO calculate more distances
# TODO do distances for nuclear sequences
align_h1 <- readr::read_rds(here::here("results", "h1-pro-alignment.Rds"))
align_h3 <- readr::read_rds(here::here("results", "h3-pro-alignment.Rds"))

# Convert the alignments to phyDat type for phangorn
phydat_h1 <- phangorn::as.phyDat(align_h1)
phydat_h3 <- phangorn::as.phyDat(align_h3)

# Load the sequence dataframes
seqs_h1 <- readr::read_rds(here::here("data", "h1-seqs-aligned.Rds"))
seqs_h3 <- readr::read_rds(here::here("data", "h3-seqs-aligned.Rds"))

# Extract the protein sequences
prot_h1 <- seqs_h1$pro_aligned
names(prot_h1) <- seqs_h1$short_name

prot_h3 <- seqs_h3$pro_aligned
names(prot_h3) <- seqs_h3$short_name

# Next we need to load the cartography data
racmacs_map_h1 <- Racmacs::read.acmap(here::here("data", "h1_post_all_2d.ace"))
racmacs_map_h3 <- Racmacs::read.acmap(here::here("data", "h3_post_all_2d.ace"))

# Calculate Hamming distance matrix ====
dist_hamming_h1 <-
	phydat_h1 |>
	phangorn::dist.hamming() |>
	as.matrix()

dist_hamming_h3 <-
	phydat_h3 |>
	phangorn::dist.hamming() |>
	as.matrix()

# Calculate p-Epitope distance matrix ====

# Calculate cartography distance matrix ====

# Do it for H1
dist_cart_h1 <- racmaps_map_to_distances(racmacs_map_h1)
# Fix the column and row names to use short names instead
colnames(dist_cart_h1) <- replace_strain_names(colnames(dist_cart_h1))
rownames(dist_cart_h1) <- replace_strain_names(rownames(dist_cart_h1))

# Dor it for H3
dist_cart_h3 <- racmaps_map_to_distances(racmacs_map_h3)

# Calculate year distance matrix
