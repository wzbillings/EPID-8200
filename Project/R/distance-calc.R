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
source(here::here("R", "p-epitope-calculator.R"))

dist_pepi_h1 <-
	prot_h1 |>
	dist.pepi(subtype = 'h1')

dist_pepi_h3 <-
	prot_h1 |>
	dist.pepi(subtype = 'h3')

test1 <- tidy_dist_mat(dist_hamming_h1) |>
	dplyr::rename(hamming = d)
test2 <- tidy_dist_mat(dist_pepi_h1) |>
	dplyr::rename(pepi = d)
testj <- dplyr::full_join(test1, test2, by = c('Var1', 'Var2'))

# Calculate cartography distance matrix ====

# Do it for H1
dist_cart_h1 <- racmaps_map_to_distances(racmacs_map_h1)
# Fix the column and row names to use short names instead
colnames(dist_cart_h1) <- replace_strain_names(colnames(dist_cart_h1))
rownames(dist_cart_h1) <- replace_strain_names(rownames(dist_cart_h1))

# Dor it for H3
dist_cart_h3 <- racmaps_map_to_distances(racmacs_map_h3)

# Calculate year distance matrix ====
dist_year_h1 <-
	prot_h1 |>
	names() |>
	dist.year()

dist_year_h3 <-
	prot_h3 |>
	names() |>
	dist.year()

# Output formatting ====
# Make a list of the distance matrices
dists_h1 <- list(
	"year" = dist_year_h1,
	"hamming" = dist_hamming_h1,
	"pepi" = dist_pepi_h1,
	"cart" = dist_cart_h1
)

dists_h3 <- list(
	"year" = dist_year_h3,
	"hamming" = dist_hamming_h3,
	"pepi" = dist_pepi_h3,
	"cart" = dist_cart_h3
)

# Tidy and combine into a nice data frame
dat_dists_h1 <-
	dists_h1 |>
	purrr::map(tidy_dist_mat) |>
	dplyr::bind_rows(.id = "method")

dat_dists_h3 <-
	dists_h3 |>
	purrr::map(tidy_dist_mat) |>
	dplyr::bind_rows(.id = "method")

dat_dists <-
	list(
		'h1' = dat_dists_h1,
		'h3' = dat_dists_h3
	) |>
	dplyr::bind_rows(.id = "subtype")

dat_dists_normalized <-
	dat_dists |>
	dplyr::group_by(subtype, method) |>
	dplyr::mutate(d = (d - min(d)) / (max(d) - min(d))) |>
	dplyr::ungroup()

# Save to file ====

readr::write_rds(
	dists_h1,
	here::here("results", "h1-dists.Rds")
)

readr::write_rds(
	dists_h3,
	here::here("results", "h3-dists.Rds")
)

readr::write_rds(
	dat_dists,
	here::here("results", "dist-df-unnormalized.Rds")
)

readr::write_rds(
	dat_dists_normalized,
	here::here("results", "dist-df-normalized.Rds")
)

#### END OF FILE ####
