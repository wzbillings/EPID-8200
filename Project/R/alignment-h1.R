###
# Multiple Sequence Alignment for H1 strains
# Zane Billings
# 2024-03-30
# Perform a multiple sequence alignment on the sequences, and calculate
# multiple pairwise distance matrices
# H1 and H3 have to be done separately or for an unknown and unfindable
# reason the MSA package stops working.
###

# ---- Setup ----
# Declare package dependencies
box::use(
	readr,
	here,
	msa[...]
)

# Load data
dat <- readr::read_rds(here::here("data", "clean-data.rds"))

# separate H1 and H3
dat_h1 <- dat |>
	dplyr::filter(subtype == "h1")

dat_h3 <- dat |>
	dplyr::filter(subtype == "h3")

# ---- Nucleotide Alignment (H1) ----
h1_nuc_seqs <- with(dat_h1, rlang::set_names(nucleotide_sequence, short_name)) |>
	# Replace T with U for MSA
	gsub(pattern = "t", replacement = "u")

h1_nuc_msa <-
	h1_nuc_seqs |>
	msa::msa(
		method = "Muscle",
		type = "rna",
		order = "input",
		verbose = TRUE
	)

alignment_to_character <- function(msa_alignment) {
	return(as.character(msa_alignment@unmasked))
}

h1_nuc_seqs_aligned <- alignment_to_character(h1_nuc_msa)

# Save the alignment to file
readr::write_rds(
	h1_nuc_msa,
	here::here("results", "h1-rna-alignment.Rds")
)

# ---- Protein Alignment (H1) ----
h1_pro_seqs <- with(dat_h1, rlang::set_names(protein_sequence, short_name))

h1_pro_msa <-
	h1_pro_seqs |>
	msa::msa(
		method = "Muscle",
		type = "protein",
		order = "input",
		verbose = TRUE
	)

h1_pro_seqs_aligned <- alignment_to_character(h1_pro_msa)

# Save the alignment to file
readr::write_rds(
	h1_pro_msa,
	here::here("results", "h1-pro-alignment.Rds")
)

# ---- Bind aligned seqs to df (H1) ----
dat_seqs_h1 <-
	tibble::tibble(
		short_name = dat_h1$short_name,
		nuc_aligned = h1_nuc_seqs_aligned,
		pro_aligned = h1_pro_seqs_aligned
	)

readr::write_rds(
	dat_seqs_h1,
	file = here::here("data", "h1-seqs-aligned.Rds")
)

#### END OF FILE ####
