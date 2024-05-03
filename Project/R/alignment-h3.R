###
# Multiple Sequence Alignment for H3 strains
# Zane Billings
# 2024-05-03
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

alignment_to_character <- function(msa_alignment) {
	return(as.character(msa_alignment@unmasked))
}

# ---- Nucleotide Alignment (H3) ----
nuc_seqs <- with(dat_h3, rlang::set_names(nucleotide_sequence, short_name)) |>
	# Replace T with U for MSA
	gsub(pattern = "t", replacement = "u")

#nuc_seqs_filtered <- h3_nuc_seqs[names(nuc_seqs != "MI/85")]

nuc_msa <-
	nuc_seqs |>
	#Biostrings::RNAStringSet() |>
	msa::msa(
		method = "Muscle",
		type = "rna",
		order = "input",
		verbose = TRUE
	)

nuc_seqs_aligned <- alignment_to_character(nuc_msa)

# Save the alignment to file
readr::write_rds(
	nuc_msa,
	here::here("results", "h3-rna-alignment.Rds")
)

# ---- Protein Alignment (H3) ----
pro_seqs <- with(dat_h3, rlang::set_names(protein_sequence, short_name))

pro_msa <-
	pro_seqs |>
	msa::msa(
		method = "Muscle",
		type = "protein",
		order = "input",
		verbose = TRUE
	)

pro_seqs_aligned <- alignment_to_character(pro_msa)

# Save the alignment to file
readr::write_rds(
	pro_msa,
	here::here("results", "h3-pro-alignment.Rds")
)

# ---- Bind aligned seqs to df (H3) ----
dat_seqs_h3 <-
	tibble::tibble(
		short_name = dat_h3$short_name,
		nuc_aligned = nuc_seqs_aligned,
		pro_aligned = pro_seqs_aligned
	)

readr::write_rds(
	dat_seqs_h3,
	file = here::here("data", "h3-seqs-aligned.Rds")
)

#### END OF FILE ####
