###
# Sequence data cleaning
# Zane Billings
# 2024-03-30
# Read in the raw sequence data that I got from GISAID and UniProt, and make
# sure everything is consistently formatted.
###

# ---- Setup ----

box::use(
	readr,
	readxl,
	here,
	dplyr,
	janitor,
	tidyr,
	msa
)

# ---- Data Loading ----
dat_raw <- readxl::read_xlsx(here::here("data", "full-sequences.xlsx"))
virus_info <- readr::read_rds(here::here("data", "viruses_used.rds"))

# ---- Data Cleaning ----
# Since the short names are in the same order as the long names, we can bind
# the two datasets together, but
# FIXME find a way to make this a real join for robustness
dat_combined <-
	dat_raw |>
	# format the column names
	janitor::clean_names() |>
	# remove the type column since Justin said we don't need to care about it
	dplyr::select(-nucleotide_sequence_type, -nucleotide_sequence_source) |>
	# Join to the virus info dataset
	dplyr::bind_cols(virus_info)

# Formatting sequences correctly
dat_clean <-
	dat_combined |>
	dplyr::mutate(
		# Make all the sequences lowercase
		dplyr::across(
			c(protein_sequence, nucleotide_sequence),
			stringr::str_to_lower
		),
		# Remove all internal whitespace
		dplyr::across(
			c(protein_sequence, nucleotide_sequence),
			\(x) stringr::str_remove_all(x, "\\s")
		),
		# Get the sequence lengths
		protein_length = nchar(protein_sequence),
		nucleotide_length = nchar(nucleotide_sequence),
		# Indicator for the sequence that is not full length
		full_length = (short_name != "MI/15")
	)

# ---- Save cleaned dataset ----
readr::write_rds(
	dat_clean,
	file = here::here("data", "clean-data.rds")
)

# ---- END OF FILE ----
