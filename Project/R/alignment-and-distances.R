###
# Multiple Sequence Alignment and distance matrix computation
# Zane Billings
# 2024-03-30
# Perform a multiple sequence alignment on the sequences, and calculate
# multiple pairwise distance matrices
###

# ---- Setup ----
# Declare package dependencies
box::use(
	readr,
	here,
	msa
)

# Load data
dat <- readr::read_rds(here::here("data", "clean-data.rds"))

# separate H1 and H3
dat_h1 <- dat |>
	dplyr::filter(subtype == "h1")

dat_h3 <- dat |>
	dplyr::filter(subtype == "h3")

# ---- Nucleotide Alignment (H1) ----
nuc_seqs <- with(dat_h1, rlang::set_names(nucleotide_sequence, short_name)) |>
	# Replace T with U for MSA
	gsub(pattern = "t", replacement = "u")

nuc_msa <-
	nuc_seqs |>
	msa::msa(
		method = "Muscle",
		type = "rna",
		order = "input",
		verbose = TRUE
	)

alignment_to_character <- function(msa_alignment) {
	return(as.character(msa_alignment@unmasked))
}

nuc_seqs_aligned <- alignment_to_character(nuc_msa)

# Save the alignment to file
readr::write_rds(
	nuc_msa,
	here::here("results", "h1-rna-alignment.Rds")
)

# ---- Protein Alignment (H1) ----
pro_seqs <- with(dat_h1, rlang::set_names(protein_sequence, short_name))

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
	nuc_msa,
	here::here("results", "h1-pro-alignment.Rds")
)

# ---- Bind aligned seqs to df (H1) ----
dat_seqs_h1 <-
	tibble::tibble(
		short_name = dat_h1$short_name,
		nuc_aligned = nuc_seqs_aligned,
		pro_aligned = pro_seqs_aligned
	)

readr::write_rds(
	dat_seqs_h1,
	file = here::here("data", "h1-seqs-aligned.Rds")
)

# ---- Nucleotide Alignment (H3) ----
nuc_seqs <- with(dat_h3, rlang::set_names(nucleotide_sequence, short_name)) |>
	# Replace T with U for MSA
	gsub(pattern = "t", replacement = "u")

nuc_msa <-
	nuc_seqs |>
	msa::msa(
		method = "Muscle",
		type = "rna",
		order = "input",
		verbose = TRUE
	)

alignment_to_character <- function(msa_alignment) {
	return(as.character(msa_alignment@unmasked))
}

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
	nuc_msa,
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

# ---- Protein distance matrices ----
pro_dists <- list()

m_vec <- c("hamming", "lv", "osa", "dl")

pro_dists <-
	purrr::map(
		m_vec,
		\(m) stringdist::stringdistmatrix(
			a = dat_seqs$nuc_aligned,
			b = dat_seqs$nuc_aligned,
			method = m
		),
		.progress = "Calculating amino acid string distances."
	)

normalize_matrix <- function(mat) {
	mmax <- max(mat)
	mmin <- min(mat)
	out <- (mat - mmin) / (mmax - mmin)
	return(out)
}

norm_hamming_test <-
	normalize_matrix(nuc_dists[[1]]) |>
	`rownames<-`(dat_h1$short_name) |>
	`colnames<-`(dat_h1$short_name)

tidy_dist_mat <- function(d) {
	out <- d |>
		tibble::as_tibble(rownames = "Var1") |>
		tidyr::pivot_longer(
			cols = -Var1,
			names_to = "Var2",
			values_to = "d"
		) |>
		# Order variable factors
		dplyr::mutate(
			Var1 = forcats::fct_inorder(Var1),
			Var2 = forcats::fct_inorder(Var2) |> forcats::fct_rev()
		)
	
	return(out)
}

# library(ggplot2)
# ggplot2::theme_set(hgp::theme_ms())
# norm_hamming_test |>
# 	tidy_dist_mat() |>
# 	ggplot() +
# 	aes(x = Var1, y = Var2, fill = d) +
# 	geom_tile() +
# 	scale_fill_viridis_c() +
# 	theme(
# 		axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
# 		legend.key.size = unit(0.06, "npc"),
# 		legend.title.position = "top"
# 	)
# 
# plt_list <- purrr::map2(
# 	nuc_dists,
# 	m_vec,
# 	\(x, m) x |>
# 		normalize_matrix() |>
# 		`rownames<-`(dat_h1$short_name) |>
# 		`colnames<-`(dat_h1$short_name) |>
# 		tidy_dist_mat() |>
# 		ggplot() +
# 		aes(x = Var1, y = Var2, fill = d) +
# 		geom_tile() +
# 		scale_fill_viridis_c(
# 			limits = c(0, 1)
# 		) +
# 		theme(
# 			axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
# 			legend.key.size = unit(0.06, "npc"),
# 			legend.title.position = "top"
# 		) +
# 		labs(
# 			x = NULL,
# 			y = NULL,
# 			title = m
# 		)
# )
# 
# library(patchwork)
# purrr::reduce(plt_list, `+`) +
# 	patchwork::plot_layout(guides = "collect")
