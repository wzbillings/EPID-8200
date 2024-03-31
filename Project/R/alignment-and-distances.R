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

# ---- Nucleotide Alignment ----
nuc_seqs <- with(dat, rlang::set_names(nucleotide_sequence, short_name)) |>
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

ss_nuc <- test_msa@unmasked
nuc_seqs_aligned <- ss_nuc |> as.character()

# ---- Protein Alignment ----
pro_seqs <- with(dat, rlang::set_names(protein_sequence, short_name))

pro_msa <-
	pro_seqs |>
	msa::msa(
		method = "Muscle",
		type = "protein",
		order = "input",
		verbose = TRUE
	)

ss_pro <- pro_msa@unmasked
pro_seqs_aligned <- ss_pro |> as.character()

# ---- Bind aligned seqs to df ----
dat_seqs <-
	tibble::tibble(
		short_name = dat_clean$short_name,
		nuc_aligned = nuc_seqs_aligned,
		pro_aligned = pro_seqs_aligned
	)

readr::write_rds(
	dat_seqs,
	file = here::here("data", "seqs-aligned.Rds")
)

# ---- Nucleotide distance matrices ----
nuc_dists <- list()

m_vec <- c("hamming", "lv", "osa", "dl")

nuc_dists <-
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


norm_hamming_test <- normalize_matrix(nuc_dists[[1]])

norm_hamming_test  |>
	`rownames<-`(dat_clean$short_name) |>
	`colnames<-`(dat_clean$short_name) |>
	tibble::as_tibble(rownames = "Var1") |>
	tidyr::pivot_longer(
		cols = -Var1,
		names_to = "Var2",
		values_to = "d"
	)

tidy_dist_mat <- function(d, names) {
	out <- d |>
		`rownames<-`(dat_clean$short_name) |>
		`colnames<-`(dat_clean$short_name) |>
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

library(ggplot2)
ggplot2::theme_set(hgp::theme_ms())
norm_hamming_test |>
	tidy_dist_mat() |>
	ggplot() +
	aes(x = Var1, y = Var2, fill = d) +
	geom_tile() +
	scale_fill_viridis_c() +
	theme(
		axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
		legend.key.size = unit(0.06, "npc"),
		legend.title.position = "top"
	)

plt_list <- purrr::map(
	nuc_dists,
	\(x) x |>
		tidy_dist_mat() |>
		ggplot() +
		aes(x = Var1, y = Var2, fill = d) +
		geom_tile() +
		scale_fill_viridis_c() +
		theme(
			axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
			legend.key.size = unit(0.06, "npc"),
			legend.title.position = "top"
		)
)

library(patchwork)
purrr::reduce(plt_list, `+`)
