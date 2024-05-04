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



norm_hamming_test <-
	normalize_matrix(nuc_dists[[1]]) |>
	`rownames<-`(dat_h1$short_name) |>
	`colnames<-`(dat_h1$short_name)



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