###
# Nicely formatted plots and tables for the manuscript
# Zane Billings
# 2024-05-05
# Read in all of the raw data and make the necessary figures and tables
# to show the results in the manuscript.
###

# Setup ====

## Package loading ====
library(ggplot2)
library(ggpubr)
library(patchwork)
library(phangorn)

## Data loading ====
dist_df <- readr::read_rds(here::here("Results", "dist-df-unnormalized.Rds"))

h1_treedist <- readr::read_rds(here::here("Results", "h1-tree-dist.Rds"))
h3_treedist <- readr::read_rds(here::here("Results", "h3-tree-dist.Rds"))

## Script loading and other setup ====
source(here::here("R", "utils.R"))
ggplot2::theme_set(hgp::theme_ms())

# Distance measure correlation plot ====
# First we need to clean up the tree distances
full_dist_df <-
	dplyr::bind_rows(
		"h1" = h1_treedist |> tidy_dist_mat(),
		"h3" = h3_treedist |> tidy_dist_mat(),
		.id = "subtype"
	) |>
	dplyr::mutate(method = "cophenetic", .after = subtype) |>
	dplyr::bind_rows(dist_df)

# Make the pairwise correlogram
cor_plot_h1 <-
	full_dist_df |>
	dplyr::filter(subtype == "h1") |>
	tidyr::pivot_wider(
		values_from = d,
		names_from = method
	) |>
	dplyr::select(cophenetic:cart) |>
	GGally::ggpairs(
		columnLabels = c("Tree", "Temporal", "Hamming", "p-Epitope", "Cartography")
	)

cor_plot_h3 <-
	full_dist_df |>
	dplyr::filter(subtype == "h3") |>
	tidyr::pivot_wider(
		values_from = d,
		names_from = method
	) |>
	dplyr::select(cophenetic:cart) |>
	GGally::ggpairs(
		columnLabels = c("Tree", "Temporal", "Hamming", "p-Epitope", "Cartography")
	)

# Make a table with the cor.test info that we can put in text.
pairwise_methods <-
	combn(unique(full_dist_df$method), 2) |>
	t() |>
	`colnames<-`(c("m1", "m2")) |>
	tibble::as_tibble()

dist_df_wide <-
	full_dist_df |>
	tidyr::pivot_wider(
		names_from = method,
		values_from = d
	)

cor_h1 <-
	purrr::pmap(
		pairwise_methods,
		\(m1, m2) cor.test(
			dist_df_wide |>
				dplyr::filter(subtype == "h1") |>
				dplyr::pull(m1),
			dist_df_wide |>
				dplyr::filter(subtype == "h1") |>
				dplyr::pull(m2)
		) |>
			broom::tidy()
	) |>
	dplyr::bind_rows() |>
	dplyr::bind_cols(pairwise_methods)

cor_h3 <-
	purrr::pmap(
		pairwise_methods,
		\(m1, m2) cor.test(
			dist_df_wide |>
				dplyr::filter(subtype == "h3") |>
				dplyr::pull(m1),
			dist_df_wide |>
				dplyr::filter(subtype == "h3") |>
				dplyr::pull(m2)
		) |>
			broom::tidy()
	) |>
	dplyr::bind_rows() |>
	dplyr::bind_cols(pairwise_methods)

# Make the cor plot with only the tree distances for the paper
methods <- unique(full_dist_df$method)
nice_methods <- c(
	"Tree distance",
	"Temporal distance",
	"Hamming distance",
	"p-Epitope distance",
	"Cartographic distance"
)
plot_list_h1 <- list()
plot_list_h3 <- list()

cor_h1_plot <-
	cor_h1[1:4, ] |>
	dplyr::mutate(
		m1, m2,
		label = paste0(
			"R = ", sprintf("%.2f", estimate), "; 95% CI (",
			sprintf("%.2f", conf.low), ", ", sprintf("%.2f", conf.high), ")"
		),
		x = c("left"),
		y = "top",
		.keep = "none"
	)

cor_h3_plot <-
	cor_h3[1:4, ] |>
	dplyr::mutate(
		m1, m2,
		label = paste0(
			"R = ", sprintf("%.2f", estimate), "; 95% CI (",
			sprintf("%.2f", conf.low), ", ", sprintf("%.2f", conf.high), ")"
		),
		x = c("left"),
		y = "top",
		.keep = "none"
	)

for (i in (2:length(methods))) {
	m <- methods[[i]]
	plt <-
		dist_df_wide |>
		dplyr::filter(subtype == "h1") |>
		ggplot() +
		aes(x = cophenetic, y = .data[[m]]) +
		geom_point(
			size = 2,
			alpha = 0.5,
			position = position_jitter(width = 0.01, height = 0.01, seed = 132413)
		) +
		ggpp::geom_label_npc(
			data = cor_h1_plot |> dplyr::filter(m2 == m),
			aes(npcx = x, npcy = y, label = label),
			fill = alpha("white", 0.5)
		) +
		labs(
			x = "Tree distance",
			y = nice_methods[[i]]
		)
	plot_list_h1[[i - 1]] <- plt
}

h1_corr_plot <-
	purrr::reduce(plot_list_h1, `+`) +
	patchwork::plot_annotation(
		title = "H1N1",
		theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
	)

for (i in (2:length(methods))) {
	m <- methods[[i]]
	plt <-
		dist_df_wide |>
		dplyr::filter(subtype == "h3") |>
		ggplot() +
		aes(x = cophenetic, y = .data[[m]]) +
		geom_point(
			size = 2,
			alpha = 0.5,
			position = position_jitter(width = 0.01, height = 0.01, seed = 132413)
		) +
		ggpp::geom_label_npc(
			data = cor_h3_plot |> dplyr::filter(m2 == m),
			aes(npcx = x, npcy = y, label = label),
			fill = alpha("white", 0.5)
		) +
		labs(
			x = "Tree distance",
			y = nice_methods[[i]]
		)
	plot_list_h3[[i - 1]] <- plt
}

h3_corr_plot <-
	purrr::reduce(plot_list_h3, `+`) +
	patchwork::plot_annotation(
		title = "H3N2",
		theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
	)

total_corr_plot <- cowplot::plot_grid(
	h1_corr_plot,
	h3_corr_plot,
	nrow = 1
)

ggsave(
	here::here("Results", "Figures", "corr-plot.png"),
	plot = total_corr_plot,
	width = 13,
	height = 7.5
)

# Tree plots ====
# First we'll plot the ML trees
ml_tree_h1 <- readr::read_rds(here::here("Results", "h1-ml-tree.Rds"))
ml_tree_h3 <- readr::read_rds(here::here("Results", "h3-ml-tree.Rds"))

library(ggtree)
h1_ml_tree_plot <-
	midpoint(ml_tree_h1$tree) |>
	ggtree(ladderize = FALSE, layout = "circular") +
	geom_tiplab() +
	labs(title = "H1N1") +
	theme(
		plot.margin = margin(-0.5, -0.5, -0.5, -0.5, "cm"),
		plot.title = element_text(hjust = 0.5, face = "bold", vjust = -15, size = 20)
	)

ggsave(
	here::here("Results", "Figures", "h1-ml-tree.png"),
	plot = h1_ml_tree_plot,
	width = 6.5, height = 6.5
)

h3_ml_tree_plot <-
	midpoint(ml_tree_h3$tree) |>
	ggtree(ladderize = FALSE, layout = "circular") +
	geom_tiplab() +
	labs(title = "H3N2") +
	theme(
		plot.margin = margin(-0.5, -0.5, -0.5, -0.5, "cm"),
		plot.title = element_text(hjust = 0.5, face = "bold", vjust = -15, size = 20)
	)

h3_ml_tree_plot <- h3_ml_tree_plot |> rotate_tree(90)

ggsave(
	here::here("Results", "Figures", "h3-ml-tree.png"),
	plot = h3_ml_tree_plot,
	width = 6.5, height = 6.5
)

both_ml_trees <-
	cowplot::plot_grid(
		h1_ml_tree_plot,
		h3_ml_tree_plot,
		nrow = 1
	) +
	theme(plot.background = element_rect(fill = "white"))

ggsave(
	here::here("Results", "Figures", "ml-trees.png"),
	plot = both_ml_trees,
	width = 13, height = 6.5
)

magick::image_read(here::here("Results", "Figures", "ml-trees.png")) |>
	magick::image_crop(
		geometry = "3861x1930",
		gravity = "center"
	) |>
	magick::image_write(here::here("Results", "Figures", "ml-trees.png"))

# Distance and Tests table ====
h1_sh_test <- readr::read_rds(here::here("Results", "h1-sh-test.Rds"))
h3_sh_test <- readr::read_rds(here::here("Results", "h3-sh-test.Rds"))

h1_treedists <- readr::read_rds(here::here("Results", "h1_treedists.Rds"))
h3_treedists <- readr::read_rds(here::here("Results", "h3_treedists.Rds"))

h1_table <-
	h1_sh_test |>
	tibble::as_tibble() |>
	dplyr::mutate(
		Tree = c(
			"Maximum Likelihood (baseline)", "Temporal Distance", "Hamming distance",
			"p-Epitope distance", "Cartographic distance"
		),
		.before = dplyr::everything()
	) |>
	dplyr::mutate(
		Tree = Tree,
		`log likelihood` = `ln L`,
		`Δll` = `Diff ln L`,
		dplyr::across(
			c(`log likelihood`, `Δll`),
			\(x) sprintf("%.1f", x)
		),
		`SH p-value` = ifelse(
			`p-value` < 0.001,
			"< 0.001",
			sprintf("%.3f", `p-value`)
		),
		`RF distance` = h1_treedists[[2]][, "ml"],
		.keep = "none"
	)

h3_table <-
	h3_sh_test |>
	tibble::as_tibble() |>
	dplyr::mutate(
		Tree = c(
			"Maximum Likelihood (baseline)", "Temporal Distance", "Hamming distance",
			"p-Epitope distance", "Cartographic distance"
		),
		.before = dplyr::everything()
	) |>
	dplyr::mutate(
		Tree = Tree,
		`log likelihood` = `ln L`,
		`Δll` = `Diff ln L`,
		dplyr::across(
			c(`log likelihood`, `Δll`),
			\(x) sprintf("%.1f", x)
		),
		`SH p-value` = ifelse(
			`p-value` < 0.001,
			"< 0.001",
			sprintf("%.3f", `p-value`)
		),
		`RF distance` = h3_treedists[[2]][, "ml"],
		.keep = "none"
	)

h1_table[1, 4] <- "N/A"
h3_table[1, 4] <- "N/A"

joined_stat_table <-
	dplyr::bind_rows(
		"H1N1" = h1_table,
		"H3N2" = h3_table,
		.id = " "
	) |>
	flextable::flextable() |>
	flextable::merge_v(j = 1) |>
	flextable::valign(j = 1, valign = "top") |>
	flextable::fix_border_issues() |>
	flextable::autofit() |>
	flextable::set_caption(paste0(
		"Log likelihood of all constructed trees, along with the decrease in",
		" log likelihood (Δll) from the ML model, the p-value of the ",
		"Shimodaira-Hasegawa test (SH p-value; evaluated on one million bootstrap",
		" resamples), and the Robinson-Foulds distance from the ML",
		" tree (RF distance)."
	))

readr::write_rds(
	joined_stat_table,
	here::here("Results", "Tables", "stat-table.Rds")
)

#### END OF FILE ####
