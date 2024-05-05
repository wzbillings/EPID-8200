###
# Tree Building
# Zane Billings
# 2024-05-02
# This script takes the distance matrices and alignment and builds the ML
# tree as well as the distance-based trees.
###

# Setup ====

set.seed(370)

## Load packages ====
library(msa)
library(phangorn)

## Load scripts and functions ====
source(here::here("R", "utils.R"))

## Load data ====
# Multiple sequence alignments
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

# Loading distance matrices
dists_h1 <- readr::read_rds(here::here("results", "h1-dists.Rds"))
dists_h3 <- readr::read_rds(here::here("results", "h3-dists.Rds"))

# The p-epitope H1N1 matrix messes up phangorn for some reason, but if we add
# tiny real numbers to it then it fixes it. So we'll do that for the h3 one
# also to be fair.
dists_h1[[3]] <- perturb_matrix(dists_h1[[3]])
dists_h3[[3]] <- perturb_matrix(dists_h3[[3]])

# ML tree ====
# First we need to build a maximum likelihood tree using the MSAs for each
# subtype.

## Model test ====
# We know we want to use the FLU substitution model but phangorn can let us
# compare a few nested models as well. So we'll compare the models on both
# H1 and H3 and pick the best compromised model.
model_test_h1 <- phangorn::modelTest(
	phydat_h1,
	model = "FLU",
	G = TRUE,
	I = TRUE,
	k = 20
)

readr::write_rds(
	data.frame(model_test_h1),
	here::here("Results", "h1-model-test.Rds")
)

model_test_h3 <- phangorn::modelTest(
	phydat_h3,
	model = "FLU",
	G = TRUE,
	I = TRUE,
	k = 20
)

readr::write_rds(
	data.frame(model_test_h3),
	here::here("Results", "h3-model-test.Rds")
)

## Tree fitting ====
# Based on the two models, it seems that the FLU + G(20) + I model would work
# well for both, so we'll use that one.
FIT_ML_TREES <- FALSE
if (isTRUE(FIT_ML_TREES)) {
	ml_tree_h1 <- phangorn::pml_bb(
		phydat_h1,
		model = "FLU+G(20)+I",
		rearrangement = "stochastic",
		method = "unrooted",
		site.rate = "gamma_quadrature"
	)
	
	readr::write_rds(
		ml_tree_h1,
		here::here("Results", "h1-ml-tree.Rds")
	)
	
	ml_tree_h3 <- phangorn::pml_bb(
		phydat_h3,
		model = "FLU+G(20)+I",
		rearrangement = "stochastic",
		method = "unrooted",
		site.rate = "gamma_quadrature"
	)
	
	readr::write_rds(
		ml_tree_h3,
		here::here("Results", "h3-ml-tree.Rds")
	)
} else {
	ml_tree_h1 <- readr::read_rds(here::here("Results", "h1-ml-tree.Rds"))
	ml_tree_h3 <- readr::read_rds(here::here("Results", "h3-ml-tree.Rds"))
}

# Distance trees ====
# We'll use the neighbor joining algorithm in order to build unrooted trees
# from all four of our distance matrices. Then we can compare each of these
# trees to the ML trees and each other.
nj_trees_h1 <-
	purrr::map(
		dists_h1,
		phangorn::NJ
	)

nj_trees_h3 <-
	purrr::map(
		dists_h3,
		phangorn::NJ
	)

readr::write_rds(
	nj_trees_h1,
	here::here("Results", "nj-trees-h1.Rds")
)

readr::write_rds(
	nj_trees_h3,
	here::here("Results", "nj-trees-h3.Rds")
)

# Tree likelihood ====
nj_tree_ll_h1 <- purrr::map(
	nj_trees_h1,
	\(tree) phangorn::pml(
		tree,
		phydat_h1,
		site.rate = "gamma_quadrature",
		bf = ml_tree_h1$bf,
		Q = ml_tree_h1$Q,
		inv = ml_tree_h1$inv,
		k = ml_tree_h1$k,
		shape = ml_tree_h1$shape,
		rate = ml_tree_h1$rate,
		model = ml_tree_h1$model,
		ASC = FALSE
	) |>
		(\(x) {message("next model"); return(x)})() |>
		optim.pml()
)

nj_tree_ll_h3 <- purrr::map(
	nj_trees_h3,
	\(tree) phangorn::pml(
		tree,
		phydat_h3,
		site.rate = "gamma_quadrature",
		bf = ml_tree_h3$bf,
		Q = ml_tree_h3$Q,
		inv = ml_tree_h3$inv,
		k = ml_tree_h3$k,
		shape = ml_tree_h3$shape,
		rate = ml_tree_h3$rate,
		model = ml_tree_h3$model,
		ASC = FALSE
	) |>
		(\(x) {message("next model"); return(x)})() |>
		optim.pml()
)

# Calculating SH test statistics ====
h1_trees_list <- c("ml" = list(ml_tree_h1), nj_tree_ll_h1)
h1_sh_test <-
	phangorn::SH.test(
		h1_trees_list,
		B = 1e7
	)

h3_trees_list <- c("ml" = list(ml_tree_h3), nj_tree_ll_h3)
h3_sh_test <-
	phangorn::SH.test(
		h3_trees_list,
		B = 1e7
	)

readr::write_rds(
	h1_sh_test,
	here::here("Results", "h1-sh-test.Rds")
)

readr::write_rds(
	h3_sh_test,
	here::here("Results", "h3-sh-test.Rds")
)

# Extract tree-based distances from ML trees====
h1_tree_dist <- ml_tree_h1$tree |> ape::cophenetic.phylo()
h3_tree_dist <- ml_tree_h3$tree |> ape::cophenetic.phylo()

readr::write_rds(
	h1_tree_dist,
	here::here("Results", "h1-tree-dist.Rds")
)

readr::write_rds(
	h3_tree_dist,
	here::here("Results", "h3-tree-dist.Rds")
)

# Delta scores ====
h1_delta <- phydat_h1 |> delta.score(arg = "all")
h3_delta <- phydat_h3 |> delta.score(arg = "all")

readr::write_rds(
	h1_delta,
	here::here("Results", "h1-delta-scores.Rds")
)
readr::write_rds(
	h3_delta,
	here::here("Results", "h3-delta-scores.Rds")
)

# Tree similarity metrics ====
# First we have to make a 'multiPhylo' because it can't coerce it itself
# (rolling eyes emoji)
multi_h1 <-
	h1_trees_list |>
	purrr::map(\(m) m$tree) |>
	do.call(what = ape:::c.phylo)
multi_h3 <-
	h3_trees_list |>
	purrr::map(\(m) m$tree) |>
	do.call(what = ape:::c.phylo)

# Now get the different tree distance metrics for H1 and H3
methods <- c(SPR.dist, RF.dist, wRF.dist, KF.dist, path.dist)
h1_treedists <-
	purrr::map(
		methods,
		\(f) f(multi_h1) |>
			as.matrix()
	)
h3_treedists <-
	purrr::map(
		methods,
		\(f) f(multi_h3) |>
			as.matrix()
	)

readr::write_rds(
	h1_treedists,
	here::here("Results", "h1_treedists.Rds")
)
readr::write_rds(
	h3_treedists,
	here::here("Results", "h3_treedists.Rds")
)

#### END OF FILE ####
