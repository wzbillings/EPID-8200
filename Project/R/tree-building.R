###
# Making trees
# 2024-04-18
# Zane Billings
# Estimates a maximum likelihood tree for the sequence data, along with
# distance-based trees for each of the computed distance matrices.
###

combs <- combn(1:length(nuc_dists), m = 2)

out <- vector(length = ncol(combs), mode = "list")
for (i in 1:ncol(combs)) {
	out[[i]] <- vegan::mantel(
		nuc_dists[[combs[1, i]]],
		nuc_dists[[combs[2, i]]],
		permutations = 9999
	)
}

purrr::map(out, \(x) x$statistic)
