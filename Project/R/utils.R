###
# Miscellaneous utility functions
# Zane Billings
# 2024-05-03
# Several functions that are needed for data cleaning and processing
###

# Convert an MSA object to a named character vector
alignment_to_character <- function(msa_alignment) {
	return(as.character(msa_alignment@unmasked))
}

# Scale all of the entries in a matrix using min-max normalization
# If the lowest entry is 0 (as a distance matrix with the diagonal) this is
# equivalent to dividing all entries by the max
normalize_matrix <- function(mat) {
	mmax <- max(mat)
	mmin <- min(mat)
	out <- (mat - mmin) / (mmax - mmin)
	return(out)
}

# Convert a distance matrix into ggplot-approved format
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

# Take a racmacs object, extract the coordinates, and return a distance matrix.
racmaps_map_to_distances <- function(map) {
	coords <- Racmacs::agCoords(racmacs_map_h1)
	strains <- rownames(coords)
	x <- coords[, 1]
	y <- coords[, 2]
	
	res <- matrix(
		nrow = length(strains),
		ncol = length(strains),
		dimnames = list(strains, strains)
	)
	
	# Calculate the lower triangle of the distance matrix
	for (i in 2:nrow(res)) {
		for (j in 1:(i-1)) {
			x_part <- (x[[j]] - x[[i]])^2
			y_part <- (y[[j]] - y[[i]])^2
			dist <- sqrt(x_part + y_part)
			res[[i, j]] <- dist
		}
	}
	
	# We know all the diagonal entries are 0 (distance of a point to itself) so
	# set those without doing the calculation.
	diag(res) <- 0
	
	# Exploit R's dist class to automatically fill in the upper triangle
	out <- res |> as.dist() |> as.matrix()
	
	return(out)
}

# Function for replacing strain names ----
# We have to load the virus info data for this to work. It's probably better to
# do this inside the function, but minimally faster and more convenient to
# load it in global scope.
virus_info <- readr::read_csv(
	here::here("data", "UGAFluVac-virus-names.csv"),
	col_types = "fcccil"
)

# This function uses that dataset to swap out strain names and order the
# strain variable correctly
replace_strain_names <- function(x, from = "analysis", to = "short",
																 drop = TRUE) {
	# Find the right column for selecting names from
	if (from == "analysis") {
		from_vec <- virus_info$analysis_name
	} else if (from == "full") {
		from_vec <- virus_info$genbank_strain_name
	} else if (from == "short") {
		from_vec <- virus_info$short_name
	} else {
		stop("'from' should be 'analysis', 'full', or 'short'.")
	}
	
	# Make sure all values of x exist in the virus info table
	if (!(all(x %in% from_vec))) {
		stop(paste0(
			"'x' should be a vector of ", from, " names that exist in the",
			'virus-info sheet.'
		))
	}
	
	# Now get the location in the virus info table for each element of x
	locs <- match(x, from_vec)
	
	# Based on the names argument, get the correct names to return.
	if (to == "analysis") {
		vals <- virus_info$analysis_name[locs]
	} else if (to == "full") {
		vals <- virus_info$genbank_strain_name[locs]
	} else if (to == "short") {
		vals <- virus_info$short_name[locs]
	} else if (to == "subtype") {
		vals <- virus_info$subtype[locs]
	}
	else {
		stop("'to' should be 'analysis', 'full', 'short', or 'subtype'.")
	}
	
	# If requested, remove unseen factor levels
	if (isTRUE(drop)) {
		vals <- forcats::fct_drop(vals)
	}
	
	return(vals)
}

# Function for extracting years from strain names ====
extract_year <- function(analysis_names) {
	year <-
		stringr::str_extract(
			analysis_names,
			"[0-9]{4}$"
		) |>
		as.integer()
	return(year)
}

# Convenience function for calculating matrix of year distances ====
dist.year <- function(names, format = "short") {
	years <-
		names |>
		replace_strain_names(from = format, to = "analysis") |>
		extract_year()
	res <- dist(years, method = "manhattan")
	
	# Remove the special dist class
	out <- as.matrix(res)
	
	# Set nice names
	short_names <- replace_strain_names(names, from = format, to = "short")
	colnames(out) <- short_names
	rownames(out) <- short_names
	
	return(out)
}

#### END OF FILE ####
