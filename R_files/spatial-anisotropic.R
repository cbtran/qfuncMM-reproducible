source("R_files/covariances.R")

split_voxels <- function(coords, direction) {
  # Extract the column corresponding to the direction
  dir_column <- coords[, direction]

  # Find the maximum and minimum
  max_coord <- max(dir_column)
  min_coord <- min(dir_column)

  # Calculate the average
  avg_coord <- (max_coord + min_coord) / 2

  # Split the coords into two sets
  set1 <- coords[dir_column <= avg_coord, ]
  set2 <- coords[dir_column > avg_coord, ]

  list(set1 = set1, set2 = set2)
}

get_partition <- function(coords) {
  x_split <- split_voxels(coords, "x")
  y_split <- split_voxels(coords, "y")
  z_split <- split_voxels(coords, "z")
  # Calculate the absolute difference in number of voxels for each split
  x_diff <- abs(nrow(x_split$set1) - nrow(x_split$set2))
  y_diff <- abs(nrow(y_split$set1) - nrow(y_split$set2))
  z_diff <- abs(nrow(z_split$set1) - nrow(z_split$set2))

  # Find the split with the smallest difference
  min_diff <- min(x_diff, y_diff, z_diff)

  if (min_diff == x_diff) {
    return(x_split)
  } else if (min_diff == y_diff) {
    return(y_split)
  } else {
    return(z_split)
  }
}

get_centroid <- function(coords) {
  # Calculate the mean of each dimension
  centroid_x <- mean(coords[, "x"])
  centroid_y <- mean(coords[, "y"])
  centroid_z <- mean(coords[, "z"])

  # Return the centroid as a vector
  c(centroid_x, centroid_y, centroid_z)
}

get_skewed_square_distances <- function(coords, skewdiag) {
  # Create the diagonal matrix
  skew_matrix <- diag(skewdiag)

  # Calculate the pairwise distances
  distances <- matrix(0, nrow = nrow(coords), ncol = nrow(coords))
  for (i in seq_len(nrow(coords))) {
    for (j in i:nrow(coords)) {
      diff_vector <- coords[i, ] - coords[j, ]
      skewed_vector <- skew_matrix %*% diff_vector
      distances[i, j] <- sqrt(sum(skewed_vector^2))
      distances[j, i] <- distances[i, j]  # The distance matrix is symmetric
    }
  }

  distances
}

anisotropic <- function(coords, rate, skewdiag1, skewdiag2) {
  stopifnot(length(skewdiag1) == 3 && length(skewdiag2) == 3)

  partition <- get_partition(coords)
  part1 <- partition$set1
  part2 <- partition$set2
  centroid1 <- get_centroid(part1)
  centroid2 <- get_centroid(part2)
  inv_dist <- function(centroid, voxel) {
    # Calculate the distance between each voxel and the centroid
    dist <- sqrt(sum((voxel - centroid)^2))
    1 / dist
  }
  cendist_p1 <- apply(coords, 1, inv_dist, centroid1)
  cendist_p2 <- apply(coords, 1, inv_dist, centroid2)

  square_dist_p1 <- get_skewed_square_distances(coords, skewdiag1)
  square_dist_p2 <- get_skewed_square_distances(coords, skewdiag2)

  skewmatern_p1 <- get_cor_mat("matern_5_2", square_dist_p1, rate)
  skewmatern_p2 <- get_cor_mat("matern_5_2", square_dist_p2, rate)

  result <- matrix(0, nrow = nrow(coords), ncol = nrow(coords))
  for (i in seq_len(nrow(coords))) {
    for (j in i:nrow(coords)) {
      p1 <- cendist_p1[i] * cendist_p1[j] * skewmatern_p1[i, j]
      p2 <- cendist_p2[i] * cendist_p2[j] * skewmatern_p2[i, j]
      result[i, j] <- p1 + p2
      result[j, i] <- p1 + p2
    }
  }

  result
}
