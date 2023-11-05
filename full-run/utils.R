kernel_dict <- function(name) {
  switch(name,
    rbf = {
      return(0L)
    },
    matern_1_2 = {
      return(1L)
    },
    matern_3_2 = {
      return(2L)
    },
    matern_5_2 = {
      return(3L)
    },
    {
      stop(paste("Invalid covariance kernel:", name))
    }
  )
}

sigmoid <- function(x) {
  -log(2 / (x + 1) - 1)
}

sigmoid_inv <- function(x) {
  2 / (1 + exp(-x)) - 1
}

softminus <- function(x) {
  log(exp(x) - 1)
}

pairCA <- function(r1avg, r2avg) {
  r1avgavg <- r1avg - mean(r1avg)
  r2avgavg <- r2avg - mean(r2avg)
  sum(r1avgavg * r2avgavg) / (sd(r1avg) * sd(r2avg) * length(r1avg))
}

computeCA <- function(signal) {
  mean_signal <- lapply(signal, \(regmat) apply(regmat, 1, mean))
  result <- c(pairCA(mean_signal$region1, mean_signal$region2),
    pairCA(mean_signal$region1, mean_signal$region3),
    pairCA(mean_signal$region2, mean_signal$region3)) |> sigmoid()
  names(result) <- c("r12", "r13", "r23")
  result
}
