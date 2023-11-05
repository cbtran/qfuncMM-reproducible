filename <- "full-run/full-run-low-low.log"
log_grad_norm <- c()

con <- file(filename, "r")
while (TRUE) {
  line <- readLines(con, n = 1)
  if (length(line) == 0) {
    break
  }
  if (grepl("^NegLL", line)) {
    lgn <- log(as.numeric(strsplit(line, ", ")[[1]][8]))
    log_grad_norm <- c(log_grad_norm, lgn)
  }
}
close(con)

hist(log_grad_norm)
abline(v = c(log(80), log(90), log(100)), col = "red")

mean(log_grad_norm > log(90), na.rm = T)
