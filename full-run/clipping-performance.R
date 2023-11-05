filename <- "full-run/high-low-temp2.log"
obj_no_clip <- c()
grad_no_clip <- c()

con <- file(filename, "r")
while (TRUE) {
  line <- readLines(con, n = 1)
  if (length(line) == 0) {
    break
  }
  if (grepl("^NegLL", line)) {
    negll <- as.numeric(strsplit(strsplit(line, ": ")[[1]][2], ", ")[[1]][1])
    grad <- as.numeric(strsplit(strsplit(line, ": ")[[1]][2], ", ")[[1]][7])
    obj_no_clip <- c(obj_no_clip, negll)
    grad_no_clip <- c(grad_no_clip, grad)
  }
}
close(con)


filename <- "full-run/high-low-temp.log"
obj_clipped <- c()
grad_clipped <- c()

con <- file(filename, "r")
while (TRUE) {
  line <- readLines(con, n = 1)
  if (length(line) == 0) {
    break
  }
  if (grepl("^NegLL:", line)) {
    negll <- as.numeric(strsplit(line, ": ")[[1]][2])
    obj_clipped <- c(obj_clipped, negll)
  } else if (grepl("^Grad norm:", line)) {
    grad <- as.numeric(strsplit(line, ": ")[[1]][2])
    grad_clipped <- c(grad_clipped, grad)
  }
}
close(con)

clip <- 100

png("full-run/clipping-performance.png", width = 800, height = 400)
par(mfrow = c(1, 2))

plot(log(obj_no_clip),
     type = "l", col = "blue", main = "Log loss",
     xlim = c(0, max(length(obj_no_clip), length(obj_clipped))),
     ylim = c(min(log(obj_no_clip), log(obj_clipped)), max(log(obj_no_clip), log(obj_clipped))))
points(log(obj_clipped), type = "l", col = "red")
legend("topright", legend = c("No clipping", "Clipping"), col = c("blue", "red"), lty = 1)

plot(log(grad_no_clip),
     type = "l", col = "blue", main = paste0("Log gradient norm,", " clip = ", clip),
     xlim = c(0, max(length(grad_no_clip), length(grad_clipped))),
     ylim = c(min(log(grad_no_clip), log(grad_clipped)), max(log(grad_no_clip), log(grad_clipped))))
points(log(grad_clipped), type = "l", col = "red")
abline(h = log(clip), lty = 2)
legend("bottomleft", legend = c("No clipping", "Clipping"), col = c("blue", "red"), lty = 1)
dev.off()