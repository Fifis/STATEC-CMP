## code to prepare `DATASET` dataset goes here

f <- list.files(pattern = "\\.csv$")
source("../R/generic.R")

# v = f[1]
new.names <- NULL
for (v in f) {
  a <- read.csv(v)
  a2 <- makeTS(a)
  n <- gsub("\\.csv$", "", gsub("-", ".", v))
  assign(n, a2)
  save.str <- paste0("usethis::use_data(", n, ", overwrite = TRUE)")
  new.names <- c(new.names, n)
  eval(parse(text = save.str))
  print(v)
}

for (v in new.names) {
  cat("#' @rdname calendars", "\n", sep = "")
  cat('"', v, '"', "\n\n", sep = "")
}
