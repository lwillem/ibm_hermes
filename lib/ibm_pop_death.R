flt <- read.fwf("fltper_1x1.txt",
                widths = c(6, 12, 13, 9, 7, 7, 8, 9, 10, 5),
                header = FALSE,
                skip = 1,
                col.names = c("Year", "Age", "mx", "qx", "ax", "lx", "dx", "Lx", "Tx", "ex"))
mlt <- read.fwf("mltper_1x1.txt",
                widths = c(6, 12, 13, 9, 7, 7, 8, 9, 10, 5),
                header = FALSE,
                skip = 1,
                col.names = c("Year", "Age", "mx", "qx", "ax", "lx", "dx", "Lx", "Tx", "ex"))

gen_pop <- function(date, age, tbl)
{
  # Search for relevant cohort
  yr <- floor(date - age)
  gpyr <- subset(tbl, Year == yr)
  fage <- floor(age) # age rounded below
  gpyrr <- subset(gpyr, Age >= fage) # remaining ages
  # Find t such that exp(-H(age + t)) = U, with U random Un(0, 1)
  thr <- -log(runif(1))
  x <- gpyrr$Age
  x[1] <- age # replace first instance (fage) by actual age
  mx0 <- gpyrr$mx[-nrow(gpyrr)] # hazard rates until age 110
  mx1 <- gpyrr$mx[nrow(gpyrr)] # hazard rate from age 110 onwards
  y <- c(0, cumsum(mx0 * diff(x)))
  idx <- max(which(y <= thr))
  x0 <- x[idx]
  if (x0 < 110) {
    y0 <- y[idx]
    y1 <- y[idx + 1]
    x1 <- x[idx + 1]
    m <- (y1 - y0) / (x1 - x0)
    res <- x0 + (thr - y0) / m
  } else { # going beyond 110, assuming constant rate mx1
    y0 <- y[idx]
    m <- mx1
    res <- x0 + (thr - y0) / m
  }
  return(res)
}
date <- 2026.5
age <- 40.3
tbl <- flt
gen_pop(date, age, tbl)

M <- 100000
gp <- rep(NA, M)
library(tictoc)
tic("Generating 100,000 population death times")
for (j in 1:M) {
  gp[j] <- gen_pop(date, age, tbl)
}
toc()
summary(gp)
plot(sort(gp), 1 - (1:M)/M, type = "s",
     xlab = "Age", ylab = "Proportion still alive")
