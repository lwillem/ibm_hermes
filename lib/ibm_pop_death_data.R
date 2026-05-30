############################################################################ #
# This file is part of the individual-based model framework called HERMES
#
# Goal: mortality data used to generate individual death times
#
# This script is distributed in the hope that it will be useful, but without
# any warranty; See the LICENCE.txt for more details.
#
# Copyright (C) 2026 hputter, ldewreede, LUMC, THE NETHERLANDS
############################################################################ #

file_dir_fltper <- paste0("data/", "fltper_1x1.txt")
file_dir_mltper <- paste0("data/", "mltper_1x1.txt")

flt <- read.fwf(file_dir_fltper,
                widths = c(6, 12, 13, 9, 7, 7, 8, 9, 10, 5),
                header = FALSE,
                skip = 1,
                col.names = c("Year", "Age", "mx", "qx", "ax", "lx", "dx", "Lx", "Tx", "ex"))
mlt <- read.fwf(file_dir_mltper,
                widths = c(6, 12, 13, 9, 7, 7, 8, 9, 10, 5),
                header = FALSE,
                skip = 1,
                col.names = c("Year", "Age", "mx", "qx", "ax", "lx", "dx", "Lx", "Tx", "ex"))