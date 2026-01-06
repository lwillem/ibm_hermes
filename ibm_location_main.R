############################################################################ #
# This file is part of the SIMID course material
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright (C) 2026 lwillem, SIMID, UNIVERSITY OF ANTWERP, BELGIUM
############################################################################ #

# clear workspace
rm(list = ls())

# load script with functions
source('./lib/lib_simid_course.R')

# # load default parameters
# get_default_parameters()
# 
# # test the creation of a population
# pop <- create_population_matrix(pop_size =  1000,
#                                 target_school_size = target_school_size,
#                                 target_school_ages = target_school_ages,
#                                 target_workplace_size =  target_workplace_size,
#                                 target_workplace_ages =  target_workplace_ages,
#                                 bool_show_demographics = TRUE)
# dim(pop)

# run the IBM
run_ibm_location(pop_size = 2000,
                 bool_show_demographics = FALSE)


