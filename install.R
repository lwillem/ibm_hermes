############################################################################ #
# This file is part of the individual-based modelling framework called HERMES
#
# Goal: download HERMES helper scripts into a local 'lib/' folder
#
# This script is distributed in the hope that it will be useful, but without
# any warranty; See the LICENCE.txt for more details.
#
# Copyright (C) 2026 lwillem, SIMID, UNIVERSITY OF ANTWERP, BELGIUM
############################################################################ #

# -------------------------------------------------------------------------- -
# USAGE NOTES
# -------------------------------------------------------------------------- -
# - Set your working directory to the desired project root before running
#   this script.
# - All helper files will be downloaded into a local 'lib/' directory.
# - Existing files in 'lib/' with the same name will be overwritten.
# -------------------------------------------------------------------------- -

# Required dependency for interacting with the GitHub API
library(jsonlite)

# GitHub API endpoint pointing to the 'lib' directory of the HERMES repository
github_api_url <- "https://api.github.com/repos/lwillem/ibm_hermes/contents/lib"

# Local destination directory for helper scripts
project_dir <- "."

# Query the GitHub API to retrieve metadata for files in the repository
repo_contents <- fromJSON(github_api_url)

# Select R scripts located in the 'lib' directory
lib_r_files <- repo_contents[
  repo_contents$type == "file" &
    grepl("\\.R$", repo_contents$name),
]

# PREPARE LOCAL DIRECTORY
# Create local project directory if it does not yet exist
dir.create(file.path(project_dir, 'lib'), showWarnings = TRUE, recursive = TRUE)

# DOWNLOAD HELPER FILES
for (i in seq_len(nrow(lib_r_files))) {
  download.file(
    url  = lib_r_files$download_url[i],
    destfile = file.path(project_dir, 'lib', lib_r_files$name[i]),
    mode = "wb"
  )
}

# DOWNLOAD MAIN and DESCRIPTION FILES
repo_meta_content <- fromJSON(dirname(github_api_url))

download.file(
  url = repo_meta_content$download_url[repo_meta_content$name == 'main.R'],
  destfile = file.path(dirname(project_dir), 'main.R'),
  mode = "wb",
)

download.file(
  url = repo_meta_content$download_url[repo_meta_content$name == 'DESCRIPTION.txt'],
  destfile = file.path(project_dir, 'DESCRIPTION.txt'),
  mode = "wb",
)

message(paste0("All HERMES files successfully installed in `", project_dir, "` directory."))
