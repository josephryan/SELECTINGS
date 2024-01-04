library(ape)

# usage: Rscript unroot.R rooted.tre > unrooted.tre

# process cmd line
args <- commandArgs(trailingOnly = TRUE)
rooted <- args[1]
if (length(args) == 0) {
  stop("usage: Rscript script.R ROOTED_NEWICK_TREEFILE")
}

# unroot the tree
tr <- read.tree(rooted)
unrooted <- unroot(tr)

# print unrooted tree
write.tree(unrooted)

# unroot.R  - script to generate an unrooted version of a rooted tree
# AUTHOR    - Joseph F. Ryan <joseph.ryan@whitney.ufl.edu>
# SYNOPSIS  - Rscript unroot.R NEWICK_TREEFILE > unrooted.tre
# BUGS      - Please report them to the author.
# COPYRIGHT - Copyright (C) 2024 Joseph F. Ryan
# LICENSE   - GNU version 3
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
