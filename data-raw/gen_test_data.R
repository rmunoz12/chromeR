# gen_test_data.R
# Test data for chrome.R
# Copyright (c) 2015 Richard Munoz
#
# This file is part of chromeR.
#
# chromeR is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# chromeR is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# chromeR If not, see http://www.gnu.org/licenses/.

# Chromosome lengths per
# https://en.wikipedia.org/wiki/Human_genome
chrome_bp <- c(249250621, 243199373, 198022430, 191154276, 180915260,
              171115067, 159138663, 146364022, 141213431, 135534747,
              135006516, 133851895, 115169878, 107349540, 102531392,
              90354753,  81195210,  78077248,  59128983,  63025520,
              48129895,  51304566)
c <- chrome_bp

# Input dataframe columns from ersa
# ------------------------
# @param id integer vector
# @param result_id integer vector
# @param indv1 character vector
# @param indv2 character vector
# @param chromosome integer vector
# @param bp_start integer vector
# @param bp_end integer vector
# @param length numeric vector
#
# Note: code assumes no NA values are in the input datafarme

# Share everything
p1 <- data.frame(result_id=rep(1), indv1=rep("A"), indv2=rep("B"),
                 chromosome=seq(22), bp_start=rep(0), bp_end=c,
                 length=runif(22), stringsAsFactors=FALSE)

# Share one segment
p2 <- data.frame(result_id=rep(2), indv1=rep("A"), indv2=rep("C"),
                 chromosome=4, bp_start=c[4]/2, bp_end=c[4]*(3/4),
                 length=runif(1), stringsAsFactors=FALSE)

# Two segments, gap of 2e5
p3 <- data.frame(result_id=rep(3), indv1=rep("A"), indv2=rep("D"),
                 chromosome=10,
                 bp_start=c(c[10]/2, c[10]/2+1e7+2e5),
                 bp_end=c(c[10]/2+1e7, c[10]/2+1e7+2e5+2e7),
                 length=runif(2), stringsAsFactors=FALSE)

# Share segments at beginning, middle, and end of chromosomes
p4 <- data.frame(result_id=rep(4), indv1=rep("A"), indv2=rep("E"),
                 chromosome=c(1, 2, 3, 4, 4, 4),
                 bp_start=c(0, 1, c[3]*(4/5), 0, c[4]/5, c[4]*(5/6)),
                 bp_end=c(c[1]/2, c[2]/5, c[3], c[4]/5-2e5, c[4]/5+1e6, c[4]),
                 length=runif(6), stringsAsFactors=FALSE)


df <- rbind(p1, p2, p3, p4)
df <- data.frame(id=seq(dim(df)[1]), df)
test_ibd_segments <- df
devtools::use_data(test_ibd_segments, chrome_bp, overwrite=TRUE)
