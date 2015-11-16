# Chromosome lengths per
# https://en.wikipedia.org/wiki/Human_genome
c <- c(249250621, 243199373, 198022430, 191154276, 180915260, 
              171115067, 159138663, 146364022, 141213431, 135534747,
              135006516, 133851895, 115169878, 107349540, 102531392, 
              90354753,  81195210,  78077248,  59128983,  63025520, 
              48129895,  51304566)

# Input dataframe columns
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
                 length=runif(1), stringsAsFactors=TRUE)

df <- rbind(p1, p2)
df <- data.frame(id=seq(dim(df)[1]), df)
rm(c, p1, p2)