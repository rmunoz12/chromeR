# chrome.R 
# Graphing chromosome IBD segments
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

library(ggplot2)
library(dplyr)
library(tidyr)

THRESHOLD <- 0

# Chromosome lengths per
# https://en.wikipedia.org/wiki/Human_genome
chrome_bp = c(249250621, 243199373, 198022430, 191154276, 180915260, 
              171115067, 159138663, 146364022, 141213431, 135534747,
              135006516, 133851895, 115169878, 107349540, 102531392, 
              90354753,  81195210,  78077248,  59128983,  63025520, 
              48129895,  51304566)


# Based on
# http://stackoverflow.com/questions/25995257/r-shift-values-in-single-column-of-dataframe-up
shift <- function(v, n) {
    # vector (i.e. column of dataframe)
    # n number to shift by:
    # negative shifts to lower row numbers
    # postive shifts to higher row numbers
    if (length(v) == 0) {
        return(v)
    }
    if (n < 0) {
        n <- -1 * n
        c(v[-(seq(n))], rep(NA, n))
    } else {
        c(rep(NA, n), v[1:(length(v) - n)])
    }
}


calc_name <- function(d, user) {
    d$name <- ""
    d$name[d$indv2 == user] <- d$indv1[d$indv2 == user]
    d$name[d$indv1 == user] <- d$indv2[d$indv1 == user]
    return(d)
}


get_gap_ahead <- function(bp_end, b_next, n, n_next, c, c_next, cend) {
    if (!is.na(n_next) && n == n_next) {
        if (!is.na(c_next) && c == c_next) {
            b_next - bp_end
        } else if (bp_end < cend) {
            cend - bp_end
        } else {
            NA
        }
    } else {
        if (bp_end < cend) {
            cend - bp_end
        } else {
            NA
        }
    }
}


get_gap_behind <- function(bp_start, n, n_prev, c, c_prev) {
    if ((!is.na(n_prev) && n == n_prev) ||
        is.na(n_prev)) {
        if ((!is.na(c_prev) && c != c_prev) || 
            is.na(c_prev)) {
            bp_start
        } else {
            NA
        }
    } else {
        NA
    }
}


fill_missing_chromes <- function(d) {
    n <- d %>% group_by(name) %>% summarise()
    c <- data.frame(name=rep(n$name, each=22), chromosome=1:22,
                    stringsAsFactors=FALSE)
    x <- d %>% 
            group_by(name, chromosome) %>% 
            summarise(count=length(chromosome))
    c <- c %>% left_join(x, by=c("name", "chromosome"))
    
    missing <- which(is.na(c$count))
    if (length(missing > 0)) {
        new_rows <- data.frame(name=character(0), 
                               chromosome=integer(0),
                               bp_start=logical(0), 
                               type=character(0),
                               len=numeric(0),
                               stringsAsFactors=FALSE)
        for (j in 1:length(missing)) {
            new <- data.frame(name=c$name[missing[j]], 
                              chromosome=c$chromosome[missing[j]],
                              bp_start=0, 
                              type="gap_ahead",
                              len=chrome_bp[c$chromosome[missing[j]]],
                              stringsAsFactors=FALSE)
            new_rows <- rbind(new_rows, new)
        }
        d <- rbind(d, new_rows)
    }
    return(d)
}


chrome_map <- function(d) {
    d <- d[c("result_id", "name", "chromosome", "cend", "bp_start", "bp_end")]
    d <- arrange(d, result_id, chromosome, bp_start)
    
    d$n_prev <- shift(d$name, 1)
    d$n_next <- shift(d$name, -1)
    d$c_prev <- shift(d$chromosome, 1)
    d$c_next <- shift(d$chromosome, -1)
    d$b_next <- shift(d$bp_start, -1)
    
    d$gap_ahead <- mapply(d$bp_end, d$b_next, d$name, d$n_next, 
                          d$chromosome, d$c_next, d$cend,
                          FUN=get_gap_ahead)
    d$gap_behind <- mapply(d$bp_start, d$name, d$n_prev, 
                          d$chromosome, d$c_prev,
                          FUN=get_gap_behind)
    d$share <- mapply(d$bp_end, d$bp_start,
                      FUN=function(end, start) end - start)
    
    d <- d[c("name", "chromosome", "bp_start", "gap_ahead", "gap_behind", "share")]
    m <- d %>% gather(type, len, -name, -chromosome, -bp_start)
    m <- na.omit(m)
    m$type <- factor(m$type, levels=c("gap_behind", "share", "gap_ahead"))
    m <- arrange(m, name, chromosome, bp_start, type)
    
    m$c_next <- shift(m$chromosome, -1)
    m$type <- as.character(m$type)
    m$type[m$type == "gap_ahead" & 
               m$chromosome == m$c_next &
               m$len <= THRESHOLD] <- "merge"
    m$type <- factor(m$type, levels=c("gap_ahead", "gap_behind", "share", "merge"))
    m$c_next <- NULL
    
    m <- fill_missing_chromes(m)
    return(m)
}


plot_chromes <- function(d, title, user, names) {
    d <- filter(d, indv1 == user || indv2 == user)
    d <- calc_name(d, user)
    d <- filter(d, name %in% names)
    d$cend = chrome_bp[d$chromosome]
    f <- chrome_map(d)
    levels(f$type) <- c("Not shared", "Not shared", "Shared", "Merge")
    p <- ggplot(f, aes(x=chromosome, y=len, fill=type)) +
        geom_bar(stat="identity", alpha=.95) + 
        coord_flip() + 
        scale_x_reverse() +
        scale_y_continuous("base pair") +
        scale_fill_manual("", values=c("#e9a3c9", "#a1d76a", "#ffffb3")) +
        theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
        ggtitle(title)
    if (length(names) > 1) {
        p <- p + facet_wrap(~name)
    }
    return(p)
}



