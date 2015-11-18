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

# library(ggplot2)
# library(dplyr)
# library(tidyr)

# Based on
# http://stackoverflow.com/questions/25995257/r-shift-values-in-single-column-of-dataframe-up
shift <- function(v, n) {
    # v vector (i.e. column of dataframe)
    # n number to shift by:
    # negative shifts to lower row numbers
    # postive shifts to higher row numbers
    v_len <- length(v)
    if (v_len < abs(n)) {
        stop("Shift length greater than vector length.")
    }
    if (v_len == 0) {
        return(v)
    }
    if (n < 0) {
        n <- -1 * n
        c(v[-(seq(n))], rep(NA, n))
    } else {
        if (v_len == n) {
            rep(NA, n)
        } else {
            c(rep(NA, n), v[1:(v_len - n)])
        }
    }
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
        if (bp_start > 0) {
            bp_start - 0
        } else {
            NA
        }
    }
}


fill_missing_chromes <- function(d) {
    # Takes input from chrome_map with columns:
    # 1: key
    # 2: chromosome
    # 3: bp_start
    # 4: type
    # 5: len

    n <- d %>% dplyr::group_by_(names(d)[1]) %>% dplyr::summarise()
    c <- data.frame(key=rep(n[[1]], each=22),
                    chromosome=1:22,
                    stringsAsFactors=FALSE)
    y <- d %>%
            dplyr::group_by_(names(d)[1], names(d)[2]) %>%
            dplyr::summarise_(count=length(names(d)[2]))
    c <- merge(c, y, by.x=names(c)[1:2], by.y=names(y)[1:2], all.x=TRUE)
    c <- subset(c, is.na(c$count))
    if (dim(c)[1] > 0) {
        c$count <- NULL
        c <- cbind(c,
                   data.frame(0, "gap_ahead",
                              chrome_bp[c[ , names(c)[2]]],
                              stringsAsFactors=FALSE))
        names(c) <- names(d)
        d <- rbind(d, c)
    }
    return(d)
}


chrome_map_helper <- function(d, threshold) {
    # Expects the following column order, which is passed
    # by build_chrome_map():
    #
    # 1: key
    # 2: chromosome
    # 3: bp_start
    # 4: bp_end

    orig_cols <- names(d)
    names(d) <- make.names(seq(ncol(d)))
    d <- d[order(d[1], d[2], d[3]), ]
    row.names(d) <- NULL

    d <- cbind(d, cend=chrome_bp[d[, 2]]) # 5

    d$n_prev <- shift(d[ , 1],  1)  # 6
    d$n_next <- shift(d[ , 1], -1)  # 7
    d$c_prev <- shift(d[ , 2],  1)  # 8
    d$c_next <- shift(d[ , 2], -1)  # 9
    d$b_next <- shift(d[ , 3], -1)  # 10

    d$gap_ahead <- mapply(d[ , 4], d$b_next, d[ , 1], d$n_next,
                          d[ , 2], d$c_next, d$cend,
                          FUN=get_gap_ahead)                 # 11
    d$gap_behind <- mapply(d[ , 3], d[ , 1], d$n_prev,
                          d[ , 2], d$c_prev,
                          FUN=get_gap_behind)                # 12
    d$share <- mapply(d[ , 4], d[ , 3],
                      FUN=function(end, start) end - start)  # 13

    d <- d[c(1:3, 11:13)]
    cols <- names(d)
    m <- d %>% tidyr::gather_("type", "len", colnames(d)[4:6])
    m <- na.omit(m)
    m$type <- factor(m$type, levels=c("gap_behind", "share", "gap_ahead"))
    m <- dplyr::arrange_(m, cols[1], cols[2], cols[3], "type")

    m$c_next <- shift(m[, cols[2]], -1)
    m$type <- as.character(m$type)
    m$type[m$type == "gap_ahead" &
               m[ , cols[2]] == m$c_next &
               m$len <= threshold] <- "merge"
    m$type <- factor(m$type, levels=c("gap_ahead", "gap_behind", "share", "merge"))
    m$c_next <- NULL

    m <- fill_missing_chromes(m)
    names(m)[1:3] <- orig_cols[1:3]
    return(m)
}

#' Calculate shared and not shared lengths.
#'
#' \code{chrome_map} returns a dataframe with shared and not shared regions.
#'
#' Based on a set of IBD segments, \code{chrome_map} calculate the lengths (in
#' basepairs) of the share segments as well as any unshared areas in each
#' chromosome (i.e., gaps). In addition, the \code{threshold} parameter can be
#' used to identifiy areas that are small (which \code{chrome_map} will then
#' classify as "merge").
#'
#' @export
#'
#' @param data Dataframe with \code{key}, \code{chromosome}, \code{bp_start},
#'      and \code{bp_end} columns. Each row describes an IBD segment.
#' @param key String identifying the column in \code{data} that uniquely
#'      identifies the pair of individuals corresponding to each IBD segment.
#' @param chromosome String indentifying the column in \code{data} that
#'      indicates the chromosome of each IBD segment.
#' @param bp_start String indentifying the column in \code{data} that
#'      indicates the starting basepair of each IBD segment.
#' @param bp_end String indentifying the column in \code{data} that
#'      indicates the ending basepair of each IBD segment.
#' @param threshold Numeric value, greater than or equal to 0 for which a gap
#'      should be classified as a "merge."
#'
#' @examples
#' data(chromeR)
#' chrome_map(test_ibd_segments, "result_id", "chromosome", "bp_start", "bp_end")
chrome_map <- function(data, key, chromosome, bp_start, bp_end,
                       threshold=2e5) {
    cols <- c(key, chromosome, bp_start, bp_end)
    data <- data[cols]
    return(chrome_map_helper(data, threshold))
}

#' IBD segments plot
#'
#' Plots the output of \code{chrome_map}.
#'
#' @export
#'
#' @param d Dataframe output from \code{chrome_map}.
#' @param key String that identifies the column in \code{d} that uniquely
#'      identifies the pair of individuals corresponding to each row.
#' @param title String for the plot's title.
#'
#' @examples
#' data(chromeR)
#' df <- chrome_map(test_ibd_segments, "result_id", "chromosome", "bp_start", "bp_end")
#' chrome_plot(df, "result_id", "Test IBD Matches")
chrome_plot <- function(d, key, title) {
    levels(d$type) <- c("Not shared", "Not shared", "Shared", "Merge")
    p <- ggplot2::ggplot(d, ggplot2::aes(x=chromosome, y=len, fill=type)) +
        ggplot2::geom_bar(stat="identity", alpha=.95) +
        ggplot2::coord_flip() +
        ggplot2::scale_x_reverse() +
        ggplot2::scale_y_continuous("base pair") +
        ggplot2::scale_fill_manual("", values=c("#e9a3c9", "#a1d76a", "#ffffb3")) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90, vjust=0.5, hjust=1)) +
        ggplot2::ggtitle(title) +
        ggplot2::facet_wrap(as.formula(paste0("~", key)))
    return(p)
}
