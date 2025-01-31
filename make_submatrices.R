require(readr)
require(lubridate)
require(dplyr)
require(ape)
require(pegas)
source("spectral.R") # https://gist.github.com/jlmelville/772060a26001d7d25d7453b0df4feff9

# Note: Metadata and sequence data should already be divided up into directories.
# Normal operations would run over all such directories, but just one example is shown here.
ff <- "20241007"
cc <- "GBR"
vv <- "B.1.177"
dd <- "2020-08-31"
filestem <- file.path(ff, cc, vv)

days_window <- 1                      # size of the sliding window; must be odd
submatrix_sizes <- c(30, 50, 75, 100) # we'll run all of these
n_submat <- 20                        # number of replicate submatrices

stopifnot(days_window %% 2 == 1)
day_range <- ymd(dd) + seq(days_window) - (days_window+1)/2
statsfile <- file.path(filestem, paste0("sumstat-day", days_window, ".csv"))

#--------------------------------------------------
# Metadata
#--------------------------------------------------
# Assemble the files that include info from GISAID and Nextstrain

ct <- cols(Name = col_character(),
           NameEPI = col_character(),
           Day = col_date(format = ""),
           Location = col_character(),
           VariantString = col_logical(),
           Country = col_character(),
           CountryCode = col_character(),
           Continent = col_character(),
           Pango = col_character(),
           VOC = col_logical(),
           Nextclade_clade = col_character(),
           Nextclade_pango = col_character(),
           qc.overallScore = col_double(),
           qc.overallStatus = col_character(),
           Nextclade_errors = col_character())

metafiles <- file.path(filestem, day_range, "metadata2.csv")
metafiles <- metafiles[file.exists(metafiles)]

meta <- bind_rows(lapply(metafiles, read_csv, col_types=ct)) %>%
        filter(qc.overallStatus %in% c("good", "mediocre")) %>%
        filter(is.na(Nextclade_errors))

# Unfortunately, sometimes seqs are duplicated in the metadata.
# If rows have the same Name, check that they have the same other relevant content.
#   If so, keep one row.  If not, discard the rows.
meta <- meta %>%
        distinct(Name, Day, CountryCode, Pango, VOC, Nextclade_clade, qc.overallStatus) %>%
        group_by(Name) %>%
        filter(n() == 1) %>%
        ungroup()

stopifnot(nrow(meta) > min(submatrix_sizes))

#--------------------------------------------------
# Alignment
#--------------------------------------------------

trim_aln <- function(aln)
{
    # drop beginning and ending sites

    drop_start <- 100
    drop_end <- 100

    keep_i <- seq((drop_start+1), (ncol(aln)-drop_end))
    aln <- aln[, keep_i]

    # drop low-quality sequences

    check_quality <- function(x)
    {
        # returns the proportion of the sequence that is "good" (a, c, g, t)
        bf <- base.freq(x, all=T)
        good <- sum(bf[c("a", "c", "g", "t")])
        return(good)
    }

    quality <- sapply(1:nrow(aln), function(i) check_quality(aln[i,]))
    keep_i <- which(quality >= 0.9)

    return(aln[keep_i, ])
}

# these were split out of the giant fasta downloaded from GISAID
alnfiles <- file.path(filestem, day_range, "seqs.fasta.gz")
alnfiles <- alnfiles[file.exists(alnfiles)]

aln <- do.call("rbind", aln)
aln <- aln[which(rownames(aln) %in% meta$Name),]

# double-check that there aren't any seqs in aln that are missing from meta
keepme <- intersect(rownames(aln), meta$Name)
aln <- aln[which(rownames(aln) %in% keepme),]
meta <- meta %>%
        filter(Name %in% keepme)
num_seqs <- length(keepme)

stopifnot(num_seqs == nrow(aln))

# record which seqs are actually available for use
if (days_window == 1)
{
    outfile <- file.path(ff, cc, vv, dd, "metadata3.csv")
    write_csv(meta, outfile)
}

#--------------------------------------------------
# Submatrices
#--------------------------------------------------

#  s = size of the submatrix
# al = alignment
get_submatrix <- function(s, al)
{
    # draw one focal sequence
    # (doesn't prevent it being used for another submatrix)
    i_focal <- sample(seq(nrow(al)), size=1)
    seq_focal <- al[i_focal,]

    # get its distance from all the other sequences
    dists <- sort(apply(al[-i_focal,], 1, function(x)
                    dist.dna(rbind(seq_focal, x), model="raw", pairwise.deletion=T)))

    # get the closest other sequences
    seq_names <- c(rownames(al)[i_focal], names(head(dists, n = s-1)))

    # make the distance matrix
    a <- al[seq_names,]
    dm <- dist.dna(a, model="raw", pairwise.deletion=T, as.matrix=T)

    return(list(distmat = dm, aln = a))
}

# sm = submatrix
# type = norm or raw
get_eigenvalues <- function(sm, type)
{
    WD <- gen_LMat(sm)
    Lmat <- lapm(WD)

    if (!(any(is.na(Lmat$L)) || any(is.na(Lmat$Lsym))))
    {
        if (type == "norm")
        {
            res <- eig(Lmat$Lsym, norm="n")
        } else if (type == "raw") {
            res <- eig(Lmat$L, norm="n")
        }
        ans <- res$values[-1]
    } else {
        ans <- NULL
    }

    return(ans)
}

submatrix_sizes <- submatrix_sizes[submatrix_sizes <= num_seqs]
stopifnot(length(submatrix_sizes) > 0)

for (s_submat in submatrix_sizes)
{
    ans <- list()
    n_success <- 0
    ev_list <- list()

    while (n_success < n_submat)
    {
        # make a sub-distance matrix and sub-alignment
        tmp <- get_submatrix(s_submat, aln)
        sm <- tmp$distmat
        sma <- tmp$aln

        # get stats for the submatrix/subalignment
        ev1 <- get_eigenvalues(sm, "norm")
        if (!is.null(ev1))
        {
            # eigenvalues
            names(ev1) <- paste0("sym", seq_along(ev1))
            ev2 <- get_eigenvalues(sm, "raw")
            names(ev2) <- paste0("raw", seq_along(ev2))

            # distance matrix mean and variance
            stats <- c("matmean" = mean(sm), "matvar" = var(as.vector(sm)))
            stats <- c(stats, "trimean" = mean(sm[lower.tri(sm)]), "trivar" = var(as.vector(sm[lower.tri(sm)])))

            # other popgen summary stats
            stats <- c(stats,
                       "nucdiv" = pegas::nuc.div(sma), # nucleotide diversity
                       "numseg" = length(ape::seg.sites(sma)), # number of segregating sites
                       "tajD" = pegas::tajima.test(sma)$D) # Tajima's D

            n_success <- n_success + 1
            ev_list[[n_success]] <- c(stats, ev1, ev2)
        }
    }

    # record the results
    ans <- data.frame(Country=cc, Pango=vv, Day=dd, do.call(rbind, ev_list))
    submatfile <- file.path(filestem, paste0("submat", s_submat, "-day", days_window, ".csv"))
    write_csv(ans, file=submatfile, append=file.exists(submatfile))
}
