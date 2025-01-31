library(tidyverse)
library(lubridate)
library(patchwork)
library(parallel)
options(mc.cores=parallel::detectCores()-4)

#--------------------------------------------------
# See which days have submatrices
#--------------------------------------------------

get_mm <- function(mm)
{
    get_ccvvdd <- function(f)
    {
        dat <- read_csv(f, show_col_types=F) %>%
               select(Country, Pango, Day) %>%
               unique
        return(dat)
    }

    smf <- Sys.glob(paste0("20241007/*/*/*/submat", mm, "-day1.csv"))
    dat <- bind_rows(lapply(smf, get_ccvvdd)) %>%
           mutate(!!mm := TRUE)
    return(dat)
}

dat30 <- get_mm("30")
dat50 <- get_mm("50")
dat75 <- get_mm("75")
dat100 <- get_mm("100")

dat <- full_join(dat30, dat50) %>% full_join(dat75) %>% full_join(dat100)
write_csv(dat, "have_submat.csv")

#--------------------------------------------------
# See which days with submatrices are the early/late dates, for out/non
#--------------------------------------------------
# nn = 3 earliest days and latest days with submatrices

get_dates <- function(dat, cc_vv, mm, nn=3)
{
    d1 <- dat %>%
          filter(Country == cc_vv[1], Pango == cc_vv[2]) %>%
          filter(!is.na(maxMM)) %>%
          filter(maxMM >= mm) %>%
          arrange(Day)

    if (nrow(d1) >= nn*2)
    {
        d2 <- d1 %>%
              slice(c(1:nn, tail(row_number(), nn))) %>%
              select(-maxMM) %>%
              add_column(MM = mm) %>%
              add_column(Type = c(rep("out", 3), rep("non", 3)))
    } else {
        d2 <- NULL
    }

    return(d2)
}

dat <- read_csv("freq_withsubmat.csv.xz")
cc_vv <- dat %>% select(Country, Pango) %>% unique

ans <- apply(cc_vv, 1, function(x) get_dates(dat, x, mm=30))
i_drop <- which(sapply(ans, is.null))
ans30 <- bind_rows(ans[-i_drop])

ans <- apply(cc_vv, 1, function(x) get_dates(dat, x, mm=50))
i_drop <- which(sapply(ans, is.null))
ans50 <- bind_rows(ans[-i_drop])

ans <- apply(cc_vv, 1, function(x) get_dates(dat, x, mm=75))
i_drop <- which(sapply(ans, is.null))
ans75 <- bind_rows(ans[-i_drop])

ans <- apply(cc_vv, 1, function(x) get_dates(dat, x, mm=100))
i_drop <- which(sapply(ans, is.null))
ans100 <- bind_rows(ans[-i_drop])

pca_dates <- bind_rows(ans30, ans50, ans75, ans100) %>%
             select(Country, Pango, Day, MM, Type)

write_csv(pca_dates, "pca_dates.csv")

#--------------------------------------------------
# Get the submatrices for the out/non dates
#--------------------------------------------------

pca_dates <- read_csv("pca_dates.csv") %>%
             mutate(File = paste0("20241007/",
                           Country, "/", Pango, "/submat", MM, "-day1.csv"))

for (mm in c(30, 50, 75, 100))
{
    pca_dates1 <- pca_dates %>% filter(MM == mm)
    dat1 <- mclapply(unique(pca_dates1$File), function(x) read_csv(x, show_col_types=F))
    dat2 <- bind_rows(dat1)

    dat3 <- dat2 %>% group_by(Country, Pango, Day) %>% sample_n(20)

    dat4 <- left_join(pca_dates1, dat3)
    write_csv(dat4, file=paste0("pca_", mm, ".csv.xz"))
}

#--------------------------------------------------
# Get submatrices for all dates between out and non
#--------------------------------------------------

pca_dates <- read_csv("pca_dates.csv")
pca_dates_lim <- pca_dates %>%
                 group_by(Country, Pango, MM) %>%
                 summarize(FirstDay = min(Day),
                           LastDay = max(Day))
have_submat <- read_csv("have_submat.csv")

mm <- 50 # repeat for 30 75 100
dat1 <- have_submat %>%
        filter(`50`) # NOTE: adjust for mm (yuck)
dat2 <- pca_dates_lim %>%
        filter(MM == mm)
getme <- left_join(dat1, dat2) %>%
         filter(Day >= FirstDay & Day <= LastDay) %>%
         arrange(Country, Pango, Day) %>%
         select(Country, Pango, Day)

files <- getme %>%
         mutate(File = paste0("20241007/",
                       Country, "/", Pango, "/submat", mm, "-day1.csv")) %>%
         distinct(File) %>%
         pull(File)
dat_sm <- mclapply(files, function(x) read_csv(x, show_col_types=F))
dat_sm <- bind_rows(dat_sm)

dat <- getme %>% left_join(dat_sm)
write_csv(dat, file=paste0("pca_", mm, "_full.csv.xz"))

#--------------------------------------------------
# Filter the variants
#--------------------------------------------------
# require the 6 submat to be within the training time range, and high enough max frequency

mm <- 50
dat_sm <- read_csv(paste0("pca_", mm, ".csv.xz"))

# Begin with Alpha
dat_sm %>% filter(Pango == "B.1.1.7") %>% group_by(Country) %>% summarize(min(Day))
# GBR     2020-11-05
# USA     2021-01-05
# End with BA.2
dat_sm %>% filter(Pango == "BA.2") %>% group_by(Country) %>% summarize(max(Day))
# GBR     2022-06-28
# USA     2022-07-11

# Max freq in the training time range, using only dates with enough data
freq <- read_csv("freq_withsubmat.csv.xz") %>%
        filter( (Day >= ymd("2020-11-05") & Day <= ymd("2022-06-28") & Country == "GBR") |
                (Day >= ymd("2021-01-05") & Day <= ymd("2022-07-11") & Country == "USA") )

# Record which Pango can be used to train, based on max freq
train1 <- freq %>%
          filter(!is.na(maxMM)) %>%
          group_by(Country, Pango) %>%
          summarize(MaxFreq = max(VarFreq)) %>%
          add_column(InDates = TRUE, .before = "MaxFreq")
train2 <- freq %>%
          select(Country, Pango) %>%
          distinct() %>%
          left_join(train1) %>%
          replace_na(list(InDates = FALSE))

# only keep pango that have both out and non, within the training time range
train3 <- dat_sm %>%
          filter( (Country == "GBR" & Type == "out" & Day >= ymd("2020-11-01")) |
                  (Country == "USA" & Type == "out" & Day >= ymd("2021-01-04")) |
                  (Country == "GBR" & Type == "non" & Day <= ymd("2022-07-01")) |
                  (Country == "USA" & Type == "non" & Day <= ymd("2022-07-15")) ) %>%
          group_by(Country, Pango) %>%
          summarize(OutNon=("out" %in% Type) & ("non" %in% Type)) %>%
          select(Country, Pango, OutNon) %>%
          distinct()

train <- full_join(train2, train3) %>%
         replace_na(list(OutNon = FALSE)) %>%
         mutate(map = OutNon & (MaxFreq >= 0.1))

write_csv(train, "pango_train.csv")

#--------------------------------------------------
# Form the PCA maps
#--------------------------------------------------

train <- read_csv("pango_train.csv")

for (mm in c(30, 50, 75, 100))
{
    dat0 <- read_csv(paste0("pca_", mm, ".csv.xz"))
    dat <- filter(train, map) %>%
           inner_join(dat0)

    pca <- prcomp(select(dat, starts_with("sym")), scale=T)

    # whether PC1 and PC2 should be sign-reversed
    dat_pca <- dat[,c(1:5, 8:18)]
    dat_pca$pc1 <- pca$x[,1]
    dat_pca$pc2 <- pca$x[,2]
    neg1 <- ((dat_pca %>% filter(Type == "out") %>% pull(pc1) %>% mean) < 0)
    neg2 <- ((dat_pca %>% filter(Type == "out") %>% pull(pc2) %>% mean) < 0)
    sign_reverse <- c(neg1, neg2)

    save(pca, dat_pca, sign_reverse, file=paste0("pca", mm, "_map.rda"))
}

mm <- 50
load(paste0("pca", mm, "_map", ".rda"))
(summary(pca)$importance)[-1, c("PC1", "PC2")]

#--------------------------------------------------
# Report PC1 for all variants
#--------------------------------------------------

for (mm in c(30, 50, 75, 100))
{
    load(paste0("pca", mm, "_map", ".rda"))

    dat_sm <- read_csv(paste0("pca_", mm, "_full.csv.xz"), show_col_types=F)
    pca_ans <- predict(pca, select(dat_sm, starts_with("sym")))
    dat <- cbind(dat_sm, pca_ans) %>%
           select(Country, Pango, Day, matmean, matvar, trimean, trivar, nucdiv, numseg, tajD, PC1, PC2) %>%
           group_by(Country, Pango, Day) %>% slice(1:20) %>% ungroup() %>%
           tibble()
    if (sign_reverse[1])
        dat$PC1 <- -dat$PC1
    if (sign_reverse[2])
        dat$PC2 <- -dat$PC2

    # separate values for submatrices
    write_csv(dat, paste0("pc1_", mm, "_map", ".csv.xz"))
}

#--------------------------------------------------
# PCA on non-normalized submatrices
#--------------------------------------------------

train <- read_csv("pango_train.csv")

for (mm in c(30, 50, 75, 100))
{
    dat0 <- read_csv(paste0("pca_", mm, ".csv.xz"))
    dat <- filter(train, map) %>%
           inner_join(dat0)

    pca <- prcomp(select(dat, starts_with("raw")), scale=T)

    # whether PC1 and PC2 should be sign-reversed
    dat_pca <- dat[,c(1:5, 7:18)]
    dat_pca$pc1 <- pca$x[,1]
    dat_pca$pc2 <- pca$x[,2]
    neg1 <- ((dat_pca %>% filter(Type == "out") %>% pull(pc1) %>% mean) < 0)
    neg2 <- ((dat_pca %>% filter(Type == "out") %>% pull(pc2) %>% mean) < 0)
    sign_reverse <- c(neg1, neg2)

    save(pca, dat_pca, sign_reverse, file=paste0("pca", mm, "_raw_map.rda"))
}
