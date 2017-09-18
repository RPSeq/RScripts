#!/usr/bin/Rscript

#################
### LIBRARIES ###
#################
library(ggplot2)
library(dplyr)
library(reshape2)
library(scales)
library(grid)
library(gridExtra)
##################
### /LIBRARIES ###
##################

#########################
### UTILITY FUNCTIONS ###
#########################
# Loads files of given ext in dir
load_dir <- function(dir, ext) {
  files <- list.files(dir, pattern = ext)
  DF <- NULL

  for (f in files) {
    dat <- read.table(paste(dir,f, sep = ""), header = F)
    dat$file <- unlist(strsplit(f,split = ".",fixed = T))[1]
    DF <- rbind(DF, dat)
  }

  return(DF)
}

# Wrapper virtually sets defaults for ggsave
saveplot <- function(outname, plot, h = 12, w = 18, d = 300) {
  return(
    ggsave(outname,
           plot,
           height = h,
           width  = w,
           dpi    = d)
  )
}

# Sorts df by GT with given order
gt_sort <- function(data, sorted) {
  return(
    within(data, GT
           <- factor(GT,
                     levels = sorted))
  )
}

# Sorts df by Sample with given order
sample_sort <- function(data, sorted) {
  return(
    within(data, Sample
           <- factor(Sample,
                     levels = sorted))
  )
}

# Replicates a control dataset for each expt sample
rep_control <- function(expts, control, name) {
  Reduce(rbind, lapply(expts,
                       function(sample, data)
                         mutate(data, Sample = sample),
                       data = control)) %>% mutate(Source = name)
}
##########################
### /UTILITY FUNCTIONS ###
##########################
# load mutant call data
input <- read.table("~/Desktop/SCNT_Caller/BLASTOS/run/pairwise/ALL5.txt", header = TRUE) %>% mutate(CASE = "BLASTOCYST")
input <- subset(input, Sample != "Old_bulk" & MGP==".")
levels(input$Sample) <- c("Bulk" ,"B1", "B2")

input <- melt(input)
colnames(input) <- c("ID", "Sample", "MGP", "CASE", "Replicate", "AB")


control <- read.table("~/Desktop/SCNT_Caller/BLASTOS/run/bulk_control.df", header = TRUE) %>% mutate(AB = ALTD/(ALTD+REFD))
levels(control$Sample) <- c("Bulk")

input <- input[c("Sample", "MGP", "AB", "CASE")]
control <- control[c("Sample", "MGP", "AB", "CASE")]

b <- rbind(input, control)


#######################
### SET COMMON VARS ###
#######################
# plot output dir
OUTDIR <- "~/Desktop/BLAST_PLOTS"


################################
### LOAD AND PREP GT/AB DATA ###
################################
#   read input
input <- "~/Desktop/SCNT_Caller/BLASTOS/run/germline_match.df"
match <- read.table(input, header = TRUE) %>% mutate(AB = ALTD / (REFD + ALTD))

#   set useful names for MGP presence
levels(match$MGP) <- c("NOVEL", "MGP")

##############################
### LOAD AND PREP COV DATA ###
##############################
# dir with histo files (*gen for genome avg hists)
covdir <- "~/Desktop/SCNT_Caller/BLASTOS/run/covsall/"

# load all .gen files in dir and subset cols
DF2 <- load_dir(covdir, "*.gen")[c("V2", "V3", "V5", "file")]

# informative colnames
colnames(DF2) <- c("COV", "COUNT", "FREQ", "Sample")

# collapse cov bins above 60 by summing FREQ and COUNT
DF2 <- rbind(
  subset(DF2, COV <= 60),
  group_by(subset(DF2, COV >= 60), Sample)
  %>% mutate(FREQ = sum(FREQ))
  %>% mutate(COUNT = sum(COUNT))
  %>% subset(COV == 60)
  %>% data.frame()
)

# limit to old samples
DF2 <- subset(DF2, substr(Sample, 1, 3) == "Old")

# trim age from sample name
DF2$Sample <- substr(DF2$Sample, 5, 9)

# define order
DF2 <- sample_sort(DF2, c("bulk",
                          "B1-1",
                          "B1-2",
                          "B2-1",
                          "B2-2"))

levels(DF2$Sample) <- c("Bulk",
                        "B1-1",
                        "B1-2",
                        "B2-1",
                        "B2-2")

# get blastocyts and set CASE var to "WGA"
BLASTS <- subset(DF2, Sample!="Bulk") %>% mutate(Source = "Blast")

blastnames <- as.character(unique(BLASTS$Sample))
# merge to add bulk dist to all sample panels
overlaid <- rbind(
  rep_control(expts = blastnames,
              control = subset(DF2, Sample == "Bulk"),
              name = "Bulk"),
  BLASTS)



######################
### PLOT FUNCTIONS ###
######################
titlesize = 22
XTITLESIZE = 20
YTITLESIZE = 20
TEXTSIZE = 16

#   Plot AB hists for the somatic mutation calls
mutant_ab_hists <- function(d, outname = NULL) {

  # set colors
  colors <- c("#0a3ece", "darkgrey")

  #make plot
  plot <- ggplot(d) +
    geom_histogram(aes(x = AB, fill = CASE, y = ..density..),
                   color = NA,

                   position = 'dodge',
                   binwidth = 0.033) +
    geom_vline(xintercept = 0.3,
               col = "red",
               linetype = "dashed") +

    labs(x = "VAF",
         y = "Percent",
         title = expression(paste(bold("(D)"), " Somatic VAFs"))
         ) +
    scale_x_continuous(limits = c(-0.05, 1.05), breaks = c(0, .25, .5, .75, 1)) +
    #     scale_y_continuous(limits = c(0, 6), expand=c(0,0), breaks = c(1,3,5)) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +

    theme(
      #       strip.text = element_blank(),
      #       strip.background = element_blank(),
      panel.spacing.y = unit(0.2, "lines"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", size = 1),
      legend.position = "none",
      #       legend.justification = "center",
      #       legend.position = c(0.2, 0.95),
      #       legend.direction = "vertical",
      #       legend.background = element_blank(),
      #       legend.key.size = unit(1, "lines"),
      #       legend.text = element_text(size = 12),
      #       legend.title = element_blank(),
      #       legend.margin = margin(t = unit(0, "lines")),
      #       legend.box.spacing = unit(0.3, "lines"),
      axis.title.x = element_text(colour = 'black', size = XTITLESIZE),
      axis.text.x = element_text(colour = 'black', size = TEXTSIZE),
      axis.ticks.x = element_line(colour = 'black'),
      axis.text.y = element_text(colour = 'black', size = TEXTSIZE),
      axis.title.y = element_text(colour = 'black', size = YTITLESIZE),
      plot.title = element_text(size = titlesize)
      # plot.margin = unit(c(2,0,0,0), "line")
    )

  if (!is.null(outname)) {
    saveplot(outname,
             plot,
             h = 3,
             w = 4,
             d = 300)
  }

  return(plot)
}

#   Plot AB hists for the het calls in each sample, MGP vs novel.
hists_ab <- function(d, outname = NULL) {

  # set colors
  colors <- c("#0a3ece", "darkgrey")

  # select old mice, MGP variants above 10 DP called het
  d <- subset(d, substr(Sample, 1, 3) == "Old")
  d <- subset(d, MGP=="MGP")
  d <- subset(d, ALTD+REFD>=10)

  # leave the errors in
#   d <- subset(d, GT == "0/1")

  # remove Old_ suffix (only plotting olds for grant)
  d$Sample <- substr(d$Sample, 5, 9)

  d <- sample_sort(d, c("bulk",
                        "B1-1",
                        "B1-2",
                        "B2-1",
                        "B2-2"))

  levels(d$Sample) <- c("Bulk",
                        "B1-1",
                        "B1-2",
                        "B2-1",
                        "B2-2")

  # get blastocyts and set CASE var to "WGA"
  BLASTS <- subset(d, Sample!="Bulk") %>% mutate(Source = "Blast")

  blastnames <- as.character(unique(BLASTS$Sample))

  over <- rbind(
    rep_control(expts = blastnames,
                control = subset(d, Sample == "Bulk"),
                name = "Bulk"),
    BLASTS)

  #make plot
  plot <- ggplot(over, aes(x = AB, fill = Source, y=(..density..)*(1.5))) +
    geom_histogram(position = 'dodge',
                   color = NA,
                   binwidth = 0.025) +
    facet_wrap(~ Sample,
               ncol = 1) +
    labs(x = "VAF",
         y = "Percent",
         title = expression(paste(bold("(B)"), " Germline VAFs"))
         ) +
    scale_x_continuous(limits = c(-0.05, 1.05), breaks = c(0, .25, .5, .75, 1)) +
#     scale_y_continuous(limits = c(0, 6), expand=c(0,0), breaks = c(1,3,5)) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +

    theme(
      strip.text = element_blank(),
      strip.background = element_blank(),
      panel.spacing.y = unit(0.2, "lines"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", size = 1),
      legend.justification = "center",
      legend.position = c(0.2, 0.95),
      legend.direction = "vertical",
      legend.background = element_blank(),
      legend.text = element_text(size = 16, face = "bold"),
#       legend.key.size = unit(2, "line"),
      legend.title = element_blank(),
      legend.margin = margin(t = unit(0, "lines")),
      legend.box.spacing = unit(0.3, "lines"),
      axis.title.x = element_text(colour = 'black', size = XTITLESIZE),
      axis.text.x = element_text(colour = 'black', size = TEXTSIZE),
      axis.ticks.x = element_line(colour = 'black'),
      axis.text.y = element_text(colour = 'black', size = TEXTSIZE),
      axis.title.y = element_blank(),
      plot.title = element_text(size = titlesize)
    )

  if (!is.null(outname)) {
    saveplot(outname,
             plot,
             h = 6.4,
             w = 2.5,
             d = 300)
  }

  return(plot)
}

#   Plot coverage dists
covplots <- function(d, outname = NULL) {

  # set colors
  colors <- c("#0a3ece", "darkgrey")

  #make plot
  plot <- ggplot(d, aes(x = COV, y = FREQ*100, fill = Source)) +
    geom_bar(stat = 'identity',
             position = 'dodge',
             width = 1) +
    geom_vline(xintercept = 10,
               color = "#e20b0b",
               linetype = "dashed") +
    geom_text(aes(label = Sample, x = 50, y = 5), size=6) +
    facet_wrap(~ Sample,
               ncol = 1) +
    labs(x = "X Coverage",
         y = "Percent",
         title = expression(paste(bold("(A)"), " Per-Base Covg."))
         ) +

    scale_x_discrete(limits = c(0, 10, 20, 30, 40, 50, 60)) +
    scale_y_continuous(limits = c(0, 6), expand=c(0,0), breaks = c(1,3,5)) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +

    theme(
      strip.text = element_blank(),
      strip.background = element_blank(),
      panel.spacing.y = unit(0.2, "lines"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", size = 1),
      legend.position = "none",
      axis.title.x = element_text(colour = 'black', size = XTITLESIZE),
      axis.text.x = element_text(size = TEXTSIZE, colour = c('black', "#e20b0b", 'black', 'black', 'black', 'black')),
      axis.ticks.x = element_line(colour = c('black', "#e20b0b", 'black', 'black', 'black', 'black')),
      axis.title.y = element_text(size = YTITLESIZE, colour = 'black'),
      axis.text.y = element_text(colour = 'black', size = TEXTSIZE),
      plot.title = element_text(size = titlesize)
    )

  if (!is.null(outname)) {
    saveplot(outname,
             plot,
             h = 6.4,
             w = 2.5,
             d = 300)
  }

  return(plot)
}

#   Plot GT concordance bars
bars_dept_NIH2_vert <- function(d, outname = NULL) {

  # set order of GT groups (affects legend order)
  d <- gt_sort(d, c("DP<10",
                    "./.",
                    "1/1",
                    "0/0",
                    "0/1"))

  # set corresponding colors
  colors <- c("#7C7D7A",
              "#e20b0b",
              "#ffcc66",
              "#1A1919",
              "#0a3ece")

  levels(d$GT) <- c("DP<10",
                    "NoCall",
                    "1/1",
                    "0/0",
                    "0/1")

  # set DP<10 gts
  d$GT[d$ALTD + d$REFD < 10] <- "DP<10"

  # filter for MGP vars only
  d <- subset(d, substr(Sample, 1, 5) != "Young" & MGP == "MGP")

  # remove Old_ suffix (only plotting olds for grant)
  d$Sample <- substr(d$Sample, 5, 9)

  d <- sample_sort(d, c("bulk",
                        "B1-1",
                        "B1-2",
                        "B2-1",
                        "B2-2"))

  levels(d$Sample) <- c("Bulk",
                        "B1-1",
                        "B1-2",
                        "B2-1",
                        "B2-2")
  mylabs <- c("DP<10  ",
              "NoCall  ",
              "1/1  ",
              "0/0  ",
              "0/1  ")

  # get counts and percents
  c <- count(d, Sample, GT) %>% mutate(SUM = sum(n), Percent = n/SUM)

  plot <- ggplot(c, aes(x = Sample, y = Percent)) +

    geom_bar(aes(col  = GT, fill = GT),
             stat = 'identity',
             width = 0.7) +

    labs(x     = "Sample",
         y = "Percent",
         title     = expression(paste(bold("(C)"), " Genotype Concordance"))
         ) +

    scale_fill_manual("Genotype", values  = colors, labels = mylabs) +
    scale_color_manual("Genotype", values = colors, labels = mylabs) +

    scale_y_continuous(expand = c(0,0),
                       labels = percent,
                       limits = c(0, 1.09),
                       breaks = c(.2, .4, .6, .8, 1.)) +

    scale_x_discrete(breaks = c("Bulk",
                                "B1-1",
                                "B1-2",
                                "B2-1",
                                "B2-2")) +

    theme(
      strip.text = element_blank(),
      strip.background = element_blank(),
#       panel.spacing.y = unit(0.2, "lines"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA, color="black", size = 1),
      axis.title.x = element_text(color="black", size = XTITLESIZE),
      axis.text.x = element_text(color="black", size = 24 ),
      axis.text.y = element_text(color = "black", size = TEXTSIZE),
      axis.title.y = element_blank(),
      legend.justification = "center",
      legend.position = c(0.5,0.965),
      legend.direction = "horizontal",
      legend.background = element_blank(),
      legend.text = element_text(size = 16, face = "bold"),
      legend.title = element_blank(),
      plot.title = element_text(size = titlesize)
    )



  if (!is.null(outname)) {
    saveplot(outname,
             plot,
             h = 6.4,
             w = 5,
             d = 300)

  }


  return(plot)
}


##################
### MAKE PLOTS ###
##################


# Mutant call AB hists
mut <- mutant_ab_hists(b)

# Hists of ABs corresponding to bars1.
# abname <- paste(OUTDIR, "abhists.png", sep = "/")
abs <- hists_ab(match)

# Per-base cov hists
# covsname <- paste(OUTDIR, "coverages2.png", sep = "/")
c <- covplots(overlaid)


# vname <- paste(OUTDIR, "BULK_MATCH_DEP10_NIH2_VERT.png", sep = "/")
vbarts <- bars_dept_NIH2_vert(match)

# set layout
lay <- rbind(c(1,2,3,3,4,4),
             c(1,2,3,3,NA,NA))

# set widths for each col
wids <- c(1.2,1.2,1,1,1,1)

# set heights for each row
hids <- c(1,1)

# arrange plots into layout
p <- arrangeGrob(c,
                 abs,
                 vbarts,
                 mut,
                 layout_matrix = lay,
                 widths = wids,
                 heights = hids
                 )

# save to PDF
mname <- paste(OUTDIR, "mergeall_DEPTH2.pdf", sep = "/")
saveplot(mname,
         p,
         h = 6,
         w = 17,
         d = 600)





# reload without factors
match <- read.table(input, header = TRUE, stringsAsFactors = FALSE) %>% mutate(AB = ALTD / (REFD + ALTD))
#plot correlation of ABs at bulk germline hets

# get bulk calls
bulks <- subset(match, Sample == "Old_bulk" & DP>10 & MGP=="1" & AB>0.30 & GT == "0/1")
