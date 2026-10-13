################################################################################
# Runs the full workflow in order, in a single R session.
# Usage (working directory = repository root):
#   Rscript run_all.R        or, from R/RStudio:  source("run_all.R")
# Tables go to Outputs/, figures to Plots/, session information to Outputs/sessionInfo.txt
################################################################################

.scripts = c("1_FUNCTIONS.R", "2_DATA.R", "3_ANALYSIS.R", "4_VISUALIZATIONS.R")
if (!all(file.exists(file.path("Scripts", .scripts)))) {
  stop("Set the working directory to the repository root (the folder containing Scripts/ and Data/).")
}
dir.create("Outputs", showWarnings = FALSE)
dir.create("Plots", showWarnings = FALSE)

# Hidden (dot) names survive the rm(list = ls()) in 1_FUNCTIONS.R
for (.s in .scripts) {
  message("Running ", .s, " ...")
  .t0 = Sys.time()
  source(file.path("Scripts", .s))
  message(.s, " finished in ", format(round(Sys.time() - .t0, 1)))
}

writeLines(capture.output(sessionInfo()), "Outputs/sessionInfo.txt")
