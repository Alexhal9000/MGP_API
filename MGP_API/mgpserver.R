# Capture the arguments from the terminal
args <- commandArgs(trailingOnly = TRUE)
# Set a default port if none is provided
port <- ifelse(length(args) > 0, as.numeric(args[1]), 6234)
plumber::plumb(file='~/Documents/MGP/jovid/MGP_API/plumber.R')$run(host="0.0.0.0", port=port)