library(tidyr)
library(dplyr)

# Process each chr*_files list in ../filelist/
fileNames <- Sys.glob("../filelist/chr*_files")

for (fileName in fileNames) {
  chrlist <- read.table(fileName, header = FALSE, stringsAsFactors = FALSE) %>%
    separate(V1, into = c("chr", "prefix", "start", "end", "suffix"), sep = "([\\.\\-])") %>%
    arrange(as.integer(end)) %>%
    unite(col1, chr, prefix, start, sep = ".") %>%
    unite(col2, col1, end, sep = "-") %>%
    unite(ordered_id, col2, suffix, sep = ".")

  write.table(chrlist, file = paste0(fileName, ".txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
}

