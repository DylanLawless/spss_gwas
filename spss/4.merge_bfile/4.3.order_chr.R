library(tidyr)
library(dplyr)

# Define input list directory and list files
fileNames <- Sys.glob("./LIST.mnt.data1.lawless.spss.Samira.SPSS.data.filelist/*_files")

for (fileName in fileNames) {
  chrlist <- read.table(fileName, header = FALSE, stringsAsFactors = FALSE) %>%
    separate(col = V1, into = c("chr", "prefix", "pos", "end", "suffix"), sep = "([\\.\\-])") %>%
    arrange(as.integer(end)) %>%
    unite(tmp1, chr, prefix, pos, sep = ".") %>%
    unite(tmp2, tmp1, end, sep = "-") %>%
    unite(sorted_id, tmp2, suffix, sep = ".")

  write.table(chrlist, file = paste0(fileName, ".txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
}

