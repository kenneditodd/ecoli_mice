# read tables
nash_meta <- read.delim2(file = "../../refs/nash_meta.tsv")
tgen_qc <- read.delim2(file = "../../refs/tgen_qc.tsv")

# reformat nash meta
colnames(nash_meta) <- tolower(colnames(nash_meta))
colnames(nash_meta) <- gsub("\\.+","_", colnames(nash_meta))
nash_meta <- nash_meta %>% rename(timepoint_days = timepoint_days_,
                                  age = group)

# add new group column
nash_meta$group <- paste0(
  nash_meta$treatment, nash_meta$timepoint_days, nash_meta$age, nash_meta$sex
)
nash_meta$group <- gsub("ecoli","E",nash_meta$group)
nash_meta$group <- gsub("saline","S",nash_meta$group)
nash_meta$group <- gsub("2","2D",nash_meta$group)
nash_meta$group <- gsub("18","18D",nash_meta$group)
nash_meta$group <- gsub("young","Y",nash_meta$group)
nash_meta$group <- gsub("old","O",nash_meta$group)
nash_meta$group <- gsub("Female","F",nash_meta$group)
nash_meta$group <- gsub("Male","M",nash_meta$group)

# add sample column
nash_meta$sample <- paste0(nash_meta$group, stringr::str_match(nash_meta$sample_id, "[oldyoung]+_[fesm]+_([0-9]+)")[,2])
