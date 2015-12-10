# read_path <- function(path, item)
#   read.delim(sprintf(path, item))
# 
# # MATS
# mats_path    <- "/genedata/NunoA/psi_files/mats_JC.%s.txt"
# mats_JC.SE   <- read_path(mats_path, "SE")
# mats_JC.A3SS <- read_path(mats_path, "A3SS")
# mats_JC.A5SS <- read_path(mats_path, "A5SS")
# mats_JC.MXE  <- read_path(mats_path, "MXE")
# mats_JC.RI   <- read_path(mats_path, "RI")
# 
# # VAST-TOOLS
# vast_all <- read.delim("/genedata/NunoA/psi_files/vast_all.tab")
# 
# # MISO
# miso_path <- "/genedata/NunoA/psi_files/miso_%s.miso_summary"
# miso_A3SS <- read_path(miso_path, "A3SS")
# miso_A5SS <- read_path(miso_path, "A5SS")
# miso_SE   <- read_path(miso_path, "SE")
# miso_AFE  <- read_path(miso_path, "AFE")
# miso_ALE  <- read_path(miso_path, "ALE")
# miso_MXE  <- read_path(miso_path, "MXE")
# miso_RI   <- read_path(miso_path, "RI")
# miso_TandemUTR <- read_path(miso_path, "TandemUTR")
# miso_isoforms  <- read_path(miso_path, "isoforms")
# 
# # MISO (annotation files)
# types <-
#   sprintf("/genedata/Resources/Annotations/MISO/hg19/%s.hg19.gff3",
#           c("AFE", "ALE", "SE", "MXE", "A5SS", "A3SS", "RI", "TandemUTR"))
# hg19.all <- NULL
# for (type in types) {
#   hg19.all <- rbind(hg19.all,
#                     read.delim(type, header = F, comment.char = "#"))
#   print(paste("MISO", type, "is ready"))
# }
# 
# # MISO event to UCSC genome browser
# ucsc <- function(df) {
#   write.table(format(df[c(1,4,5,3)]), row.names = F, col.names = F, quote = F)
# }
# 
# # SUPPA
# suppa_path   <- "/genedata/NunoA/psi_files/suppa_%s.tab.psi"
# suppa_A3.tab <- read_path(suppa_path, "A3")
# suppa_A5.tab <- read_path(suppa_path, "A5")
# suppa_AF.tab <- read_path(suppa_path, "AF")
# suppa_AL.tab <- read_path(suppa_path, "AL")
# suppa_MX.tab <- read_path(suppa_path, "MX")
# suppa_RI.tab <- read_path(suppa_path, "RI")
# suppa_SE.tab <- read_path(suppa_path, "SE")
#
# # Get SUPPA events as commands to run
# easy <- function(v) for (i in 1:length(v[1, ]))
#   cat(paste0(names(v)[i], ' = "', as.character(v[1, i]), '", '))