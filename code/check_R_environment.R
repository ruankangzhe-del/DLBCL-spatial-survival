args <- commandArgs(trailingOnly=TRUE)
expected <- read.delim(args[1], stringsAsFactors=FALSE)
cat('RUNTIME\n'); print(R.version.string); print(.libPaths())
if (as.character(getRversion()) != '4.5.2') stop('R_VERSION_MISMATCH_EXPECTED_4.5.2')
for (i in seq_len(nrow(expected))) {
  pkg <- expected$Package[i]
  actual <- as.character(packageVersion(pkg))
  cat(pkg, actual, '\n')
  if (package_version(actual) != package_version(expected$Version[i])) stop(paste('PACKAGE_VERSION_MISMATCH', pkg))
}
library(survival)
library(metafor)
print(sessionInfo())
cat('PASS_SCOPED_R_DEPENDENCIES\n')
