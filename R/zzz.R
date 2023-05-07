#' Check TMB and Matrix versions
#' 
#' Helper function to provide instructions if there is a mis-match
#' @param silent Default is TRUE. If FALSE, then gives the full instructions
#' @export
marssTMB_CheckPackageVersions <- function(silent=TRUE) {
  file1 <- paste0(system.file(package="marssTMB"),"/Matrix-version")
  file2 <- paste0(system.file(package="marssTMB"),"/TMB-version")
  cur.Matrix.version <- as.character(utils::packageVersion("Matrix"))
  cur.TMB.version <- as.character(utils::packageVersion("TMB"))
  if(!file.exists(file1)) {
    writeLines(cur.Matrix.version, con = file1)
  }
  if(!file.exists(file2)) {
    writeLines(cur.TMB.version, con = file2)
  }
  marssTMB.Matrix.version <- readLines(file1)
  marssTMB.TMB.version <- readLines(file2)
  if(!identical(marssTMB.Matrix.version, cur.Matrix.version) |
     !identical(marssTMB.TMB.version, cur.TMB.version)) {
    if(silent){
      warning(
        "marssTMB was built with TMB and Matrix versions ",
        marssTMB.TMB.version, " and ", marssTMB.Matrix.version,
        "\n",
        "TMB and Matrix versions installed on your computer are ",
        cur.TMB.version, " and ", cur.Matrix.version, "\n",
        "If you run into problems, run marssTMB_CheckPackageVersions(silent = FALSE) for instructions.\n", call.=FALSE)
    }else{
      msg <- c(
      "marssTMB was built with TMB and Matrix versions ",
      marssTMB.TMB.version, " and ", marssTMB.Matrix.version,
      "\n",
      "TMB and Matrix versions installed on your computer are ",
      cur.TMB.version, " and ", cur.Matrix.version,
      "\n",
      "What to do:\n",
      "0) Ignore the warning until TMB errors appear.",
      "1) Update your TMB and Matrix versions (and potentially R)\n",
      "  Look at the R version of marssTMB on https://atsa-es.r-universe.dev/marssTMB.\n",
      "  Update R if needed, i.e. your version of R is behind those listed.\n",
      "  update.packages()\n",
      "2) Install from source to create binaries that are built with the TMB and Matrix versions on your computer. This solution requires that you can build packages from source. Instructions\n",
      "  Run library(TMB)\n",
      "  If you get a warning about TMB and Matrix versions being out of sync, then run install.packages('TMB', type = 'source').\n",
      "  Now run install.packages('marssTMB', type = 'source').\n",
      "3) Update TMB and Matrix to match marssTMB.\n",
      "  require(devtools)\n",
      " install_version('TMB', version = 'x.x.x', repos = 'http://cran.us.r-project.org')\n",
      " install_version('TMB', version = 'x.x.x', repos = 'http://cran.us.r-project.org')\n",
      " in the code above, replace 'x.x.x' with the versions used in marssTMB.")
    cat(msg)
    }
  }
}

.onLoad <- function(lib, pkg) {
  marssTMB_CheckPackageVersions()
}


