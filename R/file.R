#' Prepare path for a new output file.
#' @param outf characters, file name with desired path
#' @param remove bool, if remove old file, default TRUE
#' @return None
#' @export
cleanOutfile <- function(outf, remove = TRUE) {
  outdir <- dirname(outf)
  if (!dir.exists(outdir)) {
    message(paste(outdir, "does not exist and will be created."))
    dir.create(outdir, recursive = T)
  }
  if (file.exists(outf) & remove) {
    message(paste(outf, "exists and will be removed."))
    file.remove(outf)
  }
}

#' If file does not exists, then prepare it.
#' @param f characters, file name
#' @return None
#' @export
prepareOutfile <- function(f) {
  if (!file.exists(f)) {
    cleanOutfile(f, remove = TRUE)
  }
}

#' If dir does not exists, then create it.
#' @param d characters, dir name
#' @return None
#' @export
prepareOutdir <- function(d) {
  if (!dir.exists(d)) {
    message(paste(d, "does not exist and will be created."))
    dir.create(d, recursive = T)
  }
}



#' load RData file like readRDS
#' @param fileRData characters, RData file.
#' @param pattern characters, to check a RData file, default is "\\.RData"
#' @param check bool, if we need to check fileRData ends with .RData, default TRUE
#' @return R object
#' @export
loadRData <- function(fileRData, pattern = "\\.RData", check = TRUE) {
  if (check & (!grepl(pattern = pattern, x = fileRData))) {
    stop(fileRData, " is not RData based on the pattern ", pattern)
  }
  env <- new.env()
  nm <- load(fileRData, env)[1]
  env[[nm]]
}

#' Check if element in argument has no defination
#' argument generated from package optparse
#' But this function can be used for a list
#' @param args list-like object with attribute names
#' @return None, if no defined attribute in args, stop will be called
#' @export
checkArgsExistOrStop <- function(args) {
  invisible(lapply(names(args), function(v) {
    if (is.null(args[[v]])) {
      stop("Args have no attribute ", v)
    }
  }))
}

#' Check if file/dir exist, if not, then stop
#' @param f characters, a file or director
#' @return None, will send stop signal if
#' f doesn't exist.
#' @export
checkFileExistOrStop <- function(f) {
  if ((!file.exists(f)) & (!dir.exists(f))) {
    stop(f, " does not exist.")
  }
}
