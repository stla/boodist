dllunload <- function(){
  dyn.unload(
    system.file("libs", "x64", "boodist.dll", package = "boodist")
  )
}

myinstall <- function(restart = FALSE) {
  try(pkgload::unload("boodist"))
  try(dllunload())
  if(restart && rstudioapi::isAvailable()) {
    rstudioapi::restartSession(
      "devtools::install(quick = TRUE, keep_source = TRUE)"
    )
  } else {
    devtools::install(quick = TRUE, keep_source = TRUE)
  }
}

mydocument <- function(restart = FALSE) {
  if(restart && rstudioapi::isAvailable()) {
    rstudioapi::restartSession(
      "roxygen2::roxygenise(load_code = roxygen2::load_installed)"
    )
  } else {
    roxygen2::roxygenise(load_code = roxygen2::load_installed)
  }
}
