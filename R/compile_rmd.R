## ADB


## Compile Rmarkdown report

# compile_rmd function ----------------------------------------------------
compile_rmd <- function(file) {
  rmarkdown::render(
    input = here::here("rmds", paste(file, ".rmd", sep = "")),
    output_file = here::here("reports", paste(file, ".html", sep = ""))
  )
}

compile_rmd_pdf <- function(file) {
  rmarkdown::render(
    input = here::here("rmds", paste(file, ".rmd", sep = "")),
    output_file = here::here("reports", paste(file, ".pdf", sep = ""))
  )
}


# quarto function ----------------------------------------------------
compile_quarto <- function(file, output_format = 'html') {
  quarto::quarto_render(input =       here::here('rmds', paste(file, '.qmd', sep ='')),
                        output_format = output_format)
                         # output_file = here::here('reports', paste(file, '.html', sep = '') ))

}
