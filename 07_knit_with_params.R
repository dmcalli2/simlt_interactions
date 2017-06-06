#07 knit with parameters

for(i in 1:6) {
  rmarkdown::render("05_Summarise_results.Rmd", params = list(
  region = i), output_file = paste0("results/text_results", i, ".html"))
}
