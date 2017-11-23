#07 knit with parameters

for(i in 1:6) {
  print(i)
  rmarkdown::render("05_Summarise_results.Rmd", params = list(
  effect_number = i), output_file = paste0("results/text_results", i, ".html"))
}
