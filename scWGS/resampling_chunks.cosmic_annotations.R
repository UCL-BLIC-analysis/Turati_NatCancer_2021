## ---- cosmic_annotations_v90
datatable(NULL)
cosmic_annotations <- tibble()
cosmic <- read.csv("../data/cosmic/grch37/cancer_gene_census.grch37.v90.csv", as.is = T)
cosmic.loc <- strsplit(cosmic$Genome.Location, ":")
cosmic$chr <- paste0("chr", sapply(cosmic.loc, function(x) {x[[1]]}))
cosmic$start <- as.integer(sapply(cosmic.loc, function(x) {strsplit(x[2], "-")[[1]][1]}))
cosmic$end <- as.integer(sapply(cosmic.loc, function(x) {strsplit(x[2], "-")[[1]][2]}))
for (i in 1:nrow(cnv_table)) {
  this_cnv <- cnv_table[i, , drop = FALSE]
  tag <- paste0("cosmic_", this_cnv$chr, "_", this_cnv$start, "_", this_cnv$end)
  location <- paste0(sub("chr", "", this_cnv$chr), ":", this_cnv$start, "-", this_cnv$end)
  label <- paste0(this_cnv$chr, ":", round(this_cnv$start / 1000000), "-", round(this_cnv$end / 1000000), "M")
  if (!is.null(opts_knit$get("output.dir"))) {
    cat(knit_child(text = paste0("\n## ", label, "\n\n"), quiet = T))
  }
  # cat(paste0('<div id="cosmic-', tag, '" class="section level2"><h2>', label, '</h2>'))
  cosmic.strict_matches <- cosmic %>%
    filter(chr == this_cnv$chr) %>%
    filter(start < this_cnv$end) %>%
    filter(end > this_cnv$start)
  cosmic.wider_matches <- cosmic %>%
    filter(chr == this_cnv$chr) %>%
    filter(start < this_cnv$end + 20000) %>%
    filter(end > this_cnv$start - 20000)
  print(htmltools::tagList(
    htmltools::h4("Cosmic"),
    datatable(elementId = label,
              cosmic.wider_matches,
              class = 'display',
              height = "auto",
              width = 800,
              options = list(columnDefs = list(
                list(width = '80px', targets = list(0,1,3,5,6,7,19,20,21)),
                list(width = '150px', targets = list(2,4,8,9,10,11,12,13,14,15,16,17,18))
              ),
              scrollX = TRUE),
              escape = FALSE),
    htmltools::h4("Ensembl"),
    htmltools::a(href = paste0("http://grch37.ensembl.org/Homo_sapiens/Location/Overview?r=", location),
                 target = "_blank",
                 paste("Link to Ensembl GRCh37 for", location))
  ))
  label <- sub(":0-0M", "", label, fixed = T)
  label <- sub("WGD.+", "WGD", label)
  label <- sub("HyDp.+", "HyDp", label)
  # cat('</div>')
  cosmic_annotations <- rbind(
    cosmic_annotations,
    tibble(lesion = label,
           genes = paste(cosmic.wider_matches %>% pull(Gene.Symbol), collapse = "/"),
           leukaemia_genes = paste(cosmic.wider_matches %>%
                                     filter(grepl("leukaemia", Tumour.Types.Germline.) |
                                              grepl("AML", Tumour.Types.Germline.) |
                                              grepl("ALL", Tumour.Types.Germline.) |
                                              grepl("amplified in other cancers", Tumour.Types.Germline.) |
                                              grepl("multiple other tumour types", Tumour.Types.Germline.) |
                                              grepl("other tumour types", Tumour.Types.Germline.) |
                                              grepl("leukaemia", Tumour.Types.Somatic.) |
                                              grepl("AML", Tumour.Types.Somatic.) |
                                              grepl("ALL", Tumour.Types.Somatic.) |
                                              grepl("amplified in other cancers", Tumour.Types.Somatic.) |
                                              grepl("multiple other tumour types", Tumour.Types.Somatic.) |
                                              grepl("other tumour types", Tumour.Types.Somatic.)
                                            ) %>%
                                     pull(Gene.Symbol), collapse = "/"))
  )
}
