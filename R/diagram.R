
library(gridDiagram)

# Draw boxes --------------------------------------------------------------

svg("figures/diagram.svg", width = 9, height = 2)

grid.box("GEO series", sum(table(ds_redline$model)), 
         x = 1/7, y = 2/3,
         name = "all")

grid.box("GEO series", table(ds_redline$model)[[1]], 
         x = 2/7, y = 2/3,
         name = "human_or_mouse")

boxCurve("all",
         "human_or_mouse", 
         sum(table(ds_redline$model)), 
         table(ds_redline$model)[[1]])

grid.box("GEO series", c(nrow(suppfilenames_present),
                         "Suppl. files", 
                         nrow(suppfilenames_present_unnested)), 
         x = 3/7, y = unit(2/3, "npc") - unit(1, "lines"),
         foreign = "Suppl. files", name = "suppf_pres")

boxCurve("human_or_mouse", "suppf_pres",
         table(ds_redline$model)[[1]], 
         nrow(suppfilenames_present))

grid.box("GEO series", 
            c(n_distinct(suppfiles_notonly_rawtar_unnested$Accession),
            "Suppl. files", 
            nrow(suppfiles_notonly_rawtar_unnested)),
            x = 4/7, y = unit(2/3, "npc") - unit(1, "lines"),
            foreign = "Suppl. files", name = "notonly_raw")

boxCurve("suppf_pres", "notonly_raw",
         nrow(suppfilenames_present_unnested),
         nrow(suppfiles_notonly_rawtar_unnested))

grid.box("GEO series", c(n_distinct(suppfiles_of_interest$Accession),
                         "Suppl. files",
                         nrow(suppfiles_of_interest)), 
         x = 5/7, y = unit(2/3, "npc") - unit(1, "lines"),
         foreign = "Suppl. files", name = "suppf_oi")

boxCurve("notonly_raw", "suppf_oi",
         nrow(suppfiles_notonly_rawtar_unnested), 
         nrow(suppfiles_of_interest))

pvalue_sets <- p_values %>% 
  filter(features >= nrowthreshold)

grid.box("GEO series", c(n_distinct(pvalue_sets$Accession),
                         "Suppl. files",
                         n_distinct(pvalue_sets$suppfiles),
                         "P value sets",
                         nrow(pvalue_sets)), 
         x = 6/7, y = unit(2/3, "npc") - unit(2, "lines"), 
         foreign = "Suppl. files", name = "pvals")

boxCurve("suppf_oi", "pvals",
         nrow(suppfiles_of_interest),
         length(pvalue_sets$pvalues),
         c("right", "down"))


dev.off()

# Add text ----------------------------------------------------------------


