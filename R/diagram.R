
library(gridDiagram)

# Draw boxes --------------------------------------------------------------

grid.box("GEO series", sum(table(ds_redline$model)), 
         x = 1/3, y = 2/3,
         name = "all")

grid.box("GEO series", table(ds_redline$model)[[1]], 
         x = 2/3, y = 2/3,
         name = "human_or_mouse")

grid.box("GEO series", c(nrow(suppfilenames_present),
                         "Suppl. files", 
                         nrow(suppfilenames_present_unnested)), 
         x = 1/5, y = 1/3,
         foreign = "Suppl. files", name = "suppf_pres")

grid.box("GEO series", c(n_distinct(suppfiles_of_interest$Accession),
                         "Suppl. files",
                         nrow(suppfiles_of_interest)), 
         x = 1/2, y = 1/3,
         foreign = "Suppl. files", name = "suppf_oi")


pvalue_sets <- p_values %>% 
  filter(features >= nrowthreshold)

grid.box("GEO series", c(n_distinct(pvalue_sets$Accession),
                         "Suppl. files",
                         nrow(pvalue_sets),
                         "P value sets",
                         length(unlist(pvalue_sets$pvalues, recursive = FALSE))), 
         x = 4/5, y = unit(1/3, "npc") - unit(1, "lines"), 
         foreign = "Suppl. files", name = "pvals")


# Draw curves -------------------------------------------------------------


boxCurve("all","human_or_mouse", sum(table(ds_redline$model)), 
         table(ds_redline$model)[[1]])

# Draw curve manually
box1 <- grid.get("human_or_mouse")
box2 <- grid.get("suppf_pres")

grid.curve(grobX(box1, 0) + unit(-2, "mm"),
           grobY(box1, 270) + unit(length(box1$cnames) - match(table(ds_redline$model)[[1]], box1$cnames) + 0.5, "lines"),
           grobX(box2, 0),
           grobY(box2, 270) + unit(length(box2$cnames) - match(nrow(suppfilenames_present), box2$cnames) + 0.5, "lines"),
           curvature = -0.01,
           inflect = TRUE,
           arrow = arrow(angle = 15,
                         length = unit(2, "mm"),
                         type = "closed"), gp = gpar(fill = "black"))

grid.circle(grobX(box1, 0) + unit(-2, "mm"),
            grobY(box1, 270) + unit(length(box1$cnames) - match(table(ds_redline$model)[[1]], box1$cnames) + 0.5, "lines"),
            r = 0.5 * unit(1, "mm"), 
            gp = gpar(fill = "white"))

boxCurve("suppf_pres", "suppf_oi",
         nrow(suppfilenames_present_unnested), 
         nrow(suppfiles_of_interest))

boxCurve("suppf_oi", "pvals",
         nrow(suppfiles_of_interest),
         length(unlist(pvalue_sets$pvalues, recursive = FALSE)),
         c("right", "down"))


# Add text ----------------------------------------------------------------


