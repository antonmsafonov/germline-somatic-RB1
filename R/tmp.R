












# Generate oncoprints
```{r generate_oncoprints, include = FALSE, echo = FALSE}
{
get_type_fun = function(x) strsplit(x, ";")[[1]]

  col_op = c("Amplification" = "red", "Deletion" = "blue", "Inframe_driver" = "brown4", "Missense_driver" =  "darkgreen", "Splice_driver" = "darkorange1", "Truncating_driver" = "black", "Fusion" = "mediumorchid1", "No_Mut" = "#CCCCCC")

alter_fun = list(
 background = function(x, y, w, h) {
   grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"),
      gp = gpar(fill = "#CCCCCC", col = NA))
 }, 
  No_Mut = function(x, y, w, h) {
   grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"),
      gp = gpar(fill = col_op["No_Mut"], col = NA))
 }, 
 Amplification = function(x, y, w, h) {
   grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"),
      gp = gpar(fill = col_op["Amplification"], col = NA))
 },
  Deletion = function(x, y, w, h) {
   grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"),
      gp = gpar(fill = col_op["Deletion"], col = NA))
 },
  Inframe_driver = function(x, y, w, h) {
   grid.rect(x, y, w-unit(2, "pt"), h*0.33,
      gp = gpar(fill = col_op["Inframe_driver"], col = NA))
 }, 
  Missense_driver = function(x, y, w, h) {
   grid.rect(x, y, w-unit(2, "pt"), h*0.33,
      gp = gpar(fill = col_op["Missense_driver"], col = NA))
 }, 
  Splice_driver = function(x, y, w, h) {
   grid.rect(x, y, w-unit(2, "pt"), h*0.33,
      gp = gpar(fill = col_op["Splice_driver"], col = NA))
 }, 
  Truncating_driver = function(x, y, w, h) {
   grid.rect(x, y, w-unit(2, "pt"), h*0.33,
      gp = gpar(fill = col_op["Truncating_driver"], col = NA))
 }, 
  Fusion = function(x, y, w, h) {
   grid.rect(x, y, w-unit(2, "pt"), h*0.33,
      gp = gpar(fill = col_op["Fusion"], col = NA))
  }
)
}

heatmap_legend_param = list(title = "Variants (Annotated Drivers)", at = c("Amplification", "Deletion", "Inframe_driver", "Missense_driver", "Splice_driver", "Truncating_driver", "Fusion"), labels = c("Amplification", "Deep Deletion", "Inframe mutation (putative driver)", "Missense mutation (putative driver)", "Splice Site mutation (putative driver)", "Truncating mutation (putative driver)", "Fusion (putative driver)")) 


oncoprint_levels = c("No_Mut", "Splice_driver", "Missense_driver", "Inframe_driver", "Truncating_driver", "Amplification", "Deletion", "Fusion")
```

```{r generate_brca1_oncoprint, include = FALSE, echo = FALSE}
brca1.op.df = read.csv("../data/Supp_Fig1b_BRCA1_ANON.csv")
brca1.op.df = brca1.op.df %>% mutate(across(everything(), ~ ifelse(. == "No_Mut", "", .)))

brca1.op.df = brca1.op.df %>% select(-X) %>% rename(Patient.ID = Anonymized.ID)

brca1.op.df = brca1.op.df %>% mutate(across(-c(Patient.ID, Receptor.Status, LOH, Sample.Type), ~ factor(.x, levels = oncoprint_levels)))

brca1.op.df = brca1.op.df  %>% arrange(-as.numeric(TP53),
-as.numeric(MYC),
-as.numeric(NF1),
-as.numeric(GATA3),
-as.numeric(CCND1),
-as.numeric(ERBB2),
-as.numeric(ARID1A),
-as.numeric(CDH1),
-as.numeric(MAP3K1),
-as.numeric(PIK3CA),
-as.numeric(MAP2K4))

ann.op.b1 = brca1.op.df %>% select(Patient.ID, Sample.Type, Receptor.Status, LOH)
b1.mat = brca1.op.df %>% select(-Sample.Type, -Receptor.Status, -LOH)

b1.samp.order = ann.op.b1$Patient.ID
b1m = b1.mat %>%
    column_to_rownames(var = "Patient.ID") %>%
    as.matrix()
brca1.op.matrix = t(b1m)

ann.op.b1 = ann.op.b1 %>% select(Sample.Type, Receptor.Status, LOH)
col.ann = list('Sample.Type' = c("Primary" = "#95aeda", "Metastasis" = "#406eb5"), 
                "Receptor.Status" = c("HR+/HER2-" = "#bd202e", "TNBC" = "#fed6e9", "HER2+" = "#f8add2"), 
                "LOH" = c("TRUE" = "#00ac5f", "FALSE" = "#7bc798"))

annot.b1 = HeatmapAnnotation(df = ann.op.b1, which = "col", col = col.ann, gap = unit(1, "mm"))

column_title_b1 = ""
ht_opt$COLUMN_ANNO_PADDING = unit(8, "mm")
op.b1 = oncoPrint(as.matrix(brca1.op.matrix), alter_fun = alter_fun, col = col_op, column_order = b1.samp.order,pct_side = "right", row_names_side = "left", top_annotation = annot.b1, column_title = column_title_b1, heatmap_legend_param = heatmap_legend_param)

```

```{r generate_brca2_oncoprint, include = FALSE, echo = FALSE}
brca2.op.df = read.csv("../data/Supp_Fig1c_BRCA2_ANON.csv")
brca2.op.df = brca2.op.df %>% mutate(across(everything(), ~ ifelse(. == "No_Mut", "", .)))
brca2.op.df = brca2.op.df %>% select(-X) %>% rename(Patient.ID = Anonymized.ID)

brca2.op.df = brca2.op.df %>% mutate(across(everything(), ~ ifelse(. == "No_Mut", "", .)))


brca2.op.df = brca2.op.df %>% mutate(across(-c(Patient.ID, Receptor.Status, LOH, Sample.Type), ~ factor(.x, levels = oncoprint_levels)))

brca2.op.df = brca2.op.df %>% arrange(-as.numeric(RB1),
-as.numeric(MYC),
-as.numeric(AURKA),
-as.numeric(PIK3CA),
-as.numeric(PTEN),
-as.numeric(GATA3),
-as.numeric(CDH1),
-as.numeric(ARID1A),
-as.numeric(ESR1),
-as.numeric(ERBB2),
-as.numeric(EGFR),
-as.numeric(NF1))
genes.op.brca2 = c("RB1", "MYC", "AURKA", "PIK3CA", "PTEN", "GATA3", "CDH1", "ARID1A", "ESR1", "ERBB2", "EGFR", "NF1") 
ann.op.b2 = brca2.op.df %>% select(Patient.ID, Sample.Type, Receptor.Status, LOH)
b2.mat = brca2.op.df %>% select(-Sample.Type, -Receptor.Status, -LOH)

b2.samp.order = ann.op.b2$Patient.ID
b2m = b2.mat %>%
    column_to_rownames(var = "Patient.ID") %>%
    as.matrix()

brca2.op.matrix = t(b2m)

ann.op.b2 = ann.op.b2 %>% select(Sample.Type, Receptor.Status, LOH)
col.ann = list('Sample.Type' = c("Primary" = "#95aeda", "Metastasis" = "#406eb5"), 
                "Receptor.Status" = c("HR+/HER2-" = "#bd202e", "TNBC" = "#fed6e9", "HER2+" = "#f8add2"), 
                "LOH" = c("TRUE" = "#00ac5f", "FALSE" = "#7bc798"))

annot.b2 = HeatmapAnnotation(df = ann.op.b2, which = "col", col = col.ann, gap = unit(1, "mm"))

column_title_b2 = ""
ht_opt$COLUMN_ANNO_PADDING = unit(8, "mm")
op.b2 = oncoPrint(
  as.matrix(brca2.op.matrix), 
  alter_fun = alter_fun, 
  col = col_op, 
  column_order = b2.samp.order,
  row_order = genes.op.brca2,
  pct_side = "right", 
  row_names_side = "left", 
  top_annotation = annot.b2, 
  column_title = column_title_b2, 
  heatmap_legend_param = heatmap_legend_param)

```

```{r generate_chek2_oncoprint, include = FALSE, echo = FALSE}
chek2.op.df = read.csv("../data/Supp_Fig1d_CHEK2_ANON.csv")
chek2.op.df = chek2.op.df %>% mutate(across(everything(), ~ ifelse(. == "No_Mut", "", .)))

chek2.op.df = chek2.op.df %>% select(-X) %>% rename(Patient.ID = Anonymized.ID)

chek2.op.df = chek2.op.df %>% mutate(across(-c(Patient.ID, Receptor.Status, LOH, Sample.Type), ~ factor(.x, levels = oncoprint_levels)))

chek2.op.df = chek2.op.df  %>% arrange(-as.numeric(PIK3CA),
-as.numeric(CCND1),
-as.numeric(GATA3),
-as.numeric(TP53),
-as.numeric(CDH1),
-as.numeric(AURKA),
-as.numeric(MDM2),
-as.numeric(RUNX1),
-as.numeric(ARID1A),
-as.numeric(MAP2K4),
-as.numeric(MAP3K1))

ann.op.ck2 = chek2.op.df %>% select(Patient.ID, Sample.Type, Receptor.Status, LOH)
c2.mat = chek2.op.df %>% select(-Sample.Type, -Receptor.Status, -LOH)

ck2.samp.order = ann.op.ck2$Patient.ID
c2m = c2.mat %>%
    column_to_rownames(var = "Patient.ID") %>%
    as.matrix()

chek2.op.matrix = t(c2m)

ann.op.ck2 = ann.op.ck2 %>% select(Sample.Type, Receptor.Status, LOH)
col.ann = list('Sample.Type' = c("Primary" = "#95aeda", "Metastasis" = "#406eb5"), 
                "Receptor.Status" = c("HR+/HER2-" = "#bd202e", "TNBC" = "#fed6e9", "HER2+" = "#f8add2"), 
                "LOH" = c("TRUE" = "#00ac5f", "FALSE" = "#7bc798"))

annot.ck2 = HeatmapAnnotation(df = ann.op.ck2, which = "col", col = col.ann, gap = unit(1, "mm"))

column_title_ck2 = ""
ht_opt$COLUMN_ANNO_PADDING = unit(8, "mm")
op.ck2 = oncoPrint(
  as.matrix(chek2.op.matrix), 
  alter_fun = alter_fun, 
  col = col_op, 
  column_order = ck2.samp.order,
  pct_side = "right", 
  row_names_side = "left", 
  top_annotation = annot.ck2, 
  column_title = column_title_ck2, 
  heatmap_legend_param = heatmap_legend_param)

```

```{r generate_atm_oncoprint, include = FALSE, echo = FALSE}
atm.op.df = read.csv("../data/Supp_Fig1e_ATM_ANON.csv")
atm.op.df = atm.op.df %>% mutate(across(everything(), ~ ifelse(. == "No_Mut", "", .)))
atm.op.df = atm.op.df %>% select(-X) %>% rename(Patient.ID = Anonymized.ID)

atm.op.df = atm.op.df %>% mutate(across(-c(Patient.ID, Receptor.Status, LOH, Sample.Type), ~ factor(.x, levels = oncoprint_levels)))

atm.op.df = atm.op.df%>% arrange(-as.numeric(PIK3CA),
-as.numeric(GATA3),
-as.numeric(MYC),
-as.numeric(CCND1),
-as.numeric(AGO2),
-as.numeric(CDH1),
-as.numeric(TP53),
-as.numeric(MDM2),
-as.numeric(ARID1A),
-as.numeric(MAP3K1),
-as.numeric(MAP2K4))

ann.op.atm = atm.op.df %>% select(Patient.ID, Sample.Type, Receptor.Status, LOH)
atm.mat = atm.op.df %>% select(-Sample.Type, -Receptor.Status, -LOH)

atm.samp.order = ann.op.atm$Patient.ID

atm.m =atm.mat %>%
    column_to_rownames(var = "Patient.ID") %>%
    as.matrix()

atm.op.matrix = t(atm.m)

ann.op.atm = ann.op.atm %>% select(Sample.Type, Receptor.Status, LOH)
col.ann = list('Sample.Type' = c("Primary" = "#95aeda", "Metastasis" = "#406eb5"), 
                "Receptor.Status" = c("HR+/HER2-" = "#bd202e", "TNBC" = "#fed6e9", "HER2+" = "#f8add2"), 
                "LOH" = c("TRUE" = "#00ac5f", "FALSE" = "#7bc798"))

annot.atm = HeatmapAnnotation(df = ann.op.atm, which = "col", col = col.ann, gap = unit(1, "mm"))

column_title_atm = ""
ht_opt$COLUMN_ANNO_PADDING = unit(8, "mm")
op.atm = oncoPrint(as.matrix(atm.op.matrix), alter_fun = alter_fun, col = col_op, column_order = atm.samp.order, pct_side = "right", row_names_side = "left", top_annotation = annot.atm, column_title = column_title_atm, heatmap_legend_param = heatmap_legend_param)
```

```{r generate_palb2_oncoprint, include = FALSE, echo = FALSE}
palb2.op.df = read.csv("../data/Supp_Fig1f_PALB2_ANON.csv")
palb2.op.df = palb2.op.df %>% mutate(across(everything(), ~ ifelse(. == "No_Mut", "", .)))
palb2.op.df = palb2.op.df %>% select(-X) %>% rename(Patient.ID = Anonymized.ID)

palb2.op.df = palb2.op.df %>% mutate(across(-c(Patient.ID, Receptor.Status, LOH, Sample.Type), ~ factor(.x, levels = oncoprint_levels)))

palb2.op.df = palb2.op.df %>% arrange(-as.numeric(TP53),
-as.numeric(PIK3CA),
-as.numeric(MYC),
-as.numeric(CCND1),
-as.numeric(MAP3K1),
-as.numeric(NF1),
-as.numeric(ERBB2),
-as.numeric(GATA3),
-as.numeric(CDH1),
-as.numeric(MAP2K4))

ann.op.palb2 = palb2.op.df %>% select(Patient.ID, Sample.Type, Receptor.Status, LOH)
p2.mat = palb2.op.df %>% select(-Sample.Type, -Receptor.Status, -LOH)

p2.samp.order = ann.op.palb2$Patient.ID

p2.m = p2.mat %>%
    column_to_rownames(var = "Patient.ID") %>%
    as.matrix()

palb2.op.matrix = t(p2.m)

ann.op.palb2 = ann.op.palb2 %>% select(Sample.Type, Receptor.Status, LOH)
col.ann = list('Sample.Type' = c("Primary" = "#95aeda", "Metastasis" = "#406eb5"), 
                "Receptor.Status" = c("HR+/HER2-" = "#bd202e", "TNBC" = "#fed6e9", "HER2+" = "#f8add2"), 
                "LOH" = c("TRUE" = "#00ac5f", "FALSE" = "#7bc798"))

annot.palb2 = HeatmapAnnotation(df = ann.op.palb2, which = "col", col = col.ann, gap = unit(1, "mm"))

column_title_palb2 = ""
ht_opt$COLUMN_ANNO_PADDING = unit(8, "mm")
op.palb2 = oncoPrint(as.matrix(palb2.op.matrix), alter_fun = alter_fun, col = col_op, column_order = p2.samp.order, pct_side = "right", row_names_side = "left", top_annotation = annot.palb2, column_title = column_title_palb2, heatmap_legend_param = heatmap_legend_param)
```

