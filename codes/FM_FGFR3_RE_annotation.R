# calling for resquired functions for annotation
source("~/Dropbox/Work_Netherlands/FGFR3/R/Functions_FM_FGFR3_RE_annotation.R")

onco_data = readxl::read_xlsx(path = "~/Dropbox/Work_Netherlands/FGFR3/R/data/FGFR3 Truncs.xlsx", sheet = "ALL variants")
onco_data = onco_data %>%
  mutate(Chr.1_RE = stringr::str_split_fixed(pos1_RE, ":|-", 3)[,1],
         Chr.2_RE = stringr::str_split_fixed(pos2_RE, ":|-", 3)[,1])

#############################################
# annotation based on breakpoint from FM
#############################################
annot_RE = matrix(NA, nrow(onco_data), 2)

# up
a = unique(onco_data$breakpoint_5p_by_Jessica)
a = a[a != "NA"]
t_u = lapply(a, FUN = function(X){
  re_annot(brkpt = X)$RE_type
})
t_u = unlist(t_u)
names(t_u) = a

# down
a = unique(onco_data$breakpoint_3p_by_Jessica)
a = a[a != "NA"]
t_d = lapply(a, FUN = function(X){
  re_annot(brkpt = X)$RE_type
})
t_d = unlist(t_d)
names(t_d) = a

for(i in 1:nrow(onco_data)){
  annot_RE[i,1] = t_u[match(onco_data$breakpoint_5p_by_Jessica[i], names(t_u))]
  annot_RE[i,2] = t_d[match(onco_data$breakpoint_3p_by_Jessica[i], names(t_d))]
}

onco_data = onco_data %>% 
  mutate(breakpoint_5p_by_Jessica_annot = annot_RE[,1],
         breakpoint_3p_by_Jessica_annot = annot_RE[,2],
         FGFR3_RE_loc = ifelse(grepl("FGFR3", breakpoint_5p_by_Jessica), "upstream", "downstream"),
         FGFR3_RE_loc = replace(FGFR3_RE_loc, breakpoint_5p_by_Jessica == "NA" & breakpoint_3p_by_Jessica == "NA", NA),
         FGFR3_brkpt_loc = ifelse(grepl("FGFR3", breakpoint_5p_by_Jessica), breakpoint_5p_by_Jessica, breakpoint_3p_by_Jessica),
         FGFR3_brkpt_loc = replace(FGFR3_brkpt_loc, breakpoint_5p_by_Jessica == "NA" & breakpoint_3p_by_Jessica == "NA", NA),
         RE_with_TACC3 = grepl("TACC3", breakpoint_3p_by_Jessica),
         TACC3_brkpt_loc = ifelse(grepl("TACC3", breakpoint_3p_by_Jessica), breakpoint_3p_by_Jessica, NA))

onco_data$FGFR3_brkpt_loc = gsub("FGFR3 ", "", onco_data$FGFR3_brkpt_loc)
onco_data$TACC3_brkpt_loc = gsub("TACC3 ", "", onco_data$TACC3_brkpt_loc)
onco_data = onco_data %>%
  mutate(TACC3_brkpt_loc_sum = ifelse(TACC3_brkpt_loc %in% names(which(table(TACC3_brkpt_loc)>=50)), TACC3_brkpt_loc, "others"),
         TACC3_brkpt_loc_sum = replace(TACC3_brkpt_loc_sum, TACC3_brkpt_loc_sum == "others" & RE_with_TACC3 == FALSE, NA))

save.image("~/Dropbox/Work_Netherlands/FGFR3/R/FM_FGFR3_RE_annotation.RData")
