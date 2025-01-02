#####################################################################################
# libraries
#####################################################################################
library(dplyr)
library(ggplot2)


#####################################################################################
# data loading 
#####################################################################################
onco_data = R.utils::loadToEnv("~/FGFR3/nvim_R/FM_FGFR3_RE_annotation.RData")[["onco_data"]]

coalt_FGFR3_trunc = c()
for(i in 1:6){
    coalt_FGFR3_trunc = rbind(coalt_FGFR3_trunc,
                              xlsx::read.xlsx2(file = "/DATA/projects/j.bhin/Daniel_FGFR3/FM/FGFR3 Trunc Concurrent Alts.xlsx", 
                                               sheetIndex = i, header = T, stringsAsFactors = F))
}

t_samp = onco_data$deidentifiedSpecimenName 
t_mut_all = unique(names(sort(table(coalt_FGFR3_trunc$gene), decreasing = T)))
t_mat_all = matrix(NA, length(t_samp), length(t_mut_all))
for(i in 1:length(t_samp)){
  for(j in 1:length(t_mut_all)){
    t_df = coalt_FGFR3_trunc %>% filter(deidentifiedSpecimenName == t_samp[i], gene == t_mut_all[j])
    t_cn = t_df %>% filter(variantType == "CN") %>% pull(amplificationDeletion_CN)
    t_sv = t_df %>% filter(variantType == "SV") %>% pull(proteinEffect_SV)
    t_mat_all[i,j] = paste(c(t_cn, t_sv), collapse = ", ")
  }
  print(i)
}
colnames(t_mat_all) = t_mut_all
t_mat_all[t_mat_all ==""] = NA

coalt_pan = data.frame(onco_data, t_mat_all)

#####################################################################################
# mutation re-labeling
#####################################################################################
# Mutation in the other genes
mut_sum_func = function(X){
  
  mut_type = function(Y){
    # splice site
    ss = grepl("splice", Y)
    # truncation (nonsense, frameshift)
    tr = grepl("[*]", Y)
    # non-frameshift mutation by indel
    nonfs = grepl("ins|del", Y) & !grepl("splice", Y) & !grepl("deletion", Y)
    # promoter mut
    pr = grepl("promoter", Y)
    # amplification
    amp = grepl("amplification", Y)
    # deletion
    del = grepl("deletion", Y)
    tmp_type = rowSums(rbind(ss, tr, nonfs, pr, amp, del))
    return(tmp_type)
  }
  
  # multiple variations
  if(grepl("[,]", X)){
    n = stringr::str_count(X, "[,]") + 1
    tmp = stringr::str_split_fixed(X, "[,]", n)
    mut_sum = mut_type(tmp) # mutation type except for missense
    
    # add missense mutation
    mut_sum = c(n-sum(mut_sum), mut_sum)
    mut_sum[mut_sum>=1] = 1
    
    # unique variations  
  } else {
    mut_sum = mut_type(X) # mutation type except for missense
    
    # add missense mutation
    if(sum(mut_sum)==0){
      mut_sum = c(1, mut_sum)
    } else {
      mut_sum = c(0, mut_sum)
    }
  }
  
  names(mut_sum) = c("missense", "splice_site", "truncation", "nonframeshift_indel", "promoter", "amplification", "deletion")
  out_code_sum = paste(names(mut_sum[mut_sum!=0]), collapse = ";")
  
  # take representative mutation in case the samples with multiple mutations (Truncation >=Missense mutation)
  mut_sum_rep = c(mut_sum[1], mut_sum[3], sum(mut_sum[c(2,4)]), mut_sum[5], mut_sum[6:7])
  mut_sum_rep[mut_sum_rep>=1] = 1
  names(mut_sum_rep) = c("missense", "truncation", "splice_nonfs",  "promoter", "amplification", "deletion")
  if(mut_sum_rep[1]==1 & mut_sum_rep[2]==1){ # Truncation > Missense mutation
    mut_sum_rep[1] =0
  }
  if(mut_sum_rep[2]==1 & mut_sum_rep[3]==1){ # Truncation > splice & nonfs
    mut_sum_rep[3] =0
  }
  if(mut_sum_rep[1]==1 & mut_sum_rep[3]==1){ # splice & nonfs > Missense mutation
    mut_sum_rep[1] =0
  }
  out_code_sum_rep = paste(names(mut_sum_rep[mut_sum_rep!=0]), collapse = ";")
  
  out_code = list(sum = out_code_sum_rep, all = out_code_sum)
  return(out_code)
}

# for all tumors
mut_sum_mat_pan = matrix(NA, nrow(t_mat_all), ncol(t_mat_all))
mut_all_mat_pan = matrix(NA, nrow(t_mat_all), ncol(t_mat_all))
for(i in 1:nrow(t_mat_all)){
  for(j in 1:ncol(t_mat_all)){
    if(!is.na(t_mat_all[i,j])){
      t = mut_sum_func(t_mat_all[i,j])
      mut_sum_mat_pan[i,j] = t$sum
      mut_all_mat_pan[i,j] = t$all
    }
  }
}
colnames(mut_sum_mat_pan) = colnames(t_mat_all)
colnames(mut_all_mat_pan) = colnames(t_mat_all)
rownames(mut_sum_mat_pan) = coalt_pan$deidentifiedSpecimenName
rownames(mut_all_mat_pan) = coalt_pan$deidentifiedSpecimenName
mut_binary_mat_pan = (!is.na(mut_sum_mat_pan)) * 1

rm(i, j, t, t_cn, t_df, t_mat_all, t_mut_all, t_samp, t_sv)

save.image("~/FGFR3/nvim_R/FM_FGFR3_coalterations.RData")

