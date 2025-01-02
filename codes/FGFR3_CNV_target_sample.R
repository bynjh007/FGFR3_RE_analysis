# loading hg19 genomic information
hg19_info = R.utils::loadToEnv("~/Dropbox/Work_Netherlands/FGFR3/R/HMF_data_purple_FGFR3.RData")[["hg19_info"]]

library(dplyr)

samp_info_meta = read.table("~/Dropbox/Work_Netherlands/FGFR3/HMF_metadata.tsv", sep = "\t", stringsAsFactors = F, header = T, comment.char = "")

# FGFR3 genomic range (uc021pzy)
ind = data.frame(hg19_info$trx) %>% filter(grepl("uc003gdu.2", tx_name)) %>% pull(tx_id)
t_exon = hg19_info$trx_exon[[ind]]
t_intron = rev(hg19_info$trx_intron[[ind]])
df1 = data.frame(data.frame(t_exon)[, c(1,2,3,5)], id = paste("Exon", 1:length(t_exon)))
df2 = rbind(data.frame(data.frame(t_intron)[, c(1,2,3,5)], id = paste("Intron", 1:length(t_intron))), NA)
FGFR3_info = (gdata::interleave(df1, df2))
FGFR3_info = GRanges(FGFR3_info[-nrow(FGFR3_info),])

purple_path = "~/Dropbox/Work_Korea/HMF/HMF_KWF/somatics"
# loading the purple annotation data - data was processed in Turing and downloaded to a local folder
FGFR3_SV_PURPLE_ALL = R.utils::loadToEnv("~/Dropbox/Work_Netherlands/FGFR3/R/HMF_data_purple_FGFR3.RData")[["FGFR3_PURPLE_ALL"]] 
SAMPLE_LIST_ALL = R.utils::loadToEnv("~/Dropbox/Work_Netherlands/FGFR3/R/HMF_data_purple_FGFR3.RData")[["target_samp"]]
FGFR3_SV_PURPLE_ALL = FGFR3_SV_PURPLE_ALL %>% mutate(C_trunc = replace(C_trunc, (grepl("Exon", up.loc_in_trx) & C_trunc == T), FALSE))

# samples to be analyzed : samples with non-SGL SVs from purple
FGFR3_SV_PURPLE_TARGET = FGFR3_SV_PURPLE_ALL %>% filter(RE_type != "SGL")

target_samp = unique(FGFR3_SV_PURPLE_TARGET %>% pull(sample_name))
samp_info_target = samp_info_meta[match(target_samp, samp_info_meta$setName),]

# CNV
FGFR3_CNV = rep(NA, length(target_samp))
FGFR3_CNV_brkpt = c()
for(i in 1:length(target_samp)){
  t_path = file.path(purple_path, target_samp[i], "purple")
  t_files = list.files(t_path)
  t_df = read.table(file.path(t_path, t_files[grepl("cnv.gene", t_files)]), header = T, stringsAsFactors = F, sep = "\t")
  FGFR3_CNV[i] = t_df %>% filter(gene == "FGFR3") %>% pull(maxCopyNumber)
  
  # somatic CNV segment
  t_cnv = read.table(file.path(t_path, t_files[grepl("cnv.somatic", t_files)]), header = T, stringsAsFactors = F, sep = "\t")
  t_start = t_cnv[,c(1,2,2)] %>% mutate(chromosome = paste("chr", chromosome, sep = ""))
  colnames(t_start) = colnames(t_start) = c("chromosome", "start", "end")
  t_start = GRanges(t_start)
  t_start_id = rep(NA, length(t_start))
  t_start_id[queryHits(findOverlaps(t_start, FGFR3_info))] = FGFR3_info$id[subjectHits(findOverlaps(t_start, FGFR3_info))]
  
  t_end = t_cnv[,c(1,3,3)] %>% mutate(chromosome = paste("chr", chromosome, sep = ""))
  colnames(t_end) = colnames(t_end) = c("chromosome", "start", "end")
  t_end = GRanges(t_end)
  t_end_id = rep(NA, length(t_end))
  t_end_id[queryHits(findOverlaps(t_end, FGFR3_info))] = FGFR3_info$id[subjectHits(findOverlaps(t_end, FGFR3_info))]
  
  FGFR3_CNV_brkpt = rbind(FGFR3_CNV_brkpt, 
                          t_cnv %>% 
                            mutate(start_id = t_start_id, end_id = t_end_id, sample_id = target_samp[i]) %>% 
                            filter(!is.na(start_id) | !is.na(end_id)))
  
}
names(FGFR3_CNV) = target_samp

# for samples with entire E1-E17 CNV status
t_samp = unique(FGFR3_CNV_brkpt %>% filter(start_id %in% c("Intron 17", "Exon 18") & is.na(end_id)) %>% pull(sample_id))
samp_with_brkpt_I17_E18 = c()
for(i in 1:length(t_samp)){
  E1_E17 = FGFR3_CNV_brkpt %>% filter(sample_id == t_samp[i], end_id %in% c("Intron 17", "Exon 18") & is.na(start_id)) %>% pull(copyNumber)
  E18 = FGFR3_CNV_brkpt %>% filter(sample_id == t_samp[i], start_id %in% c("Intron 17", "Exon 18") & is.na(end_id)) %>% pull(copyNumber)
  if(length(E1_E17) == 1 & length(E18) == 1){
    samp_with_brkpt_I17_E18 = rbind(samp_with_brkpt_I17_E18, data.frame(E1_E17 = E1_E17, E18 = E18, CN_diff = E1_E17-E18, sample_id = t_samp[i]))
  }
}

save(FGFR3_CNV, samp_with_brkpt_I17_E18, file = "~/Dropbox/Work_Netherlands/FGFR3/R/FGFR3_CNV_target_sample.RData")
