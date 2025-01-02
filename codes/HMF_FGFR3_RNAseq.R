################################################
# genome information
################################################

# loading the same genomic information used for STAR
txdb = AnnotationDbi::loadDb(file = '~/Dropbox/Work_Netherlands/FGFR3/R/data/txdb.GRCh38_genecode_v32_CTAT_lib.sqlite')
hg38_v32 = read.table("~/Dropbox/Work_Netherlands/FGFR3/R/data/GRCh38_genecode_v32_trx.txt", header = T, sep = "\t", stringsAsFactors = F)
hg38_v32_gene = unique(read.table("~/Dropbox/Work_Netherlands/FGFR3/R/data/GRCh38_genecode_v32_trx.txt", header = T, sep = "\t", stringsAsFactors = F)[, c("ENSG_ID", "GENE_SYMBOL")])

# list of transcript in the DB
trx = transcripts(txdb)
gene_exon = exonsBy(txdb, by = "gene")

# intron and exon information
trx_intron = intronsByTranscript(txdb)
trx_exon = exonsBy(txdb)

# TXID is same with index of the "trx", "trx_gene"
trx_gene = AnnotationDbi::select(txdb, key = as.character(trx$tx_id), keytype = "TXID", columns = c("TXNAME", "GENEID"))
exon_info = AnnotationDbi::select(txdb, key = as.character(trx$tx_id), keytype = "TXID", columns = c("EXONID", "EXONNAME", "TXID",  "TXNAME", "GENEID"))

# adding SYMBOL and STRAND information
trx_gene$SYMBOL = hg38_v32$GENE_SYMBOL[match(trx_gene$TXNAME, hg38_v32$ENST_ID)]
trx_gene$STRAND = as.character(strand(trx))

# protein coding region in each transcript (GRange object)
G_cds = cdsBy(txdb, by="tx")

# list of non- scaffold and mitochondrial chromosome
ind_canonical = which(!(grepl("_|M", as.character(seqnames(trx)))))
G_cds = G_cds[which(as.numeric(names(G_cds)) %in% ind_canonical)]

# DNA and AA sequence
cds_seqs = extractTranscriptSeqs(BSgenome.Hsapiens.UCSC.hg38::Hsapiens, G_cds)
cds_seqs_aa = suppressWarnings(Biostrings::translate(cds_seqs))

# filtering out the non-canonical transcript (stop codon in the middle of protein sequence which is caused by non-canonical translation start site)
ind_canonical = which(!(grepl("\\*", substr(as.character(cds_seqs_aa), 1, nchar(as.character(cds_seqs_aa))-1))))
G_cds = G_cds[ind_canonical]
cds_seqs = cds_seqs[ind_canonical]
cds_seqs_aa = cds_seqs_aa[ind_canonical]

# filtering out the non-canonical transcript (stop codon in the middle of protein sequence which is caused by non-canonical translation start site)
ind_canonical = which(!(grepl("\\*", substr(as.character(cds_seqs_aa), 1, nchar(as.character(cds_seqs_aa))-1))))
G_cds = G_cds[ind_canonical]
cds_seqs = cds_seqs[ind_canonical]
cds_seqs_aa = cds_seqs_aa[ind_canonical]

library(dplyr)

# adding the length of the AA and DNA of each transcript and sorting based on 1) AA and 2) DNA length
trx_gene_AA = trx_gene %>% mutate(aa_length = width(cds_seqs_aa)[match(TXID, names(cds_seqs_aa))],
                                  trx_length = width(trx)) %>% 
  arrange(desc(aa_length), desc(trx_length))

# selecting the major transcript for each gene based on the length of AA (longest transcript)
# trx with non-canonical protein sequence is also annotated as NA
total_gene_list = unique(trx_gene_AA$SYMBOL)
trx_gene_AA_major = trx_gene_AA[match(total_gene_list, trx_gene_AA$SYMBOL),] %>% arrange(TXID)
trx_gene_AA_major = trx_gene_AA_major %>% arrange(TXID)

# GRange object for the major transcript
trx_major = trx[trx_gene_AA_major$TXID]

# re-ordering trx table based on the TXID
trx_gene_AA = trx_gene_AA %>% arrange(TXID)

trx_info = list(trx = trx, gene_exon = gene_exon, 
                exon_info = exon_info, trx_intron = trx_intron, trx_exon = trx_exon, 
                trx_major = trx_major, trx_gene_AA = trx_gene_AA, 
                trx_gene_AA_major = trx_gene_AA_major, 
                G_cds = G_cds, cds_seqs = cds_seqs, cds_seqs_aa = cds_seqs_aa)


################################################
# libraries
################################################
# before the REs are annotated, the output of STAR should be refined by the reads mapping to the target genes
# cd /DATA/projects/j.bhin/Daniel_FGFR3/HMF/
# ./preprocessing_target.sh -p RNAseq/star_fusion/ -s samples.HMF.csv -t FGFR3 -g /DATA/projects/j.bhin/reference/GRCh38_gencode_v32_CTAT_lib_Dec062019.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf


#renv::install("bynjh007/fusion.detection")
library(fusion.detection)

################################
# file loading
################################
files = list.files("~/Dropbox/Work_Netherlands/FGFR3/HMF/star_fusion/")

FGFR3_RE = c()
for(i in 1:length(files)){
  chimeric_before_filter = paste("~/Dropbox/Work_Netherlands/FGFR3/HMF/star_fusion/",
                                 files[i],
                                 "/star-fusion.junction_breakpts_to_genes.txt.FGFR3",
                                 sep = "")
  chimeric_after_filter = paste("~/Dropbox/Work_Netherlands/FGFR3/HMF/star_fusion/",
                                files[i],
                                "/star-fusion.junction_breakpts_to_genes.txt.FGFR3.filtered",
                                sep = "")
  star_log = read.table(paste("~/Dropbox/Work_Netherlands/FGFR3/HMF/Log/star_align/", files[i], ".Log.final.out", sep = ""),
                        sep = "\t", stringsAsFactors = F, fill = T)
  total_number_reads = as.numeric(star_log$V2[5])
  t_chimeric_annot = breaks_to_gene_table(chimeric_before_filter, chimeric_after_filter, "FGFR3")
  if(!is.null(t_chimeric_annot)){
    t_summary_brkpt = summarizing_fusions(t_chimeric_annot, total_number_reads)
    t_summary_brkpt = t_summary_brkpt %>% filter(count>=3) %>% mutate(sample_id = files[i])
    FGFR3_RE = rbind(FGFR3_RE, t_summary_brkpt)
  }
  print(i)
}

FGFR3_RE = FGFR3_RE %>% filter(!(up.symbol == "FGFR3" & down.symbol == "FGFR3"))



################################
# STAR-fusion output
################################
# paired-end data
FGFR3_star_fusion = c()
for(i in 1:length(files)){
  t_fus = read.table(paste("~/Dropbox/Work_Netherlands/FGFR3/HMF/star_fusion/", files[i], "/star-fusion.fusion_predictions.abridged.coding_effect.tsv", sep = ""),
                     header = T, sep = "\t", stringsAsFactors = F, comment.char = "")
  ind = grep("FGFR3", t_fus$X.FusionName)
  if(length(ind)>0){
    FGFR3_star_fusion = rbind(FGFR3_star_fusion, data.frame(t_fus[ind,], sample_id = files[i], stringsAsFactors = F))
  }
  print(i)
}

# annotation of the breakpoint
for(i in 1:nrow(FGFR3_star_fusion)){
  # transcripts for upstream and downstream fusion partner
  g_l = strsplit(unlist(strsplit(FGFR3_star_fusion$LeftGene[i], "^", 2))[2], "[.]")[[1]][1]
  g_r = strsplit(unlist(strsplit(FGFR3_star_fusion$RightGene[i], "^", 2))[2], "[.]")[[1]][1]
  
  # brkpt
  t_brkpt_l = unlist(strsplit(FGFR3_star_fusion$LeftBreakpoint[i], ":", 3))
  #t_pos = ifelse(t_brkpt_l[3] == "-", as.numeric(t_brkpt_l[2])-1, as.numeric(t_brkpt_l[2])+1)
  t_brkpt_l = GRanges(seqnames = t_brkpt_l[1], ranges = IRanges(as.numeric(t_brkpt_l[2]), as.numeric(t_brkpt_l[2])), strand = t_brkpt_l[3])
  
  t_brkpt_r = unlist(strsplit(FGFR3_star_fusion$RightBreakpoint[i], ":", 3))
  #t_pos = ifelse(t_brkpt_r[3] == "-", as.numeric(t_brkpt_r[2])+1, as.numeric(t_brkpt_r[2])-1)
  t_brkpt_r = GRanges(seqnames = t_brkpt_r[1], ranges = IRanges(as.numeric(t_brkpt_r[2]), as.numeric(t_brkpt_r[2])), strand = t_brkpt_r[3])
  
  t_range_l = data.frame(trx_info$gene_exon[which(stringr::str_split_fixed(names(trx_info$gene_exon), "[.]", 2)[,1] == g_l)])
  t_range_r = data.frame(trx_info$gene_exon[which(stringr::str_split_fixed(names(trx_info$gene_exon), "[.]", 2)[,1] == g_r)])
  
  ########################################
  # left
  ########################################
  # overlapping gene is in "-" strand --> check whether junction is at the end of known exon
  if(t_range_l$strand[1] == "-"){
    t_exon_l = t_range_l$exon_name[which((t_range_l$seqnames %in% seqnames(t_brkpt_l)) & (t_range_l$start %in% start(t_brkpt_l)))]
    # overlapping gene is in "+" strand --> check whether junction is at the end of known exon
  } else {
    t_exon_l = t_range_l$exon_name[which((t_range_l$seqnames %in% seqnames(t_brkpt_l)) & (t_range_l$end %in% start(t_brkpt_l)))]
  }
  
  # reference transcript to annotate breakpoint
  t_trx_l = trx_info$trx_gene_AA %>% filter(TXNAME %in% trx_info$exon_info$TXNAME[which(exon_info$EXONNAME %in% t_exon_l)]) %>%
    arrange(desc(aa_length), desc(trx_length))
  # junctions from the end of known exon - select the most representative transcript (aa length, trx length)
  # this is for the case whether the junctions from the exon of non-major transcript
  # identify the location of the breakpoint based on the transcript
  if(nrow(t_trx_l)>0){
    t_trx_id_l = t_trx_l$TXID[1]
    splice_type_l = "ref:in-strand"
    # if breakpoint is somewhere in the exon or intron, major transcript will be used as reference
  } else{
    # sorting the transcripts of target gene based on protein and trx length
    t_trx_id_l = (trx_info$trx_gene_AA[queryHits(findOverlaps(trx_info$trx, t_brkpt_l)),] %>%
                    arrange(desc(aa_length), desc(trx_length)) %>% pull(TXID))[1]
    
    # in case of out-of-strand, no trx is mapped
    if(is.na(t_trx_id_l)){
      splice_type_l = "non_ref:out-of-strand"
      t_trx_id_l = (trx_info$trx_gene_AA[queryHits(findOverlaps(trx_info$trx, t_brkpt_l, ignore.strand = T)),] %>%
                      arrange(desc(aa_length), desc(trx_length)) %>% pull(TXID))[1]
    } else {
      splice_type_l = "non_ref:in-strand"
    }
  }
  G_ref_l = create_Gref(genome_info = trx_info, TXID = t_trx_id_l)
  loc_in_trx_l = G_ref_l[queryHits(findOverlaps(G_ref_l, t_brkpt_l, ignore.strand = T))]
  
  # if fusion of the downstream occurs at 5' of the gene, brkpt is not overlapped with reference (brkpt + 1/-1)
  # target range is labeled as just first exon of the transcript as we don't know the exact breakpoint
  if(length(loc_in_trx_l)==0){
    loc_in_trx_l = tail(G_ref_l,1)
    if(t_range_l$strand[1] == "+" & (start(t_brkpt_l)-1) == tail(end(G_ref_l),1)){
      loc_in_trx_r$id = "3P"
    } else if (t_range_l$strand[1] == "-" & (start(t_brkpt_l)+1) == tail(start(G_ref_l),1)){
      loc_in_trx_l$id = "3P"
    } else {
      stop("undefined breakpoint")
    }
  }
  
  
  ########################################
  # right
  ########################################
  # overlapping gene is in "-" strand --> check whether junction is at the end of known exon
  if(t_range_r$strand[1] == "-"){
    t_exon_r = t_range_r$exon_name[which((t_range_r$seqnames %in% seqnames(t_brkpt_r)) & ((t_range_r$end) %in% start(t_brkpt_r)))]
    # overlapping gene is in "+" strand --> check whether junction is at the end of known exon
  } else {
    t_exon_r = t_range_r$exon_name[which((t_range_r$seqnames %in% seqnames(t_brkpt_r)) & ((t_range_r$start) %in% start(t_brkpt_r)))]
  }
  
  # reference transcript to annotate breakpoint
  t_trx_r = trx_info$trx_gene_AA %>% filter(TXNAME %in% trx_info$exon_info$TXNAME[which(exon_info$EXONNAME %in% t_exon_r)]) %>%
    arrange(desc(aa_length), desc(trx_length))
  # junctions from the end of known exon - select the most representative transcript (aa length, trx length)
  # this is for the case whether the junctions from the exon of non-major transcript
  # identify the location of the breakpoint based on the transcript
  if(nrow(t_trx_r)>0){
    t_trx_id_r = t_trx_r$TXID[1]
    splice_type_r = "ref:in-strand"
    # if breakpoint is somewhere in the exon or intron, major transcript will be used as reference
  } else{
    # sorting the transcripts of target gene based on protein and trx length
    t_trx_id_r = (trx_info$trx_gene_AA[queryHits(findOverlaps(trx_info$trx, t_brkpt_r)),] %>%
                    arrange(desc(aa_length), desc(trx_length)) %>% pull(TXID))[1]
    
    # in case of out-of-strand, no trx is mapped
    if(is.na(t_trx_id_r)){
      splice_type_r = "non_ref:out-of-strand"
      t_trx_id_r = (trx_info$trx_gene_AA[queryHits(findOverlaps(trx_info$trx, t_brkpt_r, ignore.strand = T)),] %>%
                      arrange(desc(aa_length), desc(trx_length)) %>% pull(TXID))[1]
    } else {
      splice_type_r = "non_ref:in-strand"
    }
  }
  G_ref_r = create_Gref(genome_info = trx_info, TXID = t_trx_id_r)
  loc_in_trx_r = G_ref_r[queryHits(findOverlaps(G_ref_r, t_brkpt_r, ignore.strand = T))]
  
  # if fusion of the downstream occurs at 5' of the gene, brkpt is not overlapped with reference (brkpt + 1/-1)
  # target range is labeled as just first exon of the transcript as we don't know the exact breakpoint
  if(length(loc_in_trx_r)==0){
    loc_in_trx_r = G_ref_r[1]
    if(t_range_r$strand[1] == "+" & (start(t_brkpt_r)+1) == head(start(G_ref_r),1)){
      loc_in_trx_r$id = "5P"
    } else if (t_range_r$strand[1] == "-" & (start(t_brkpt_l)-1) == head(end(G_ref_r),1)){
      loc_in_trx_r$id = "5P"
    } else {
      stop("undefined breakpoint")
    }
  }
  
  loc_in_trx_l$symbol = trx_info$trx_gene_AA_major %>% filter(stringr::str_split_fixed(GENEID, "[.]", 2)[,1] %in% g_l) %>% pull(SYMBOL)
  loc_in_trx_r$symbol = trx_info$trx_gene_AA_major %>% filter(stringr::str_split_fixed(GENEID, "[.]", 2)[,1] %in% g_r) %>% pull(SYMBOL)
  
  loc_in_trx_l$txname = trx_info$trx_gene_AA %>% filter(TXID == t_trx_id_l) %>% pull(TXNAME)
  loc_in_trx_r$txname = trx_info$trx_gene_AA %>% filter(TXID == t_trx_id_r) %>% pull(TXNAME)
  
  loc_in_trx_l = data.frame(loc_in_trx_l) %>% mutate(label = paste(symbol, txname, id, splice_type_l, sep = ":"))
  loc_in_trx_r = data.frame(loc_in_trx_r) %>% mutate(label = paste(symbol, txname, id, splice_type_r, sep = ":"))
  
  FGFR3_star_fusion$left_annot[i] = loc_in_trx_l$label
  FGFR3_star_fusion$right_annot[i] = loc_in_trx_r$label
  print(i)
}

save(FGFR3_star_fusion, FGFR3_RE, file = "~/Dropbox/Work_Netherlands/FGFR3/R/HMF_FGFR3_RNAseq.RData")
