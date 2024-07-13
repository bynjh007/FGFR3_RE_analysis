library(GenomicFeatures)
library(dplyr)
options(stringsAsFactors = F)

hg_id_19_ENST = read.table("~/Dropbox/Work_Netherlands/FGFR3/R/data/UCSC.hg.19_ID_to_SYMBOL.txt", header = T, stringsAsFactors = F, sep = "\t", comment.char = "")
hg_id_19 = hg_id_19_ENST[,c(1,4)]

##############################################################################
# function for generating required information from genome annotation database
##############################################################################

generating_genome_info = function(txdb, id_sym_tab, genome){
  # id_sym_tab: col1-TXNAME, col2-SYMBOL
  
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
  trx_gene$SYMBOL = id_sym_tab[match(trx_gene$TXNAME, id_sym_tab[,1]), 2]
  trx_gene$STRAND = as.character(strand(trx))
  
  # protein coding region in each transcript (GRange object)
  G_cds = cdsBy(txdb, by="tx")
  
  # list of non- scaffold and mitochondrial chromosome
  ind_canonical = which(!(grepl("_|M", as.character(seqnames(trx)))))
  G_cds = G_cds[which(as.numeric(names(G_cds)) %in% ind_canonical)]
  
  # DNA and AA sequence
  cds_seqs = extractTranscriptSeqs(genome, G_cds)
  cds_seqs_aa = suppressWarnings(Biostrings::translate(cds_seqs))
  
  # filtering out the non-canonical transcript (stop codon in the middle of protein sequence which is caused by non-canonical translation start site)
  ind_canonical = which(!(grepl("\\*", substr(as.character(cds_seqs_aa), 1, nchar(as.character(cds_seqs_aa))-1))))
  G_cds = G_cds[ind_canonical]
  cds_seqs = cds_seqs[ind_canonical]
  cds_seqs_aa = cds_seqs_aa[ind_canonical]
  
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
  
  return(trx_info)
  
}


##############################################################
# generation of information required for breakpoint annotation
##############################################################
hg19_info = generating_genome_info(txdb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene, 
                                   id_sym_tab = hg_id_19, genome = BSgenome.Hsapiens.UCSC.hg19::Hsapiens)

##############################################################
# function
##############################################################
brk_annot_framework = function(G_up, up_orient, G_down, down_orient, genome_info, BSgenome){
  ################################################## 
  # find the overlapping transcript
  ##################################################
  index_ovp = list(up = queryHits(suppressWarnings(findOverlaps(genome_info$trx_major, G_up))),
                   down = queryHits(suppressWarnings(findOverlaps(genome_info$trx_major, G_down))))
  
  G_target = list(up = G_up, down = G_down)
  
  # if there is no overlap for the representative transcript
  # 1 is upstream and 2 is downstream
  target_annot_all = c()
  for(i in 1:2){
    
    if(i == 1){
      target_orient = up_orient
      up_down = "up"
    } else {
      target_orient = down_orient
      up_down = "down"
    }
    
    if(length(index_ovp[[i]]) == 0){
      # find the target region in the all the transcripts
      
      index = queryHits(suppressWarnings(findOverlaps(genome_info$trx, G_target[[i]])))
      if(length(index) !=0){
        index = genome_info$trx_gene_AA[index,] %>% arrange(desc(aa_length), desc(trx_length)) %>% pull(TXID)
        target_annot = brk_annot(ref_trx = genome_info$trx_gene_AA, index = index[1], target_region = G_target[[i]], target_orient = target_orient, 
                                 up_down = up_down, genome_info = genome_info, BSgenome = BSgenome)
        
        # no overlapping and SGL
      } else if(as.character(seqnames(G_target[[i]])) == "chr-1"){
        target_annot = data.frame(G_target[[i]])[,1:2] %>% mutate(gene_id = NA,
                                                                  gene_symbol = NA,
                                                                  trx_name = NA,
                                                                  aa_length = NA,
                                                                  trx_length = NA,
                                                                  loc_in_trx = NA,
                                                                  down_exon = NA,
                                                                  down_exon_frame = "no partner",
                                                                  RE_category = "SGL",
                                                                  orientation = target_orient)
        
      } else {
        target_annot = data.frame(G_target[[i]])[,1:2] %>% mutate(gene_id = NA,
                                                                  gene_symbol = NA,
                                                                  trx_name = NA,
                                                                  aa_length = NA,
                                                                  trx_length = NA,
                                                                  loc_in_trx = NA,
                                                                  down_exon = NA,
                                                                  down_exon_frame = "no partner",
                                                                  RE_category = "intergenic",
                                                                  orientation = target_orient)
      }
      
    } else {
      index = index_ovp[[i]][(genome_info$trx_gene_AA_major[index_ovp[[i]],] %>% mutate(ind = 1:length(index_ovp[[i]])) %>% arrange(desc(aa_length), desc(trx_length)) %>% pull(ind))]
      target_annot = brk_annot(ref_trx = genome_info$trx_gene_AA_major, index = index[1], target_region = G_target[[i]], target_orient = target_orient, up_down = up_down,
                               genome_info = genome_info, BSgenome = BSgenome)
    }
    
    target_annot_all = rbind(target_annot_all, target_annot)
  }
  
  # if both up and downstream are out-of-strand --> this is acually in-strand fusion
  if(sum(target_annot_all$RE_category == "out-of-strand")==2){
    target_annot_all = target_annot_all[c(2,1),]
    target_annot_all$RE_category = "in-strand"
  }
  
  return(target_annot_all)
}


# function to define the breakpoint
brk_annot = function(ref_trx, index, target_region, target_orient, up_down, genome_info, BSgenome){
  
  target_gene_trx = ref_trx[index,]
  target_gene_trx = target_gene_trx %>% mutate(Protein_coding = ifelse(!is.na(aa_length), TRUE, FALSE))
  
  if(nrow(target_gene_trx)>1){
    stop("more than one gene in the brkpt")
  }
  
  # target region is upstream of the RE
  if(up_down == "up"){
    
    if((target_gene_trx$STRAND == "+" & target_orient == 1) | (target_gene_trx$STRAND == "-" & target_orient == (-1))){
      re_categ = "in-strand"
    } else {
      re_categ = "out-of-strand"
    }
    # target region is downstream of the RE
  } else if (up_down == "down") {
    if((target_gene_trx$STRAND == "+" & target_orient == (-1)) | (target_gene_trx$STRAND == "-" & target_orient == 1)){
      re_categ = "in-strand"
    } else {
      re_categ = "out-of-strand"
    }
  } else {
    stop("'up_down' parameter should be either 'up' or 'down'")
  }
  
  
  ################################################## 
  # build the GRange object for the target transcript
  ##################################################
  ind = as.numeric(target_gene_trx$TXID)
  t_exon = genome_info$trx_exon[[ind]]
  t_intron = genome_info$trx_intron[[ind]]
  t_str = unique(as.character(strand(t_exon)))
  
  if(length(t_intron)!=0){
    if(t_str == "-"){
      t_intron = rev(t_intron)
    }
    df1 = data.frame(data.frame(t_exon)[, c(1,2,3,5)], id = paste("Exon", 1:length(t_exon)))
    df2 = rbind(data.frame(data.frame(t_intron)[, c(1,2,3,5)], id = paste("Intron", 1:length(t_intron))), NA)
    df_comb = (gdata::interleave(df1, df2))
    df_comb = df_comb[-nrow(df_comb),]
    G_ref = makeGRangesFromDataFrame(df_comb)
    G_ref$id = df_comb$id
  } else {
    df1 = data.frame(data.frame(t_exon)[, c(1,2,3,5)], id = paste("Exon", 1:length(t_exon)))
    G_ref = makeGRangesFromDataFrame(df1)
    G_ref$id = df1$id
  }
  
  #######################################################################################
  # Find the location of the breakpoint and downstream exons within the target transcript
  #######################################################################################
  # Compare the Ranges between the target transcript and target region (breakpoint)
  temp_ovp = findOverlaps(G_ref, target_region)
  loc_in_trx = G_ref[queryHits(temp_ovp)]
  
  ind = queryHits(temp_ovp)
  
  # following exon from the breakpoint
  if(ind == length(G_ref)){
    down_exon = G_ref[length(G_ref)]
  } else {
    down_exon = G_ref[(ind+1):length(G_ref)]
  }
  down_exon = down_exon[grepl("Exon", down_exon$id)]
  
  #######################################################################################
  # check whether the corresponding trx is protein-coding or not
  #######################################################################################
  if(target_gene_trx$Protein_coding == T){
    # protein coding region in the transcript
    G_cds_target = genome_info$G_cds[[which(names(genome_info$G_cds) == target_gene_trx$TXID)]]
    # Protein sequence for the transcript
    target_cds_seq_AA = as.character(genome_info$cds_seqs_aa[names(genome_info$cds_seqs_aa) == target_gene_trx$TXID][[1]])
    
    #######################################################################################
    # check the following exon is in-frame or out-of-frame or unknown
    #######################################################################################
    # the nearby exon is completely overlapped with exon of full CDS on the transcript (exon in the middle of transcript)
    if(sum(G_cds_target == down_exon[1]) == 1){
      # check whether the protein sequence of the neigbor exon (previous exon for upstream and next exon for downstream) is in-frame or not
      # this can be identified by the aa sequence of the next exon of the breakpoint
      t_seq = as.character(suppressWarnings(Biostrings::translate(Biostrings::getSeq(BSgenome, G_cds_target[which(G_cds_target == down_exon[1])])))[[1]])
      # target exon is inframe --> fully overlapped with protein sequence of the transcript
      if(grepl(t_seq, target_cds_seq_AA, fixed = T)){
        down_exon_frame = "inframe:CDS_full"
      } else {
        down_exon_frame = "out_of_frame:CDS_full"
      }
      
      # the following exon is partly overlapped with CDS of the transcript (exons with translation start (5') or end (3'))
    } else if(length(queryHits(findOverlaps(down_exon[1], G_cds_target))) == 1){
      
      # if the following exon include translation start site -> the following exon should include the first CDS exon
      if(length(queryHits(findOverlaps(down_exon[1], G_cds_target[1])))==1){
        down_exon_frame = "unknown:CDS_start"
        
        # if the following exon include translation end site -> the following exon should include the last CDS exon
      } else if(length(queryHits(findOverlaps(down_exon[1], G_cds_target[length(G_cds_target)])))==1){
        # protein sequence of the last CDS exon
        t_seq = as.character(Biostrings::translate(Biostrings::getSeq(BSgenome, G_cds_target[length(G_cds_target)]))[[1]])
        # check whether the last CDS exon is in-frame or not
        if(grepl(t_seq, target_cds_seq_AA, fixed = T)){
          down_exon_frame = "inframe:CDS_end"
        } else {
          down_exon_frame = "out_of_frame:CDS_end"
        }
        
      } else {
        stop("unexpected REs with partial CDS")
      }
      
      # the following exon is completely UTR
    } else if(length(queryHits(findOverlaps(down_exon[1], G_cds_target))) == 0){
      
      # if there is remaining CDS exons (5'-UTR)
      if(length(queryHits(findOverlaps(down_exon[1], G_cds_target)))!=0){
        down_exon_frame = "unknown:UTR_5"
        # if there is no remaining CDS exons (3'-UTR)
      } else {
        down_exon_frame = "unknown:UTR_3"
      }
      
    } else {
      stop("unexpected RE scenario")
    }
    
    # if target transcript is non-coding  
  } else {
    down_exon_frame = "unknown:noncoding trx"
  }
  down_exon = as.character(down_exon$id[1])
  
  # output
  if(as.character(loc_in_trx$id) != "NA"){
    t_loc = apply(data.frame(loc_in_trx), 1, paste, collapse = ":")
  } else {
    t_loc = "out_of_trx"
  }
  
  target_annot = data.frame(target_region)[,1:2] %>% mutate(gene_id = target_gene_trx$GENEID,
                                                            gene_symbol = target_gene_trx$SYMBOL,
                                                            trx_name = target_gene_trx$TXNAME,
                                                            aa_length = target_gene_trx$aa_length,
                                                            trx_length = target_gene_trx$trx_length,
                                                            loc_in_trx = t_loc,
                                                            down_exon = down_exon,
                                                            down_exon_frame = down_exon_frame,
                                                            RE_category = re_categ,
                                                            orientation = target_orient)
  return(target_annot)
}


# annotation for the vcf file
get_orient = function(purple_vcf){
  if(endsWith(purple_vcf$V5, "[")){
    orient = c(1, -1)
  } else if (startsWith(purple_vcf$V5, "]")){
    orient = c(-1, 1)
  } else if (startsWith(purple_vcf$V5, "[")){
    orient = c(-1, -1)
  } else if (endsWith(purple_vcf$V5, "]")){
    orient = c(1, 1)
    # SGL
  } else {
    if(endsWith(purple_vcf$V5, ".")){
      orient = c(1, 0)
    } else {
      orient = c(-1, 0)
    }
  }
  return(orient)
}


##############################################################
# data processing
##############################################################
purple_path = "/DATA/projects/j.bhin/Daniel_FGFR2/HMF/somatics/"
target_samp = list.files(purple_path)

FGFR3_SV = c()
for(i in 1:length(target_samp)){

  ########################################
  # purple
  ########################################
  t = list.files(paste(purple_path, target_samp[i], sep = "/"))
  purple_SV = read.table(gzfile(paste(purple_path, target_samp[i], t[grepl("purple.sv.ann", t) & !grepl("tbi", t)], sep = "/")), sep = "\t", stringsAsFactors = F, header = F)
  purple_SV_FGFR3 = purple_SV %>% filter(V1 == 4, V2 >=1795039, V2 <=1810599)
  
  if(nrow(purple_SV_FGFR3)!=0){
    t_partner = stringr::str_split_fixed(purple_SV_FGFR3$V5, "\\]|\\[|\\:",4)[,2:3]
    if(nrow(purple_SV_FGFR3)==1){
      t_partner = data.frame(t(t_partner), stringsAsFactors = F)
    } else {
      t_partner = data.frame(t_partner, stringsAsFactors = F)
    }
    t_partner[t_partner[,1] == "", 1] = "-1"
    t_partner[t_partner[,2] == "", 2] = 0
    colnames(t_partner) = c("V1", "V2")

    # annotation
    purple_SV_annot = c()
    for(j in 1:nrow(purple_SV_FGFR3)){
      orient = get_orient(purple_SV_FGFR3[j,])
      G_1 = data.frame(purple_SV_FGFR3[j, c("V1", "V2", "V2")], ".") %>%
        mutate(Orient = orient[1], V1 = paste("chr", V1, sep = ""))
      colnames(G_1) = c("seqnames", "tx_start", "tx_end", "strand", "Orient")
      G_1 = makeGRangesFromDataFrame(G_1)

      G_2 = data.frame(t_partner[j,1:2], t_partner[j,2],  ".") %>%
        mutate(Orient = orient[2], V1 = paste("chr", V1, sep = ""))
      colnames(G_2) = c("seqnames", "tx_start", "tx_end", "strand", "Orient")
      G_2 = makeGRangesFromDataFrame(G_2)

      t_out = suppressWarnings(brk_annot_framework(G_up = G_1, up_orient =  orient[1],
                                                   G_down = G_2, down_orient =  orient[2],
                                                   BSgenome = BSgenome.Hsapiens.UCSC.hg19::Hsapiens,
                                                   genome_info = hg19_info))
      t_out = data.frame(up = t_out[1,], down = t_out[2,])
      purple_SV_annot = rbind(purple_SV_annot, t_out)
    }
    # remove redundant SVs
    purple_SV_annot = purple_SV_annot %>% arrange(up.RE_category)
    t_sv = t(apply(purple_SV_annot %>% select(up.start, down.start), 1, sort))
    t_id = apply(t_sv, 1, paste, collapse = "-")
    purple_SV_annot = purple_SV_annot[match(unique(t_id), t_id),]
    purple_SV_annot = data.frame(purple_SV_annot, sample_name = target_samp[i])
    
    FGFR3_SV = rbind(FGFR3_SV, purple_SV_annot)
  }
  print(i)
}


# function to define in-frame or out-of-frame from purple annotation
RE_type = function(purple_annot){
  # up and down are FGFR3
  if((purple_annot$up.gene_symbol %in% "FGFR3") & (purple_annot$down.gene_symbol %in% "FGFR3")){
    RE = "internal"
    # one of the partner is intergenic
  } else if(is.na(purple_annot$up.gene_symbol) | is.na(purple_annot$down.gene_symbol)){
    if(purple_annot$up.RE_category == "SGL" | purple_annot$down.RE_category == "SGL"){
      RE = "SGL"
    } else {
      RE = "intergenic"
    }
    # in-strand
  } else if((purple_annot$up.RE_category %in% "in-strand") & (purple_annot$down.RE_category %in% "in-strand")){
    # in-frame: both up and downstream is inframe
    if(grepl("inframe", purple_annot$up.down_exon_frame) & grepl("inframe", purple_annot$down.down_exon_frame)){
      if(grepl("Exon", purple_annot$up.loc_in_trx) | grepl("Exon", purple_annot$down.loc_in_trx)){
        RE = "out-of-frame (brkpt is in exon)"
      } else {
        RE = "in-frame"
      }
    } else {
      RE = "out-of-frame"
    }
  } else if((purple_annot$up.RE_category %in% "out-of-strand") | (purple_annot$down.RE_category %in% "out-of-strand")){
    RE = "out-of-strand"
  } else {
    RE = "unknown"
  }
  return(RE)
}

################################################################################
# combined results from PURPLE
################################################################################
t_RE = c()
for(i in 1:nrow(FGFR3_SV)){
  t_RE = c(t_RE, RE_type(FGFR3_SV[i,]))
}


FGFR3_PURPLE_ALL = FGFR3_SV %>%
  mutate(C_trunc = (up.down_exon == "Exon 18" & up.RE_category == "in-strand" & up.gene_symbol == "FGFR3")) %>%
  mutate(Exon_up = ifelse(up.RE_category == "in-strand", 
                          as.numeric(stringr::str_split_fixed(up.down_exon, " ", 2)[,2])-1,
                          as.numeric(stringr::str_split_fixed(up.down_exon, " ", 2)[,2])),
         Exon_down = ifelse(down.RE_category == "in-strand", 
                            as.numeric(stringr::str_split_fixed(down.down_exon, " ", 2)[,2]),
                            as.numeric(stringr::str_split_fixed(down.down_exon, " ", 2)[,2])-1)) %>%
  mutate(RE_type = factor(t_RE, levels = c("in-frame", "out-of-frame (brkpt is in exon)", "out-of-frame", "intergenic", "out-of-strand", "internal", "SGL")), 
         FGFR3_is_upstream = (up.gene_symbol == "FGFR3")) %>%
  arrange(desc(C_trunc), desc(FGFR3_is_upstream), RE_type)


save(target_samp, FGFR3_PURPLE_ALL, hg19_info, file = "~/FGFR3/R_project/HMF_data_purple.RData")
