library(GenomicFeatures)
library(Biostrings)
library(dplyr)
options(stringsAsFactors = F)

hg_id_19_ENST = read.table("~/gene_set/UCSC.hg.19_ID_to_SYMBOL.txt", header = T, stringsAsFactors = F, sep = "\t", comment.char = "")
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
  # find the overlapping representative protein-coding transcript (genome_info$trx_major)
  ##################################################
  index_ovp = list(up = queryHits(suppressWarnings(findOverlaps(genome_info$trx_major, G_up))),
                   down = queryHits(suppressWarnings(findOverlaps(genome_info$trx_major, G_down))))
  
  G_target = list(up = G_up, down = G_down)
  
  re_strand = rep(NA, 2)
  for(i in 1:2){
    target_orient = ifelse(i == 1, up_orient, down_orient)
    
    if(target_orient !=0){ # 
      # no overlapping protein-coding transcript --> mapping to entire transcriptome 
      if(length(index_ovp[[i]])==0){ 
        index = queryHits(suppressWarnings(findOverlaps(genome_info$trx, G_target[[i]])))
        if(length(index) == 0){
          re_strand[i] = "intergenic"
        } else {
          re_strand[i] = brk_orientation(ref_trx = genome_info$trx, index = index, target_orient = target_orient, up_down = ifelse(i==1, "up", "down"))
        }
      } else {
        index = index_ovp[[i]][(genome_info$trx_gene_AA_major[index_ovp[[i]],] %>% mutate(ind = 1:length(index_ovp[[i]])) %>% arrange(desc(aa_length), desc(trx_length)) %>% pull(ind))]
        re_strand[i] = brk_orientation(ref_trx = genome_info$trx_major, index = index, target_orient = target_orient, up_down = ifelse(i==1, "up", "down"))
      }
    } else { # SGL
      re_strand[i] = "SGL"
    }
  }
  
  # if both up and downstream are out-of-strand --> this is actually in-strand fusion
  if(sum(re_strand == "out-of-strand")==2){
    G_target = list(up = G_target[[2]],down = G_target[[1]])
    re_strand = c("in-strand", "in-strand")
    index_ovp = list(up = index_ovp[[2]], down = index_ovp[[1]])
    re_orient = c(down_orient, up_orient)
  } else {
    re_orient = c(up_orient, down_orient)
  }
  
  # if there is no overlap for the representative transcript
  # 1 is upstream and 2 is downstream
  target_annot_all = c()
  for(i in 1:2){
    
    up_down = ifelse(i==1, "up", "down")
    if(length(index_ovp[[i]]) == 0){
      # find the target region in the all the transcripts
      
      index = queryHits(suppressWarnings(findOverlaps(genome_info$trx, G_target[[i]])))
      if(length(index) !=0){
        index = genome_info$trx_gene_AA[index,] %>% arrange(desc(aa_length), desc(trx_length)) %>% pull(TXID)
        options(warn = -1)
        target_annot = data.frame(brk_annot(ref_trx = genome_info$trx_gene_AA, index = index[1], target_region = G_target[[i]],
                                            up_down = up_down, genome_info = genome_info, BSgenome = BSgenome),
                                  RE_category = re_strand[i], orientation = re_orient[i])
        options(warn = 0)
        # no overlapping and SGL
      } else if(as.character(seqnames(G_target[[i]])) == "chr-1"){
        target_annot = data.frame(G_target[[i]])[,1:2] %>% mutate(gene_id = NA,
                                                                  gene_symbol = NA,
                                                                  trx_name = NA,
                                                                  aa_length = NA,
                                                                  trx_length = NA,
                                                                  loc_in_trx = NA,
                                                                  base_left = NA,
                                                                  AA_seq = "no partner",
                                                                  AA_length = NA,
                                                                  RE_category = "SGL",
                                                                  orientation = target_orient)
        
      } else {
        target_annot = data.frame(G_target[[i]])[,1:2] %>% mutate(gene_id = NA,
                                                                  gene_symbol = NA,
                                                                  trx_name = NA,
                                                                  aa_length = NA,
                                                                  trx_length = NA,
                                                                  loc_in_trx = NA,
                                                                  base_left = NA,
                                                                  AA_seq = "no partner",
                                                                  AA_length = NA,
                                                                  RE_category = "intergenic",
                                                                  orientation = target_orient)
      }
      
    } else {
      index = index_ovp[[i]][(genome_info$trx_gene_AA_major[index_ovp[[i]],] %>% mutate(ind = 1:length(index_ovp[[i]])) %>% arrange(desc(aa_length), desc(trx_length)) %>% pull(ind))]
      target_annot = data.frame(brk_annot(ref_trx = genome_info$trx_gene_AA_major, index = index[1], target_region = G_target[[i]], 
                                          up_down = up_down, genome_info = genome_info, BSgenome = BSgenome),
                                RE_category = re_strand[i], orientation = re_orient[i])
    }
    
    target_annot_all = rbind(target_annot_all, target_annot)
  }
  
  return(target_annot_all)
}


# function to define the breakpoint
brk_orientation = function(ref_trx, index, target_orient, up_down){
  target_gene_trx = ref_trx[index,]
  
  # target region is upstream of the RE
  if(up_down == "up"){
    
    if((as.character(strand(target_gene_trx)) == "+" & target_orient == 1) | (as.character(strand(target_gene_trx)) == "-" & target_orient == (-1))){
      re_categ = "in-strand"
    } else {
      re_categ = "out-of-strand"
    }
    # target region is downstream of the RE
  } else if (up_down == "down") {
    if((as.character(strand(target_gene_trx)) == "+" & target_orient == (-1)) | (as.character(strand(target_gene_trx)) == "-" & target_orient == 1)){
      re_categ = "in-strand"
    } else {
      re_categ = "out-of-strand"
    }
  } else {
    stop("'up_down' parameter should be either 'up' or 'down'")
  }
  return(re_categ)
}

# function to define the breakpoint
brk_annot = function(ref_trx, index, target_region, up_down, genome_info, BSgenome){
  
  target_gene_trx = ref_trx[index,]
  target_gene_trx = target_gene_trx %>% mutate(Protein_coding = ifelse(!is.na(aa_length), TRUE, FALSE))
  
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
  
  if(target_gene_trx$Protein_coding == T){
    G_cds_target = genome_info$G_cds[[which(names(genome_info$G_cds) == target_gene_trx$TXID)]]
    target_cds_seq_AA = as.character(genome_info$cds_seqs_aa[names(genome_info$cds_seqs_aa) == target_gene_trx$TXID][[1]])
    
    if(up_down == "up"){
      
      # breakpoint is in intronic region
      cds_target_region = G_cds_target[queryHits(findOverlaps(query = G_cds_target, subject = G_ref[1:ind]))]
      if(grepl("Exon", loc_in_trx$id)){
        if(length(queryHits(findOverlaps(query = G_cds_target, subject = target_region)))!=0){ # breakpoint in UTR region in the exon covering both UTR and coding region
          if(as.character(strand(G_ref))[1] == "+"){
            end(ranges(cds_target_region[length(cds_target_region)])) = end(ranges(target_region))
          } else {
            start(ranges(cds_target_region[length(cds_target_region)])) = start(ranges(target_region))
          }
        } else {
          cds_target_region = NULL
        }
      }
      
    } else { # downstream of the breakpoint
      
      cds_target_region = G_cds_target[queryHits(findOverlaps(query = G_cds_target, subject = G_ref[ind:length(G_ref)]))]
      if(grepl("Exon", loc_in_trx$id)){
        if(length(queryHits(findOverlaps(query = G_cds_target, subject = target_region)))!=0){ # breakpoint in UTR region in the exon covering both UTR and coding region
          if(as.character(strand(G_ref))[1] == "+"){
            start(ranges(cds_target_region[1])) = start(ranges(target_region))
          } else {
            end(ranges(cds_target_region[1])) = end(ranges(target_region))
          }
        } else {
          cds_target_region = NULL
        }
      }
    }
    
    # breakpoint at UTR region
    if(length(cds_target_region)!=0){
      # AA sequence before (if target is upstream) or after (if target is downstream) the breakpoint
      cds_target_region_DNA = do.call(xscat, getSeq(BSgenome, cds_target_region))
      cds_target_region_DNA_res = nchar(as.character(cds_target_region_DNA)) %% 3
      
      options(warn = -1)
      if(cds_target_region_DNA_res == 0 | up_down == "up"){
        cds_target_region_AA = translate(cds_target_region_DNA)
      } else {# if downstream of a breakpoint
        cds_target_region_AA = translate(cds_target_region_DNA[-c(1:cds_target_region_DNA_res)])
      }
      options(warn = 0)
      
      if(grepl(cds_target_region_AA, target_cds_seq_AA, fixed = T)){
        annot_target_region = data.frame("base_left" = cds_target_region_DNA_res, "AA_seq" = as.character(cds_target_region_AA), "AA_length" = nchar(cds_target_region_AA))
      } else {
        annot_target_region = data.frame("base_left" = cds_target_region_DNA_res, "AA_seq" = "not_consistent with reference", "AA_length" = nchar(cds_target_region_AA))
      }
      
    } else {
      annot_target_region = data.frame("base_left" = NA, "AA_seq" = "UTR", "AA_length" = NA)
    }
    
  } else { # breakpoint at non-coding gene
    annot_target_region = data.frame("base_left" = NA, "AA_seq" = "non-coding RNA", "AA_length" = NA)
  }
  
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
                                                            base_left = annot_target_region$base_left,
                                                            AA_seq = annot_target_region$AA_seq,
                                                            AA_length = annot_target_region$AA_length)
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
    
    base_left_sum = purple_annot$up.base_left + purple_annot$down.base_left
    if((base_left_sum %% 3) == 0 & !is.na(base_left_sum)){
      # if both up and downstream sequences are in-frame & up and down breakpoints are both exons or introns
      n_exons = sum(c(grepl("Exon", purple_annot$up.loc_in_trx), grepl("Exon", purple_annot$down.loc_in_trx)))
      # (in case up is in exon and down is in intron, it is frame-unknown)
      if(n_exons !=1){ 
        RE = "in-frame"
      } else {
        RE = "frame-unknown"
      }
    } else {
      RE = "frame-unknown"
    }
  } else if((purple_annot$up.RE_category %in% "out-of-strand") | (purple_annot$down.RE_category %in% "out-of-strand")){
    RE = "out-of-strand"
  } else {
    RE = "unknown"
  }
  return(RE)
}


##############################################################
# data loading
##############################################################
list_linx = list.files("/DATA/projects/j.bhin/KWF_Daniel/HMF/somatics")
temp = c()
for(i in 1:length(list_linx)){
  t = list.files(paste("/DATA/projects/j.bhin/KWF_Daniel/HMF/somatics", list_linx[i], "linx", sep = "/"))
  t = t[grep("T", t)][1]
  temp = c(temp, stringr::str_split_fixed(t, "[.]", 2)[1])
}
list_linx = data.frame(file = list_linx, sample = temp, stringsAsFactors = F)


##############################################################
# data processing
##############################################################
data_path = "/DATA/projects/j.bhin/KWF_Daniel/HMF/somatics"

target_samp = list.files(data_path)

FGFR3_SV = vector(mode = "list", length(target_samp))
names(FGFR3_SV) = target_samp
options(warn=2)
for(i in 1:length(target_samp)){
  ######################
  # linx fusion & SV clusters
  ######################
  # fusions defined by linx
  t = list.files(paste(data_path, target_samp[i], "linx", sep = "/"))
  linx_fus = read.table(paste(data_path, target_samp[i], "linx", t[grep("linx.fusion.tsv", t)], sep = "/"), sep = "\t", stringsAsFactors = F, header = T)
  linx_fus_FGFR3 = linx_fus %>% filter(grepl("FGFR3", name))
  
  # breakends corresponding to SVID of target transcript
  brkends = read.table(paste(data_path, target_samp[i], "linx", t[grep("linx.breakend.tsv", t)], sep = "/"), sep = "\t", stringsAsFactors = F, header = T)
  # for some cases, first column name is "id" not "Id"
  t_brkends = brkends %>% filter(gene == "FGFR3")
  brkends_FGFR3 = t_brkends %>%
    mutate(sample_name = target_samp[i], svid_samp = paste(svId, sample_name, sep = "_"), id_samp = paste(t_brkends[,1], sample_name, sep = "_"))
  
  # extracting entire chains including the corresponding SVs
  sv_data = read.table(paste(data_path, target_samp[i], "linx", t[grep("vis_sv_data.tsv", t)], sep = "/"), sep = "\t", stringsAsFactors = F, header = T)
  sv_data_FGFR3 = sv_data %>% filter(SvId %in% brkends_FGFR3$svId) %>% 
    mutate(sample_name = target_samp[i], svid_samp = paste(SvId, sample_name, sep = "_"),
           Cluster_Chain = paste(ClusterId, ChainId, sep = "_"))
  
  ########################################
  # purple
  ########################################
  t = list.files(paste(data_path, target_samp[i], "purple", sep = "/"))
  purple_SV = read.table(gzfile(paste(data_path, target_samp[i], "purple",  t[grepl("purple.sv.vcf.gz", t) & !grepl("tbi", t)], sep = "/")), sep = "\t", stringsAsFactors = F, header = F)
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
      
      # allele frequency of SVs
      t_out_gene = t_out %>% mutate(gene_symbol = replace(gene_symbol, is.na(gene_symbol), "NA")) %>% pull(gene_symbol)
      if(grepl("PURPLE_AF", purple_SV_FGFR3$V8[j])){
        t_af = as.numeric(strsplit(grep("PURPLE_AF", strsplit(purple_SV_FGFR3$V8[j], ";")[[1]], value = T), "=|,")[[1]][-1])
        if(t_out_gene[1] == "FGFR3" & t_out_gene[2] == "FGFR3"){
          t_af = mean(t_af)
        } else {
          t_af = t_af[grep("FGFR3", t_out_gene)]
        }
      } else {
        t_af = NA
      }
      t_out = data.frame(up = t_out[1,], down = t_out[2,], PURPLE_AF = t_af)
      purple_SV_annot = rbind(purple_SV_annot, t_out)
    }
    
    # remove redundant SVs
    purple_SV_annot = purple_SV_annot %>% arrange(up.RE_category)
    t_sv = t(apply(purple_SV_annot %>% select(up.start, down.start), 1, sort))
    t_id = apply(t_sv, 1, paste, collapse = "-")
    purple_SV_annot = purple_SV_annot[match(unique(t_id), t_id),]
    
    # check whether the RE was detected by LINX
    t_la = do.call("paste", c(purple_SV_annot[, c("up.gene_symbol", "down.gene_symbol")], sep = "_"))
    t_linx_lab = c(); t_sf_lab = c()
    for(j in 1:length(t_la)){
      t_la_2 = t_la[j]
      
      # LINX-Fusion
      if(!("NA" %in% (unlist(strsplit(t_la_2, "_"))))){
        # LINX-Fusion
        if(!is.na(match(t_la_2, linx_fus_FGFR3$name))){
          t_linx_lab_2 = do.call("paste", c(linx_fus[match(t_la_2, linx_fus$name), c("name", "fusedExonUp", "fusedExonDown")], sep = "_"))
        } else {
          t_linx_lab_2 = NA
        }
      } else {
        t_linx_lab_2 = NA
      }
      t_linx_lab = c(t_linx_lab, t_linx_lab_2)
    }
    purple_SV_annot$LINX_Fusion = t_linx_lab
    
  } else {
    purple_SV_annot = data.frame(NULL)
  }
  
  ##################################################
  # check the chain of the identified SVs
  ##################################################
  if(nrow(purple_SV_annot)!=0){
    pos_id = apply(purple_SV_annot[, c(2, 15)], 1, function(x){paste(sort(as.numeric(x)), collapse = "_")})
    sv_data_FGFR3_id = apply(sv_data_FGFR3[, 10:11], 1, function(x){paste(sort(as.numeric(x)), collapse = "_")})
    purple_linx = matrix(NA, nrow(purple_SV_annot), 5)
    for(k in 1:length(pos_id)){
      ind = which(sv_data_FGFR3_id %in% pos_id[k])
      if(length(ind)==0){
        ind = c(which(sv_data_FGFR3$PosStart %in% purple_SV_annot$up.start[k]),
                which(sv_data_FGFR3$PosEnd %in% purple_SV_annot$up.start[k]),
                which(sv_data_FGFR3$PosStart %in% purple_SV_annot$down.start[k]),
                which(sv_data_FGFR3$PosEnd %in% purple_SV_annot$down.start[k]))
        ind = ind[!is.na(ind)]
      }
      t_id = sv_data_FGFR3[ind, c("JunctionCopyNumber", "ClusterId", "ChainId", "SvId", "Cluster_Chain")]
      purple_linx[k,] = c(t_id$JunctionCopyNumber[1], t_id$ClusterId[1], paste(t_id$ChainId, collapse = "/"), t_id$SvId[1], paste(t_id$Cluster_Chain, collapse = "/"))
    }
    purple_SV_annot = purple_SV_annot %>% mutate(Ploidy = as.numeric(purple_linx[,1]), ClusterID = as.numeric(purple_linx[,2]), 
                                                 ChainID = purple_linx[,3], SVID = as.numeric(purple_linx[,4]),
                                                 Cluster_Chain = purple_linx[,5]) %>%
      mutate(ChainID = ifelse(ChainID == "", NA, ChainID), Cluster_Chain = ifelse(Cluster_Chain == "", NA, Cluster_Chain))
  }
  
  FGFR3_SV[[i]] = list(linx_fus = linx_fus_FGFR3, purple_SV_annot = purple_SV_annot)
  print(i)
}
FGFR3_SV = FGFR3_SV[which(unlist(lapply(FGFR3_SV, function(x){nrow(x[[2]])}))!=0)]
target_samp = names(FGFR3_SV)

FGFR3_SV_n = c()
for(i in 1:length(FGFR3_SV)){
  t = unlist(lapply(FGFR3_SV[[i]], nrow))[c(1,2)]
  FGFR3_SV_n = rbind(FGFR3_SV_n, t)
}
rownames(FGFR3_SV_n) = names(FGFR3_SV)


################################################################################
# combined results from PURPLE
################################################################################
ind_tar = which(rowSums(FGFR3_SV_n) >0)
FGFR3_PURPLE_ALL = c()
for(i in 1:length(ind_tar)){
  
  # Linx fusion
  t_linx = FGFR3_SV[[ind_tar[i]]]$linx_fus %>% arrange(desc(fusedExonUp)) %>% select(name, fusedExonUp, fusedExonDown)
  
  t_RE = c()
  for(j in 1:nrow(FGFR3_SV[[ind_tar[i]]]$purple_SV_annot)){
    t_RE = c(t_RE, RE_type(FGFR3_SV[[ind_tar[i]]]$purple_SV_annot[j,]))
  }
  
  # Purple SV
  t_purple = FGFR3_SV[[ind_tar[i]]]$purple_SV_annot %>% 
    mutate(C_trunc = (grepl("Exon 18|Intron 17", up.loc_in_trx) & up.RE_category == "in-strand" & up.gene_symbol == "FGFR3")) %>%
    mutate(RE_type = factor(t_RE, levels = c("in-frame", "frame-unknown", "intergenic", "out-of-strand", "internal", "SGL")), 
           FGFR3_is_upstream = (up.gene_symbol == "FGFR3")) %>%
    mutate(sample_name = target_samp[ind_tar[i]], sample_id = ind_tar[i]) %>%
    arrange(desc(C_trunc), desc(FGFR3_is_upstream), RE_type)
  FGFR3_PURPLE_ALL = rbind(FGFR3_PURPLE_ALL, t_purple)
}

FGFR3_PURPLE_ALL$promiscuous = FGFR3_PURPLE_ALL$sample_id %in% names(which(table(FGFR3_PURPLE_ALL$sample_id)>10))

################################################################################
# combined results from LINX
################################################################################
FGFR3_LINX_ALL = c()
for(i in 1:length(ind_tar)){
  t_linx = FGFR3_SV[[ind_tar[i]]]$linx_fus
  if(is.null(t_linx)){
    t_linx = data.frame(NULL)
  }
  
  if(nrow(t_linx)!=0){
    t_linx = t_linx %>% arrange(desc(fusedExonUp)) %>% select(name, fusedExonUp, fusedExonDown)
    t_linx = do.call("paste", c(t_linx, sep = "_"))
    
    t_purple = FGFR3_SV[[ind_tar[i]]]$purple_SV_annot
    
    t_df = c()
    for(j in 1:length(t_linx)){
      # comparing with PURPLE
      if(nrow(t_purple)!=0){
        t_lab = t_linx[j]
        t_df_2 = data.frame(t_lab, t_purple[match(t_linx[j], t_purple$LINX_Fusion), c(8:11, 20:23)])
      } else {
        t_df_2 = data.frame(t_linx[j], rep(NA, 8))
      }
      colnames(t_df_2) = c("LINX_Fusion", "up.loc_in_trx", "up.down_exon", "up.down_exon_frame", "up.RE_category", 
                           "down.loc_in_trx", "down.down_exon", "down.down_exon_frame", "down.RE_category")
      
      t_df_2 = data.frame(sample_id = ind_tar[[i]], t_df_2)
      t_df = rbind(t_df, t_df_2)
    }
    FGFR3_LINX_ALL = rbind(FGFR3_LINX_ALL, t_df)
  }
}
save(target_samp, FGFR3_SV, FGFR3_SV_n, hg19_info, FGFR3_PURPLE_ALL, FGFR3_LINX_ALL, file = "~/FGFR3/nvim_R/HMF_data_purple_FGFR3.RData")

