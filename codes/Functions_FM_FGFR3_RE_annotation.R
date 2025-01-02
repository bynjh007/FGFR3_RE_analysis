library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

hg_id = read.table("~/Dropbox/Work_Netherlands/FGFR3/R/data/UCSC.hg.19_ID_to_SYMBOL.txt", header = T, stringsAsFactors = F, sep = "\t", comment.char = "")

# required information to run the function
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
trx = transcripts(txdb)
trx_intron = intronsByTranscript(txdb)
trx_exon = exonsBy(txdb)

trx_gene = select(txdb, key = as.character(trx$tx_id), keytype = "TXID", columns = c("TXNAME", "GENEID"))
# TXID is same with index of the "trx", "trx_gene"

annots <- select(org.Hs.eg.db, keys=trx_gene$GENEID, columns=c("SYMBOL","GENENAME"), keytype="ENTREZID")
trx_gene$SYMBOL = hg_id$hg19.kgXref.geneSymbol[match(trx_gene$TXNAME, hg_id$X.hg19.knownGene.name)]
trx_gene$STRAND = as.character(strand(trx))

# coding DNA sequence
G_cds = cdsBy(txdb, by="tx")
cds_seqs = extractTranscriptSeqs(Hsapiens, cdsBy(txdb, by="tx"))
cds_seqs_aa = translate(cds_seqs)

library(dplyr)

# length of the amino acid and transcript 
trx_gene_AA = trx_gene %>% mutate(aa_length = width(cds_seqs_aa)[match(TXID, names(cds_seqs_aa))],
                                  trx_length = width(trx),
                                  GENEID = replace(GENEID, is.na(GENEID), TXNAME[is.na(GENEID)])) %>% 
  arrange(desc(aa_length), desc(trx_length))

total_gene_list = unique(trx_gene_AA$GENEID)

trx_gene_AA_major = trx_gene_AA[match(total_gene_list, trx_gene_AA$GENEID),]
trx_gene_AA_major = trx_gene_AA_major %>% arrange(TXID)

trx_major = trx[trx_gene_AA_major$TXID]
trx_gene_AA = trx_gene_AA %>% arrange(TXID)

# function
re_annot = function(brkpt){
  
  target_gene = stringr::str_split_fixed(brkpt, " ", 2)[1]
  target_loc_1 = stringr::str_split_fixed(brkpt, " ", 2)[2]
  target_loc = gsub("_UTR", "", target_loc_1)
  
  target_gene_trx = trx_gene_AA_major %>% filter(SYMBOL == target_gene)
  target_gene_trx = target_gene_trx %>% mutate(Protein_coding = ifelse(!is.na(aa_length), TRUE, FALSE))
  
  ################################################## 
  # build the GRange object for the target transcript
  ##################################################
  ind = as.numeric(target_gene_trx$TXID)
  t_exon = trx_exon[[ind]]
  t_intron = trx_intron[[ind]]
  t_str = unique(as.character(strand(t_exon)))
  if(t_str == "-"){
    t_intron = rev(t_intron)
  }
  df1 = data.frame(data.frame(t_exon)[, c(1,2,3,5)], id = paste("exon", 1:length(t_exon)))
  if(length(t_intron)!=0){
    df2 = rbind(data.frame(data.frame(t_intron)[, c(1,2,3,5)], id = paste("intron", 1:length(t_intron))), NA)
  } else {
    df2 = c()
  }
  
  df_comb = (gdata::interleave(df1, df2))
  df_comb = df_comb[!is.na(df_comb$seqnames),]
  
  G_ref = makeGRangesFromDataFrame(df_comb)
  G_ref$id = df_comb$id
  
  #######################################################################################
  # Find the location of the breakpoint and downstream exons within the target transcript
  #######################################################################################
  # Compare the Ranges between the target transcript and target region (breakpoint)
  loc_in_trx = G_ref[G_ref$id == target_loc]
  
  # target region is either exon or intron
  if(length(loc_in_trx) == 1){
    ind = which(G_ref$id == target_loc)
    
    # following exon from the breakpoint
    if(ind == length(G_ref)){
      down_exon = G_ref[length(G_ref)]
    } else {
      down_exon = G_ref[(ind+1):length(G_ref)]
    }
    down_exon = down_exon[grepl("exon", down_exon$id)]
    
    #######################################################################################
    # check whether the corresponding trx is protein-coding or not
    #######################################################################################
    if(target_gene_trx$Protein_coding == T){
      
      if(grepl("intron", target_loc)){
        # protein coding region in the transcript
        G_cds_target = G_cds[[which(names(G_cds) == target_gene_trx$TXID)]]
        # Protein sequence for the transcript
        target_cds_seq_AA = as.character(cds_seqs_aa[names(cds_seqs_aa) == target_gene_trx$TXID][[1]])
        
        #######################################################################################
        # check the following exon is in-frame or out-of-frame or unknown
        #######################################################################################
        # the nearby exon is completely overlapped with exon of full CDS on the transcript (exon in the middle of transcript)
        if(sum(G_cds_target == down_exon[1]) == 1){
          # check whether the protein sequence of the neigbor exon (previous exon for upstream and next exon for downstream) is in-frame or not
          # this can be identified by the aa sequence of the next exon of the breakpoint
          t_seq = as.character(suppressWarnings(translate(getSeq(Hsapiens, G_cds_target[which(G_cds_target == down_exon[1])])))[[1]])
          # target exon is inframe --> fully overlapped with protein sequence of the transcript
          if(grepl(t_seq, target_cds_seq_AA, fixed = T)){
            re_type = "inframe:CDS_full"
          } else {
            re_type = "out_of_frame:CDS_full"
          }
        
          # the following exon is partly overlapped with CDS of the transcript (exons with translation start (5') or end (3'))
        } else if(length(queryHits(findOverlaps(down_exon[1], G_cds_target))) == 1){
        
          # if the following exon include translation start site -> the following exon should include the first CDS exon
          if(length(queryHits(findOverlaps(down_exon[1], G_cds_target[1])))==1){
            re_type = "unknown:CDS_start"
          
          # if the following exon include translation end site -> the following exon should include the last CDS exon
          } else if(length(queryHits(findOverlaps(down_exon[1], G_cds_target[length(G_cds_target)])))==1){
            # protein sequence of the last CDS exon
            t_seq = as.character(translate(getSeq(Hsapiens, G_cds_target[length(G_cds_target)]))[[1]])
            # check whether the last CDS exon is in-frame or not
            if(grepl(t_seq, target_cds_seq_AA)){
              re_type = "inframe:CDS_end"
            } else {
              re_type = "out_of_frame:CDS_end"
            }
          
          } else {
            stop("unexpected REs with partial CDS")
          }
        
        # the following exon is completely UTR
        } else if(length(queryHits(findOverlaps(down_exon[1], G_cds_target))) == 0){
        
          # if there is remaining CDS exons (5'-UTR)
          if(length(queryHits(findOverlaps(down_exon[1], G_cds_target)))!=0){
            re_type = "unknown:UTR_5"
            # if there is no remaining CDS exons (3'-UTR)
          } else {
            re_type = "unknown:UTR_3"
          }
        
        } else {
          stop("unexpected RE scenario")
        }
      
      # breakpoint is located in the exonic region  
      } else {
        re_type = "brkpt_in_exon"
      }
      
    # if target transcript is non-coding  
    } else {
      re_type = "unknown:noncoding trx"
    }
    
    # output
    if(as.character(loc_in_trx$id) != "NA"){
      t_loc = apply(data.frame(loc_in_trx), 1, paste, collapse = ":")
    } else {
      t_loc = "out_of_trx"
    }
  
  # target region spans both exon and intron
  } else if(length(temp_ovp)>1){
    
    if(target_gene_trx$Protein_coding == T){
      re_type = "unknown:spanning exon and intron:coding"
    } else {
      re_type = "unknown:spanning exon and intron:noncoding"
    }
    
    t_loc = paste(apply(data.frame(loc_in_trx), 1, paste, collapse = ":"), collapse = " / ")
    ind = tail(queryHits(temp_ovp), 1)
    
    if(ind != length(G_ref)){
      down_exon = G_ref[(ind + 1):length(G_ref)]
    } else {
      down_exon = G_ref[ind:length(G_ref)]
    }
    # all the downstream exons
    down_exon = down_exon[grepl("Exon", down_exon$id)]
    
  } else {
    stop("target region spans nothing")
  }
  
  # if brkpt is in UTR
  if(grepl("UTR", target_loc)){
    loc_in_trx = paste(re_type, ":UTR")
  }
  
  target_annot = data.frame(loc_in_trx) %>% mutate(gene_id = target_gene_trx$GENEID,
                                                 gene_symbol = target_gene_trx$SYMBOL,
                                                 trx_name = target_gene_trx$TXNAME,
                                                 loc_in_trx = t_loc,
                                                 down_exon = as.character(down_exon$id[1]),
                                                 RE_type = re_type)
  return(target_annot)
}


re_annot_framework = function(G_up, G_down){
  target_annot_all = c()
  for(i in 1:2){
    
    if(length(index_ovp[[i]]) == 0){
      # find the target region in the all the transcripts
      index = queryHits(findOverlaps(trx, G_target[[i]]))
      if(length(index) !=0){
        index = trx_gene_AA[index,] %>% arrange(desc(aa_length), desc(trx_length)) %>% pull(TXID)
        target_annot = re_annot(ref_trx = trx_gene_AA, index = index[1], target_region = G_target[[i]])
      } else {
        target_annot = data.frame(G_target) %>% mutate(gene_id = target_gene_trx$GENEID,
                                                 gene_symbol = target_gene_trx$SYMBOL,
                                                 trx_name = target_gene_trx$TXNAME,
                                                 loc_in_trx = t_loc,
                                                 down_exon = as.character(down_exon$id[1]),
                                                 RE_type = re_type)
      }
    
    } else {
      index = index_ovp[[i]][(trx_gene_AA_major[index_ovp[[i]],] %>% mutate(ind = 1:length(index_ovp[[i]])) %>% arrange(desc(aa_length), desc(trx_length)) %>% pull(ind))]
      target_annot = re_annot(ref_trx = trx_gene_AA_major, index = index[1], target_region = G_target[[i]])
    }
    
    target_annot_all = rbind(target_annot_all, target_annot)
  }
  return(target_annot_all)
}


