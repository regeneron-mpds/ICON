#'######################################################################################################
#' Copyright 2021 Regeneron Pharmaceuticals Inc. All rights reserved.
#' License for Non-Commercial Use of ICON code
#' All files in this repository (“source code”) are licensed under the following terms below:
#' “You” refers to an academic institution or academically employed full-time personnel only.
#' “Regeneron” refers to Regeneron Pharmaceuticals, Inc.
#' Regeneron hereby grants You a right to use, reproduce, modify, or distribute the source code to the TCRAI and ICON algorithms, in whole or in part, whether in original or modified form, for academic research purposes only.  The foregoing right is royalty-free, worldwide, revocable, non-exclusive, and non-transferable.
#' Prohibited Uses:  The rights granted herein do not include any right to use by commercial entities or commercial use of any kind, including, without limitation, any integration into other code or software that is used for further commercialization, any reproduction, copy, modification or creation of a derivative work that is then incorporated into a commercial product or service or otherwise used for any commercial purpose, or distribution of the source code not in conformity with the restrictions set forth above, whether in whole or in part and whether in original or modified form, and any such commercial usage is not permitted.
#' Except as expressly provided for herein, nothing in this License grants to You any right, title or interest in and to the intellectual property of Regeneron (either expressly or by implication or estoppel).  Notwithstanding anything else in this License, nothing contained herein shall limit or compromise the rights of Regeneron with respect to its own intellectual property or limit its freedom to practice and to develop its products and product candidates.
#' If the source code, whole or in part and in original or modified form, is reproduced, shared or distributed in any manner, it must (1) identify Regeneron Pharmaceuticals, Inc. as the original creator, and (2) include the terms of this License.
#' UNLESS OTHERWISE SEPARATELY AGREED UPON, THE SOURCE CODE IS PROVIDED ON AN AS-IS BASIS, AND REGENERON PHARMACEUTICALS, INC. MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE SOURCE CODE, IN WHOLE OR IN PART AND IN ORIGINAL OR MODIFIED FORM, WHETHER EXPRESS, IMPLIED, STATUTORY, OR OTHER REPRESENTATIONS OR WARRANTIES. THIS INCLUDES, WITHOUT LIMITATION, WARRANTIES OF TITLE, MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NON-INFRINGEMENT, ABSENCE OF LATENT OR OTHER DEFECTS, ACCURACY, OR THE PRESENCE OR ABSENCE OF ERRORS, WHETHER OR NOT KNOWN OR DISCOVERABLE.
#' In no case shall Regeneron be liable for any loss, claim, damage, or expenses, of any kind, which may arise from or in connection with this License or the use of the source code. You shall indemnify and hold Regeneron and its employees harmless from any loss, claim, damage, expenses, or liability, of any kind, from a third-party which may arise from or in connection with this License or Your use of the source code.
#' You agree that this License and its terms are governed by the laws of the State of New York, without regard to choice of law rules or the United Nations Convention on the International Sale of Goods.
#' Please reach out to Regeneron Pharmaceuticals Inc./Administrator relating to any non-academic or commercial use of the source code.
#'######################################################################################################

#' Normalize dextramer signals and identify antigen specific TCRs
#'
#' @param Dex_UMI a dextramer UMI matrix of QC passed T cells
#' @param TCR_Pairs a data frame of T cells with paired alpha and beta chains
#' @param Bg_cutoff threshold of estimated background
#' @param out_dir path to the output/result directory
#'
#' @return a matrix of antigen specific TCRs
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @importFrom stats density
#' @importFrom stats heatmap
#' @importFrom utils write.csv
#' @export
Identify_binders <- function(Dex_UMI = NULL, TCR_Pairs = NULL, Bg_cutoff = 10, out_dir = NULL) {
  #Dex_UMI: matrix of test dextramers
  if(is.matrix(Dex_UMI) == FALSE) stop("'Dex_UMI' must be a data matrix")
  if(is.data.frame(TCR_Pairs) == FALSE) stop("'TCR_Pairs' must be a data frame")
  if(is.numeric(Bg_cutoff) == FALSE) stop("'Bg_cuoff' must be numeric")
  if (missing(out_dir))
    stop("Need to input a path to a result directory")

  ### subtract background noise
  message("Subtracting background noise...")
  sorted.clean <- as.matrix(Dex_UMI) - Bg_cutoff
  sorted.clean[which(sorted.clean < 0)] <- 0

  #### Estimat the ratio of a dextramer signal within a cell to the
  #### sum of all dextramers binding to the cell
  message("Estimating RC...")
  sum.sorted.clean<- colSums(sorted.clean)
  ratio<- sweep(sorted.clean, 2, sum.sorted.clean+1, FUN = '/')
  sorted.clean <- gdata::unmatrix(sorted.clean)
  ratio <- gdata::unmatrix(ratio)
  cell <- sapply(strsplit(names(sorted.clean), ":"), "[", 2)
  pMHC <- sapply(strsplit(names(sorted.clean), ":"), "[", 1)
  sorted.clean <- data.frame(cell, pMHC, as.vector(sorted.clean), as.vector(ratio))
  colnames(sorted.clean)[3:4]<-c("UMI", "RC")
  sorted.clean <- sorted.clean[sorted.clean$cell %in% TCR_Pairs[,1],]

  #### Map to the T cells with paired alpha and beta chains
  tcr.map <- TCR_Pairs[match(sorted.clean$cell, TCR_Pairs$barcode),]
  tcr.map<-data.frame(tcr.map, sorted.clean[,2:4])
  tcr.map$unique<-paste0(tcr.map[,3],tcr.map[,4],tcr.map[,5],tcr.map[,7],tcr.map[,8],tcr.map[,9])

  #### Estimate the fraction of T cells within a clone binding to a particular dextramer
  message("Estimating RT...")
  temp<-tcr.map[which(tcr.map$UMI>0),]
  uniq.tcr <- unique(temp$unique)
  TCR.clone <- rep(0,length(uniq.tcr))
  for(m in 1:length(uniq.tcr)){
    peptide.all <- temp[which(temp$unique == uniq.tcr[m]), 10]
    peptide.uniq <- unique(peptide.all)
    TCR.clone[m] <- length(peptide.all)
    sub.count <- rep(0,length(peptide.uniq))
    if(length(peptide.uniq) > 1){
      for(n in 1:length(peptide.uniq)){
        sub.count[n] <- length(which(peptide.all == peptide.uniq[n]))
      }
    }
    else{sub.count[1] <- 1}
    Pout<-data.frame(uniq.tcr[m], peptide.uniq, TCR.clone[m], sub.count)
    colnames(Pout)<-c("TCR", "pMHC", "TCR_clone_size", "TCR_subClone_size")
    if(m == 1){portion <- Pout}
    else{portion <- rbind(portion, Pout)}
  }
  portion <- cbind(paste(portion$TCR, portion$pMHC), portion)
  tcr.map <- cbind(paste(tcr.map$unique, tcr.map$pMHC), tcr.map)
  portion <- portion[match(tcr.map[,1], portion[,1]),]
  tcr.map <- cbind(tcr.map, portion[,4:5])
  tcr.map$TCR_clone_size[which(is.na(tcr.map$TCR_clone_size))]<-0
  tcr.map$TCR_subClone_size[which(is.na(tcr.map$TCR_subClone_size))]<-0
  tcr.map$RT <- tcr.map[,16]/tcr.map[,15]
  tcr.map$RT[which(is.na(tcr.map$RT))] <- 0
  tcr.map <- tcr.map[,c(-1,-15,-16)]

  #### Dextramer signal correction
  message("Correcting dextramer signal...")
  cor.UMI<-log(tcr.map[,11]+0.01)*tcr.map[,12]^2*tcr.map[,14]
  cor.UMI[which(is.na(cor.UMI))]<-0
  cor.UMI[which(cor.UMI<1)]<-0
  tcr.map<-cbind(tcr.map,cor.UMI)
  rm(temp)

  #### Cell- and pMHC-wise normalization
  message("Normalizing data...")
  normalized<-reshape2::acast(tcr.map, pMHC ~ barcode, value.var = "cor.UMI")
  sum.tcr.map <- colSums(normalized)
  sum.tcr.map[which(sum.tcr.map==0)]<-1
  normalized <- sweep(normalized, 2, sum.tcr.map, FUN = '/')
  normalized<-t(normalized)
  normalized <- scale(normalized)
  normalized <- t(normalized)
  temp<-gdata::unmatrix(normalized)
  smallest<-min(temp[which(!is.na(temp))])
  normalized[which(is.na(normalized))]<-smallest
  temp[which(is.na(temp))]<-smallest
  saveRDS(normalized,file=paste0(out_dir,"/Normalized_data.RDS"))
  tcr.map <- cbind(paste0(tcr.map$pMHC,":",tcr.map$barcode), tcr.map)
  temp<-temp[match(tcr.map[,1],names(temp))]
  tcr.map<-cbind(tcr.map,as.vector(temp))
  tcr.map<-tcr.map[,-1]
  tcr.map<-tcr.map[,c(1:9,13,10:12,14:16)]
  colnames(tcr.map)<-c("barcode", "chainType_1", "beta_v_gene", "beta_j_gene",
                       "TCRB_cdr3", "chainType_2", "alpha_v_gene", "alpha_j_gene",
                       "TCRA_cdr3", "unique_TCR", "pMHC", "UMI","RC","RT",
                       "Bg_corrected_data","normalized_data")
  #saveRDS(tcr.map,file="ICON_output.RDS")

  #### Identify and report binders
  message("Identifying binders...")
  pdf(paste0(out_dir,"/normalized_data_distribution.pdf"))
  temp[which(temp<0)]<-0
  den<-density(log(temp+0.1,2))
  plot(den,main="normalized_data_distribution",ylim=c(0,0.5))
  dev.off()
  cutoff = 0 ## user defined
  out <- tcr.map[which(tcr.map[,16] > cutoff),]
  write.csv(out, file = paste0(out_dir,"/ICON_identified_binders.csv"),row.names = FALSE)

  #### Generate stats report
  message("Generating binder stats report...")
  pMHC <- unique(out$pMHC)
  total <- rep(0,length(pMHC))
  unique <- rep(0,length(pMHC))
  for(k in 1:length(pMHC)){
    total.tcr <- out[which(out[,11] == pMHC[k]), 10]
    total[k] <- length(total.tcr)
    unique[k] <- length(unique(total.tcr))
  }
  stats <- data.frame(pMHC, total, unique)
  write.csv(stats,file=paste0(out_dir,"/Binder_summary.csv"),row.names=FALSE)

  #### plot binder heatmap
  message("Plotting binder heatmap...\n")
  pdf(paste0(out_dir,"/binder_heatmap.pdf"))
  normalized[which(normalized<0)]<-0
  normalized<-normalized[,which(colSums(normalized)!=0)]
  normalized<-normalized[which(rowSums(normalized)!=0),]
  rbPal<- colorRampPalette(c("white","#cb181d","#99000d"))
  high<-log(max(temp),2)
  colors<-rbPal(7)[cut(c(0:high), breaks=7)]
  heatmap(log(normalized+0.01,2), col = colors, scale = "none",margins = c(5,10))
  dev.off()

  return(out)
}
