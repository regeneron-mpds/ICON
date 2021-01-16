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

#' Select T cells with paired alpha and beta chains
#'
#' @param filtered_contig a dataframe containing TCR sequence contig information
#'
#' @return a data frame of T cells with paired alpha and beta chains
#' @export
pair_TCR <- function(filtered_contig = NULL) {
  if(is.data.frame(filtered_contig) == FALSE) stop("filtered_contig is a dataframe")

  filtered_contig<-filtered_contig[which(toupper(filtered_contig$high_confidence)=="TRUE"&toupper(filtered_contig$productive)=="TRUE"),]
  filtered_contig<-filtered_contig[,c(1,6,7,9,13,15)]
  entity<- length(filtered_contig[,1])
  message(paste("Number of entity (alpha + beta):",entity))

  tcrb<-filtered_contig[which(filtered_contig[,2]=="TRB"),]
  Dup.tcrb<-tcrb[which(duplicated(tcrb[,1])),1]
  single.tcrb<-tcrb[!tcrb[,1] %in% Dup.tcrb,]
  Dup.tcrb<-tcrb[tcrb[,1] %in% Dup.tcrb,]
  tcra<-filtered_contig[which(filtered_contig[,2]=="TRA"),]
  Dup.tcra<-tcra[which(duplicated(tcra[,1])),1]
  single.tcra<-tcra[!tcra[,1] %in% Dup.tcra,]
  Dup.tcra<-tcra[tcra[,1] %in% Dup.tcra,]

  tcrb.cell<-unique(Dup.tcrb[,1])
  tcra.cell<-unique(Dup.tcra[,1])

  message("Process alpha chain...")
  for(i in 1:length(tcra.cell)){
    #print(paste("tcra",i))
    temp<-Dup.tcra[which(Dup.tcra[,1]==tcra.cell[i]),]
    temp<-temp[which(temp[,6]==max(temp[,6])),]
    if(i==1){out.tcra<-temp}
    else{out.tcra<-rbind(out.tcra,temp)}
  }

  message("Process beta chain...")
  for(i in 1:length(tcrb.cell)){
    #print(paste("tcrb",i))
    temp<-Dup.tcrb[which(Dup.tcrb[,1]==tcrb.cell[i]),]
    temp<-temp[which(temp[,6]==max(temp[,6])),]
    if(i==1){out.tcrb<-temp}
    else{out.tcrb<-rbind(out.tcrb,temp)}
  }
  tcrb<-rbind(single.tcrb,out.tcrb)
  tcra<-rbind(single.tcra,out.tcra)
  pair.tcrb<-tcrb[tcrb[,1] %in% tcra[,1],]
  pair.tcra<-tcra[tcra[,1] %in% pair.tcrb[,1],]
  pair.tcra<-pair.tcra[match(pair.tcrb[,1],pair.tcra[,1]),]
  paired<-data.frame(pair.tcrb,pair.tcra)
  paired<-paired[,-c(6,7,12)]
  colnames(paired)<-c("barcode","b_chain","b_v_gene","b_j_gene","b_cdr3","a_chain","a_v_gene","a_j_gene","a_cdr3")
  return(paired)
}
