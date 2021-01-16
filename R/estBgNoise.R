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

#' Plot background noise and signal distributions using negative control and test dextramers
#' Return max dextramer UMI matrixes for both test dextramers and negative controls
#'
#' @param dex a data frame or matrix of test dextramers signals in UMI
#' @param nc a data frame or matrix of negative control dextramer signals in UMI
#' @param plot logic flag for plotting the UMI distributions
#'
#' @return a list of maximum UMIs of test dextramers and negative controls per cell
#' @importFrom graphics legend
#' @importFrom graphics lines
#' @importFrom graphics plot
#' @importFrom stats density
#' @export
estBgNoise <- function (dex = NULL,nc = NULL,plot = TRUE){
  if(is.matrix(dex) == FALSE) stop("'dex' must be a dataframe or matrix")
  if(is.matrix(nc) == FALSE) stop("'nc' must be a dataframe or matirx")
  colMax <- function(data) sapply(data, max, na.rm = TRUE)

  max.dex<-colMax(as.data.frame(dex))
  max.nc<-colMax(as.data.frame(nc))
  den.max.dex<-density(log(max.dex+1,10))
  den.max.nc<-density(log(max.nc+1,10))

  if(plot){
    plot(den.max.nc, ylim=c(0,2),xlim=c(0,3),col="red",main="Distribution of Dextramer UMI",cex.main=0.8,cex.lab=0.8, cex.axis=0.8, cex.sub=0.8,xlab="Log(Dextramer UMI+1, 10)")
    lines(den.max.dex,ylim=c(0,2),xlim=c(0,3),col="blue")
    legend("topright", legend = c("Neg_control","Test_dex"),text.col=c("red","blue"))
   }
  data_list<-list(den.max.dex,den.max.nc)
  return(data_list)
}
