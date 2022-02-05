library(TCGA2STAT)
library(plyr)
library(DescTools)
library(BSDA)
tumortypes <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIPAN","KIRC","KIRP",
              "LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD",
              "TGCT","THCA","THYM","UCEC","UCS","UVM")
fdr.cutoff <- 0.05
for (t in 1:length(tumortypes)){
  tumor=tumortypes[t]
  rnaseqdata <- getTCGA(disease=tumor, data.type="RNASeq2")
  rsem.log <- log2(rnaseqdata$dat + 0.01)
  granzyme <- rnaseqdata$dat[which(rownames(rnaseqdata$dat)=="GZMA"),] + 0.01
  perforin <- rnaseqdata$dat[which(rownames(rnaseqdata$dat)=="PRF1"),] + 0.01
  cyt <- log2(sqrt(granzyme*perforin))
  nsamples <- length(cyt)
  p.vals <- NULL
  corr.vals <- NULL
  for (i in 1:nrow(rsem.log)){
    exp <- rsem.log[i,]
    corr <- cor.test(cyt,exp, method="pearson")
    rval <- corr$estimate
    if (is.na(rval)){
      rval <- 0
    }
    corr.vals <- rbind(corr.vals, rval)
    rprime <- FisherZ(rval)
    zval <- rprime*sqrt(nsamples-3)
    if (zval > 0){
      p.vals <- rbind(p.vals, (1-pnorm(zval))*2)
    }
    else{
      p.vals <- rbind(p.vals, pnorm(zval)*2)
    }
  }
  p.vals[is.na(p.vals)] <- 1
  fdr <- p.adjust(p.vals, "BH")
  corrdf <- as.data.frame(cbind(rownames(rsem.log), p.vals, fdr, corr.vals))
  colnames(corrdf) <- c("Gene", "pvals", "FDR", "corr")
  write.csv(corrdf, file=paste0("TCGA_", tumor, "_CYT_pearson_fdr.csv"), row.names=FALSE)
}
