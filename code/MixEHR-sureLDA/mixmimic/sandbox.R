setwd('~/Desktop/mixehr-surelda/mixmimic')
library(foreach)

mimic_trainData = read.table('mimic_trainData.txt',sep=' ')
mimic_pats = unique(mimic_trainData[,1])
n = length(mimic_pats); k = 100
mimic_prior = foreach(pat=mimic_pats, .combine=rbind) %do% {
	topics = which(rbinom(k,1,0.1)==1)
	cbind(pat,topics,rbeta(length(topics),1,4))
}
write.table(mimic_prior,file="mimic_prior_compressed.txt",row.names=FALSE,col.names=FALSE)

phecode_icd9_map = read.csv('phecode_icd9_map_unrolled.csv')
phecode_icd9_map = phecode_icd9_map[!duplicated(phecode_icd9_map),]

icds <- unique(phecode_icd9_map$icd9)
icd9_phecode_mapping <- sapply(icds,function(cod){
  max(phecode_icd9_map$phecode[phecode_icd9_map$icd9==cod])
})
phecodes <- sort(unique(icd9_phecode_mapping))

icd9_phecode_mapping_1d <- sort(icd9_phecode_mapping)
for (i in 2:length(icd9_phecode_mapping_1d)){
  if (round(icd9_phecode_mapping_1d[i],digits=1) != icd9_phecode_mapping_1d[i] &
      round(icd9_phecode_mapping_1d[i],digits=1) == icd9_phecode_mapping_1d[i-1]){
    names(icd9_phecode_mapping_1d)[i] <- names(icd9_phecode_mapping_1d)[i-1]
    icd9_phecode_mapping_1d[i] <- icd9_phecode_mapping_1d[i-1]
  }
}
phecodes_1d <- sort(unique(icd9_phecode_mapping_1d))

save(icd9_phecode_mapping,file='ICD9_PheCode_mapping.RData')
save(icd9_phecode_mapping_1d,file='ICD9_PheCode_mapping_1decimal.RData')

icd_events <- read.csv('DIAGNOSES_ICD.csv')
icd_events$ICD9_CODE <- sapply(icd_events$ICD9_CODE,function(cod){
  if (substr(cod,1,1) == 'E'){
    if (nchar(cod) > 4){
      prefix <- substr(cod,1,4)
      suffix <- substr(cod,5,nchar(cod))
      paste0(prefix,'.',suffix)
    }
    else{
      cod
    }
  }
  else if (substr(cod,1,1) == 'V'){
    if (nchar(cod) > 3){
      prefix <- substr(cod,1,3)
      suffix <- substr(cod,4,nchar(cod))
      paste0(prefix,'.',suffix)
    }
    else{
      cod
    }
  }
  else{
    if (nchar(cod) > 3){
      prefix <- substr(cod,1,3)
      suffix <- substr(cod,4,nchar(cod))
      if (suffix == '0' | suffix == '00'){
        prefix
      }
      else{
        paste0(prefix,'.',suffix)
      }
    }
    else{
      cod
    }
  }
})

patIDs <- unique(icd_events$HADM_ID)
mimic_prior_phecodes <- t(sapply(1:length(patIDs),function(i){
  if (i %% 1000 == 0){
    print(i)
  }
  
  pat <- patIDs[i]
  row <- rep(0,length(phecodes))
  names(row) <- phecodes
  
  icds <- icd_events$ICD9_CODE[icd_events$HADM_ID==pat]
  cods <- as.character(icd9_phecode_mapping[icds])
  for (it in 1:2){
    if (anyNA(cods)){
      cods[is.na(cods)] <- as.character(icd9_phecode_mapping[substr(icds[is.na(cods)],1,nchar(icds[is.na(cods)])-1)])
    }
  }
  cods <- cods[!is.na(cods)]
  row[cods] <- 1
  
  names(row) <- NULL
  row
}))
rownames(mimic_prior_phecodes) <- patIDs
colnames(mimic_prior_phecodes) <- phecodes
save(mimic_prior_phecodes, file='MIMIC_PheCode_prior.RData')

mimic_prior_phecodes_1d <- t(sapply(1:length(patIDs),function(i){
  if (i %% 1000 == 0){
    print(i)
  }
  
  pat <- patIDs[i]
  row <- rep(0,length(phecodes_1d))
  names(row) <- phecodes_1d
  
  icds <- icd_events$ICD9_CODE[icd_events$HADM_ID==pat]
  cods <- as.character(icd9_phecode_mapping_1d[icds])
  for (it in 1:2){
    if (anyNA(cods)){
      cods[is.na(cods)] <- as.character(icd9_phecode_mapping_1d[substr(icds[is.na(cods)],1,nchar(icds[is.na(cods)])-1)])
    }
  }
  cods <- cods[!is.na(cods)]
  row[cods] <- 1
  
  names(row) <- NULL
  row
}))
rownames(mimic_prior_phecodes_1d) <- patIDs
colnames(mimic_prior_phecodes_1d) <- phecodes_1d
save(mimic_prior_phecodes_1d, file='MIMIC_PheCode1d_prior.RData')

