setwd('~/Documents/HMS/Biostatistics\ PhD/Research/mixEHR-sureLDA')
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

load('ICD9_PheCode_mapping_1decimal.RData')
# icd9_phecode_mapping_1d <- read.csv('mixmimic/ICD9_PheCode_mapping_1decimal.csv')

icd_events <- read.csv('DIAGNOSES_ICD.csv')
icd_events$ICD9_CODE <- as.character(icd_events$ICD9_CODE)
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

patIDs <- sort(unique(icd_events$SUBJECT_ID))
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

mimic_prior_phecodes_1d <- foreach(i=1:length(patIDs),.combine=rbind) %dopar% {
  if (i %% 1000 == 0){
    print(i)
  }
  
  pat <- patIDs[i]

  icds <- icd_events$ICD9_CODE[icd_events$SUBJECT_ID==pat]
  cods <- as.character(icd9_phecode_mapping_1d[icds])
  for (it in 1:2){
    if (anyNA(cods)){
      cods[is.na(cods)] <- as.character(icd9_phecode_mapping_1d[substr(icds[is.na(cods)],1,nchar(icds[is.na(cods)])-1)])
    }
  }
  cods <- cods[!is.na(cods)]
  if (length(cods)>0){
    codIdx <- sapply(cods,function(cod){which(phecodes_1d==cod)}) - 1
    if (length(codIdx)>0){
      cbind(pat,codIdx,1)
    }
    else{
      c(pat,length(phecodes_1d),1)
    }
  }
  else{
    c(pat,length(phecodes_1d),1)
  }
}

save(mimic_prior_phecodes_1d, file='MIMIC_PheCode1d_prior.RData')



## Interpret MIMIC results

featIDs <- get(load('ehrFeatId.RData'))
phecodes_1d <- read.csv('phecodes.csv')[,1]
phi <- read.csv("mimic_trainData_JCVB0_nmar_K1515_iter100_phi_normalized.csv")
cods <- c(411.2,585.3,295.1,296.1,428.3,428.4,250.1,250.2)
codIdx <- sapply(cods,function(cod){which(phecodes_1d==cod)})
topFeatures <- foreach (i=3:ncol(phi), .combine=cbind, .packages=c('foreach','doParallel')) %dopar% {
  print(i)
  topIdx = order(phi[,i],decreasing=TRUE)[1:50]
  foreach(j=1:50, .combine=rbind) %do% {
    typeId = phi[topIdx[j],1]
    pheId = phi[topIdx[j],2]
    c(featIDs[[typeId]][featIDs[[typeId]][,'pheId']==pheId, 'pheName'], phi[topIdx[j],i], names(featIDs)[typeId])
  }
}
topFeatures <- array(topFeatures,dim=c(50,3,1515))
topFeatures <- aperm(topFeatures,c(1,3,2))
rownames(topFeatures) <- 1:50
colnames(topFeatures) <- phecodes_1d
save(topFeatures,file='mixEHR_sureLDA_MIMIC_topFeatures.RData')

nameMap <- c('#0000FF','#00FF00','#FF0000','#FFA500','#800080','#FFC0CB')
names(nameMap) <- c('icd_cm','drg','icd_cpt','presc','lab','notes')
wordClouds <- lapply(codIdx,function(cod){
  names <- sub(',', '', sapply(topFeatures[,cod,1],function(nam){
    if (strsplit(nam,'_')[[1]][1]=='MS'){strsplit(nam,'_')[[1]][3]}
    else{strsplit(nam,'_')[[1]][2]}
  }))
  scores <- as.numeric(topFeatures[,cod,2])
  scores <- round(100*scores/max(scores))
  type <- nameMap[topFeatures[,cod,3]]
  data.frame(scores,names,type)
})
for (i in 1:length(cods)){
  write.csv(wordClouds[[i]][!is.na(wordClouds[[i]][,2]),],file=paste0('wordCloud_',cods[i],'.csv'),row.names=FALSE)
}
  

# Heat Map
library(ggplot2)
library(xtable)
library(plyr)
library(scales)
library(gridExtra)
library(tidytext)
library(ComplexHeatmap)
library(circlize)
library(reshape2)

phi <- as.matrix(phi[phi[,1]==2,-c(1:2)])
phiCorr <- cor(t(phi),method='pearson')
phiCorr <- as.matrix(read.csv('phiCorr.csv'))
icd9s <- as.character(read.csv('mimic_icd9s.csv')[,1])
rownames(phiCorr) <- colnames(phiCorr) <- icd9s

icd9Codes <- as.character(read.csv('mimic_icd9_codes.csv')[,2])
lookup <- data.frame('category'=c('Infectious','Oncologic','Endocrine','Hematologic','Psychologic',
                                  'Nervous','Cardiovascular','Respiratory','Gastrointestinal','Genitourinary',
                                  'Pregnancy+Childbirth','Dermatologic','Musculoskeletal','Congenital',
                                  'Perinatal','Symptoms+Signs','Injury+Poisoning','External+Supplement'),
                     'lower'=c(1,140,240,280,290,320,390,460,520,580,630,680,710,740,760,780,800,NA),
                     'upper'=c(139,239,279,289,319,389,459,519,579,629,679,709,739,759,779,799,999,NA))
category <- sapply(icd9Codes,function(cod){
  if (substr(cod,1,1)=='E' | substr(cod,1,1)=='V'){18}
  else{
    codint <- as.integer(substr(cod,1,3))
    which(lookup$lower <= codint & lookup$upper >= codint)
  }
})
colors <- rainbow(18)

hmp <- heatmap(phiCorr, col=cm.colors(256), RowSideColors=colors[category],
               ColSideColors=colors[category], symm=TRUE, margins=c(5,5), scale=NULL)

tiff('mimic_mixehr_surelda_heatmap.tiff')
hmp
dev.off()



## Preprocess data
library(foreach)
library(doParallel)

icds <- read.csv('/home/yueli/data/services.csv')
icds_list <- as.character(read.csv('/home/yueli/data/icd9_qc.csv')[-c(1:2),1])
actes_list <- read.csv('/home/yueli/data/actes.csv')
save(icds_list,file='icd9_list.RData')
# icds_sorted <- sort(unique(as.character(icds$icd)))
# icds_sorted = icds_sorted[-1]
# save(icds_sorted,file='ICDs_unique_sorted.RData')
# load('ICDs_unique_sorted.RData')
icds$codIdx <- sapply(1:nrow(icds),function(i){
  if (i%%10000==0){print(paste0(i,'/',nrow(icds)))}
  match <- which(icds_list==icds$icd[i])
  if (length(match)==0){NA}
  else{match}
})
icds_mixehr_surelda <- icds[,c('id','date','codIdx')]
icds_mixehr_surelda <- icds_mixehr_surelda[!is.na(icds_mixehr_surelda[,'codIdx']),]
colnames(icds_mixehr_surelda) <- c('patId','date','pheId')
save(icds_mixehr_surelda,file='icds_mixehr_surelda.RData')
metadata <- read.csv('icd9_qc.csv')[3:7237,1:3]
metadata <- cbind(seq(7235),metadata)
colnames(metadata)[1] <- 'pheId'
save(metadata,file='metadata.RData')

drugs <- read.csv('drugs.csv')
cod_list <- sort(as.numeric(unique(drugs[,4])))
save(cod_list, file='drug_code_list.RData')
drugs$codIdx <- sapply(1:nrow(drugs),function(i){
  if (i%%10000==0){print(paste0(i,'/',nrow(drugs)))}
  match <- which(cod_list==drugs$din[i])
  if (length(match)==0){NA}
  else{match}
})
drugs_mixehr_surelda <- drugs[,c('id','date','codIdx')]
colnames(drugs_mixehr_surelda) <- c('patId','date','pheId')
save(drugs_mixehr_surelda,file='drugs_mixehr_surelda.RData')
metadata2 <- read.csv('drug_product.csv')
metadata2 <- metadata2[,c('drug_code','drug_identification_number','brand_name')]
metadata2$drug_identification_number <- as.numeric(as.character(metadata2$drug_identification_number))
metadata2$pheId <- sapply(metadata2$drug_identification_number,function(id){
  if (id %in% cod_list){which(cod_list==id)}
  else{NA}
})
metadata2 <- metadata2[!is.na(metadata2$pheId),c('pheId','drug_code','drug_identification_number','brand_name')]
metadata2 <- metadata2[order(metadata2$pheId),]


drugs <- read.csv('/home/yueli/data/drugs.csv')
drugs$date <- as.character(drugs$date)
drug_product <- read.csv('/home/yueli/data/drug_product.csv')
drug_product$drug_identification_number <- as.integer(as.character(drug_product$drug_identification_number))
drug_ingredients <- read.csv('/home/yueli/data/active_ingredients.csv')
n <- nrow(drugs)

drugs <- merge(drugs,drug_product[,c('drug_code','drug_identification_number')],by.x='din',by.y='drug_identification_number')
drugs <- merge(drugs,drug_ingredients[,c('drug_code','active_ingredient_code','ingredient')],by='drug_code')
drugs <- drugs[,c('id','date','active_ingredient_code','ingredient')]
save(drugs,file='drugs_byIngredient.RData')
colnames(drugs) <- c('id','date','pheId','name')

cod_list <- sort(as.numeric(unique(drugs$pheId)))
save(cod_list, file='drug_ingredient_code_list.RData')

drugs_mixehr_surelda <- drugs[,c('patId','date','pheId')]
save(drugs_mixehr_surelda,file='drugs_mixehr_surelda.RData')
metadata2 <- drugs[,c('pheId','name')]
metadata2 <- metadata2[order(metadata2$date),]
metadata2 <- metadata2[order(metadata2$pheId),]
metadata2 <- metadata2[!duplicated(metadata2),]
save(metadata2,file='metadata_drugs.RData')
metadata$Drug_byIngredient <- metadata2

# Reconcile ICD, drug files, add ACT codes
load('metadata.RData')
actes_metadata <- actes_list[,2:3]
colnames(actes_metadata) <- c('pheId','description')
metadata <- list('ICD9'=metadata, 'Drug'=metadata2, 'Act'=actes_metadata)
save(metadata,file='metadata.RData')

load('icds_mixehr_surelda.RData')
icds_mixehr_surelda$typeId <- 1
drugs_mixehr_surelda$typeId <- 2
icds_mixehr_surelda$state <- drugs_mixehr_surelda$state <- 1
combined_data <- rbind(icds_mixehr_surelda, drugs_mixehr_surelda)
pheIds <- sapply(1:nrow(combined_data),function(i){
  if (i %% 100000 == 0){print(paste0(i,'/',nrow(combined_data)))}
  combined_data$pheId[[i]][1]
})
actes_mixehr_surelda <- data.frame('patId'=icds$id,'date'=icds$date,'pheId'=icds$act_code,
                                   'typeId'=3,'state'=1)
combined_data <- rbind(combined_data,actes_mixehr_surelda)
combined_data <- combined_data[order(combined_data$pheId),]
combined_data <- combined_data[order(combined_data$typeId),]
combined_data <- combined_data[order(combined_data$patId),]
save(combined_data,file='combined_data.RData')

## Using drugs by ingredient in lieu of drugs by DIN code
load('drugs_mixehr_surelda.RData')
drugs_mixehr_surelda$typeId <- 4; drugs_mixehr_surelda$state <- 1
combined_data <- combined_data[combined_data$typeId!=2,]
combined_data <- rbind(combined_data,drugs_mixehr_surelda)
combined_data <- combined_data[order(combined_data$pheId),]
combined_data <- combined_data[order(combined_data$typeId),]
combined_data <- combined_data[order(combined_data$patId),]
save(combined_data,file='combined_data_drugsByIngredient.RData')
##

load('combined_data_drugsByIngredient.RData')
firstHalf <- as.Date(combined_data$date) - as.Date('2007-01-01') < 0
combined_data_firstHalf <- combined_data[firstHalf,]
combined_data_firstHalf <- combined_data_firstHalf[order(combined_data_firstHalf$pheId),]
combined_data_firstHalf <- combined_data_firstHalf[order(combined_data_firstHalf$typeId),]
combined_data_firstHalf <- combined_data_firstHalf[order(combined_data_firstHalf$patId),]
# save(combined_data_firstHalf,file='combined_data_firstHalf.RData')
save(combined_data_firstHalf,file='combined_data_firstHalf_drugsByIngredient.RData')

duplicates <- duplicated(combined_data[,c('patId','typeId','pheId')])
uniqueCases <- combined_data[!duplicates,c('patId','typeId','pheId','state')]

duplicates <- duplicated(combined_data_firstHalf[,c('patId','typeId','pheId')])
uniqueCases_firstHalf <- combined_data_firstHalf[!duplicates,c('patId','typeId','pheId','state')]
uniqueCases_firstHalf <- uniqueCases_firstHalf[order(uniqueCases_firstHalf$pheId),]
uniqueCases_firstHalf <- uniqueCases_firstHalf[order(uniqueCases_firstHalf$typeId),]
uniqueCases_firstHalf <- uniqueCases_firstHalf[order(uniqueCases_firstHalf$patId),]


library(Rcpp)
library(RcppArmadillo)
sourceCpp('create_datamat_david.cpp')

freqs <- findCounts(as.matrix(uniqueCases[,1:3]),as.matrix(combined_data[,c(1,4,3)]))
uniqueCases$freq <- freqs
save(uniqueCases,file='combinedData_mixEHRready_drugsByIngredient.RData')
write.table(uniqueCases,file='combinedData_drugsByIngredient.txt',sep=' ',row.names=FALSE,col.names=FALSE)

freqs <- findCounts(as.matrix(uniqueCases_firstHalf[,1:3]),as.matrix(combined_data_firstHalf[,c(1,4,3)]))
uniqueCases_firstHalf$freq <- freqs
save(uniqueCases_firstHalf,file='combinedData_firstHalf_mixEHRready_drugsByIngredient.RData')
write.table(uniqueCases_firstHalf,file='combinedData_firstHalf_drugsByIngredient.txt',sep=' ',row.names=FALSE,col.names=FALSE)

combined_metadata <- uniqueCases[,c('typeId','pheId')]
combined_metadata <- combined_metadata[order(combined_metadata$pheId),]
combined_metadata <- combined_metadata[order(combined_metadata$typeId),]
combined_metadata <- combined_metadata[!duplicated(combined_metadata),]
combined_metadata$stateCnt <- 1
write.table(combined_metadata,file='combinedMetaData_withactes_drugsByIngredient.txt',sep=' ',row.names=FALSE,col.names=FALSE)

combined_metadata_firstHalf <- uniqueCases_firstHalf[,c('typeId','pheId')]
combined_metadata_firstHalf <- combined_metadata_firstHalf[order(combined_metadata_firstHalf$pheId),]
combined_metadata_firstHalf <- combined_metadata_firstHalf[order(combined_metadata_firstHalf$typeId),]
combined_metadata_firstHalf <- combined_metadata_firstHalf[!duplicated(combined_metadata_firstHalf),]
combined_metadata_firstHalf$stateCnt <- 1
write.table(combined_metadata_firstHalf,file='combinedMetaData_firstHalf_withactes_drugsByIngredient.txt',sep=' ',row.names=FALSE,col.names=FALSE)

## MAKE SURE TO REMOVE PATIENTS FROM DATA WITH ALL PRIORS=0; WILL RESULT IN ERRORS


# Filter out rare codes
uniqueCodes <- uniqueCases[!duplicated(uniqueCases[,c('typeId','pheId')]), c('typeId','pheId')]
nPats <- length(unique(uniqueCases$patId))
codeCounts <- table(uniqueCases[,c('typeId','pheId')])
codeFreqs <- codeCounts/nPats
save(codeFreqs,file='MTL_codeFrequencies.RData')

stopWords <- as.data.frame(t(sapply(which(codeFreqs>=0.25),function(x){
  i <- (x-1)%%3 + 1
  j <- ceiling(x/3)
  as.integer(c(rownames(codeFreqs)[i], colnames(codeFreqs)[j]))
})))
colnames(stopWords) <- c('typeId','pheId')
for (i in 1:nrow(stopWords)){
  print(metadata[[stopWords$typeId[i]]][metadata[[stopWords$typeId[i]]]$pheId==stopWords$pheId[i],])
}

uniqueCases_stopIdx <- which((1000000*uniqueCases$typeId + uniqueCases$pheId) %in% (1000000*stopWords$typeId + stopWords$pheId))
uniqueCases <- uniqueCases[-uniqueCases_stopIdx,]
keep <- which(uniqueCases$patId %in% map_prior_compressed[,1])
uniqueCases <- uniqueCases[keep,]
save(uniqueCases,file='combinedData_mixEHRready_withactes_noStopWords.RData')
write.table(uniqueCases,file='combinedData_withactes_noStopWords.txt',sep=' ',row.names=FALSE,col.names=FALSE)

uniqueCases_firstHalf_stopIdx <- which((1000000*uniqueCases_firstHalf$typeId + uniqueCases_firstHalf$pheId) %in% (1000000*stopWords$typeId + stopWords$pheId))
uniqueCases_firstHalf <- uniqueCases_firstHalf[-uniqueCases_firstHalf_stopIdx,]
keep <- which(uniqueCases_firstHalf$patId %in% map_prior_firstHalf_compressed[,1])
uniqueCases_firstHalf <- uniqueCases_firstHalf[keep,]
save(uniqueCases_firstHalf,file='combinedData_firstHalf_mixEHRready_withactes_noStopWords.RData')
write.table(uniqueCases_firstHalf,file='combinedData_firstHalf_withactes_noStopWords.txt',sep=' ',row.names=FALSE,col.names=FALSE)


combined_metadata <- uniqueCases[,c('typeId','pheId')]
combined_metadata <- combined_metadata[order(combined_metadata$pheId),]
combined_metadata <- combined_metadata[order(combined_metadata$typeId),]
combined_metadata <- combined_metadata[!duplicated(combined_metadata),]
combined_metadata$stateCnt <- 1
write.table(combined_metadata,file='combinedMetaData_withactes_noStopWords.txt',sep=' ',row.names=FALSE,col.names=FALSE)

combined_metadata_firstHalf <- uniqueCases_firstHalf[,c('typeId','pheId')]
combined_metadata_firstHalf <- combined_metadata_firstHalf[order(combined_metadata_firstHalf$pheId),]
combined_metadata_firstHalf <- combined_metadata_firstHalf[order(combined_metadata_firstHalf$typeId),]
combined_metadata_firstHalf <- combined_metadata_firstHalf[!duplicated(combined_metadata_firstHalf),]
combined_metadata_firstHalf$stateCnt <- 1
write.table(combined_metadata_firstHalf,file='combinedMetaData_firstHalf_withactes_noStopWords.txt',sep=' ',row.names=FALSE,col.names=FALSE)



## Generate MAP prior
library(foreach)
library(doParallel)
library(MAP)

icds_list <- get(load('icd9_list.RData'))
phecode_icd_mappings <- read.csv('phecode_icd9_map_unrolled.csv',row.names=NULL)
phecode_icd_mappings$phecode <- round(phecode_icd_mappings$phecode,digits=1)
phecode_icd_mappings <- phecode_icd_mappings[!duplicated(phecode_icd_mappings),]
phecode_icd_mappings$icd9 <- as.character(phecode_icd_mappings$icd9)

for (i in 1:(nrow(phecode_icd_mappings)-1)){
  if (i %% 1000 == 0){print(i)}
  if (sum(phecode_icd_mappings$icd9 == phecode_icd_mappings$icd9[i]) > 1 &&
      round(phecode_icd_mappings$phecode[i])==phecode_icd_mappings$phecode[i]){
    phecode_icd_mappings$phecode[i] <- NA
  }
}
phecode_icd_mappings <- phecode_icd_mappings[!is.na(phecode_icd_mappings$phecode),]
write.csv(phecode_icd_mappings,'phecode_icd9_mappings_fixed102920.csv',row.names=FALSE)

phecodes <- sort(unique(phecode_icd_mappings$phecode))
ints <- 1:11031
Es <- which(sapply(phecode_icd_mappings$icd9,function(cod){substr(cod,1,1)=='E'}))
Vs <- which(sapply(phecode_icd_mappings$icd9,function(cod){substr(cod,1,1)=='V'}))

pheCode_dat <- uniqueCases[uniqueCases$typeId==1,]
pheCode_dat$ICD <- icds_list[pheCode_dat$pheId]
pheCode_dat$phecode <- sapply(1:nrow(pheCode_dat),function(i){
  if (i %% 10000 == 0){print(paste(i,'/',nrow(pheCode_dat)))}
  icd <- pheCode_dat$ICD[i]
  if (substr(icd,1,1)=='E'){
    'E'
  }
  else if (substr(icd,1,1)=='V'){
    substr(icd,1,2)
  }
  else{
    code <- as.character(as.numeric(paste0(substr(icd,1,3),'.',substr(icd,4,4))))
    phecode <- phecode_icd_mappings$phecode[phecode_icd_mappings$icd9==code]
    if (length(phecode) == 0){
      code <- as.character(floor(as.numeric(code)))
      phecode <- phecode_icd_mappings$phecode[phecode_icd_mappings$icd9==code]
    }
    phecode
  }
})

all_phecode <- sort(unique(unlist(pheCode_dat$phecode)))
all_pats <- unique(pheCode_dat$patId)
pat_HU <- pheCode_dat[,c('patId','freq')]
HU <- pat_HU %>%
  group_by(patId) %>%
  transmute(Total=sum(freq))
HU <- HU[!duplicated(HU$patId),]
save(HU,file='HU_byPat.RData')

lens <- sapply(pheCode_dat$phecode,length)
pheCode_dat <- pheCode_dat[lens!=0,c('patId','phecode','freq')]
lens <- sapply(pheCode_dat$phecode,length)
pheCode_dat_1cod <- pheCode_dat[lens==1,]
pheCode_dat_2cod <- pheCode_dat[lens==2,]
pheCode_dat_2cod <- data.frame('patId'=rep(pheCode_dat_2cod$patId,each=2),
                               'phecode'=unlist(pheCode_dat_2cod$phecode),
                               'freq'=rep(pheCode_dat_2cod$freq,each=2))
pheCode_dat <- rbind(pheCode_dat_1cod,pheCode_dat_2cod)
pheCode_dat$phecode <- unlist(pheCode_dat$phecode)
save(pheCode_dat,file='pheCode_dat.RData')

library(dplyr)

pheCode_dat_summed <- pheCode_dat %>%
  group_by(patId,phecode) %>%
  transmute(Total=sum(freq))
pheCode_dat_summed <- pheCode_dat_summed[!duplicated(pheCode_dat_summed[,1:2]),]
pheCode_dat_summed <- as.data.frame(pheCode_dat_summed)
save(pheCode_dat_summed,file='pheCode_dat_summed.RData')

library(doParallel)

logfile <- "map_prior_predictor.txt"
writeLines(c(""), file(logfile,'w'))
clust <- makeCluster(20, outfile=logfile)
registerDoParallel(clust)

map_prior <- foreach(i=1:length(all_phecode), .combine=cbind, .packages=c('Matrix','MAP')) %dopar% {
  print(paste('Phenotype',i,'of',length(all_phecode)))
  tryCatch({
    phe <- all_phecode[i]
    matches <- which(pheCode_dat_summed$phecode == phe)
    pats <- as.character(pheCode_dat_summed$patId[matches])
    note <- Matrix(HU[pats,'Total'], sparse=TRUE)
    mat <- Matrix(data=pheCode_dat_summed$Total[matches], sparse=TRUE)
    colnames(mat) <- 'ICD'
    res <- MAP(mat=mat, note=note)
    out <- rep(0,length(all_pats)); names(out) <- all_pats
    out[pats] <- res$scores
    out
  }, error=function(e){
    rep(NA,length(all_pats))
  })
}
save(map_prior,file='map_prior.RData')

nas = which(is.na(map_prior[1,]))
for (i in nas){
  print(i)
  phe <- all_phecode[i]
  matches <- which(pheCode_dat_summed$phecode == phe)
  pats <- as.character(pheCode_dat_summed$patId[matches])
  out <- rep(0,length(all_pats)); names(out) <- all_pats
  out[pats] <- 1
  map_prior[,i] <- out
}

map_prior_compressed <- foreach(i=1:nrow(map_prior), .combine=rbind) %do% {
  if (i %% 10000 == 0){print(paste(i,'/',length(all_pats)))}
  nonZeros <- which(map_prior[i,]>0)
  cbind(rep(all_pats[i],length(nonZeros)), all_phecode[nonZeros], map_prior[i,nonZeros])
}
observed_phecodes <- sort(unique(map_prior_compressed[,2]))
observed_phecodes <- intersect(all_phecode,observed_phecodes)
pheId_metadata <- data.frame('pheId'=seq(length(observed_phecodes))-1, 'phecodes'=observed_phecodes)
save(pheId_metadata,file='pheId_metadata.RData')

rownames(pheId_metadata) <- observed_phecodes
pheIds <- pheId_metadata[as.character(map_prior_compressed[,2]),'pheId']
map_prior_compressed[,2] <- pheIds


save(map_prior_compressed, file='map_prior_compressed.RData')
write.table(map_prior_compressed, file='map_prior_compressed.txt', sep=' ', row.names=FALSE, col.names=FALSE)


# MAP prior first half
uniqueCases_firstHalf <- uniqueCases_firstHalf[uniqueCases_firstHalf$freq!=0,]
pheCode_dat_firstHalf <- uniqueCases_firstHalf[uniqueCases_firstHalf$typeId==1,]
pheCode_dat_firstHalf$ICD <- icds_list[pheCode_dat_firstHalf$pheId]
pheCode_dat_firstHalf$phecode <- sapply(1:nrow(pheCode_dat_firstHalf),function(i){
  if (i %% 10000 == 0){print(paste(i,'/',nrow(pheCode_dat_firstHalf)))}
  icd <- pheCode_dat_firstHalf$ICD[i]
  if (substr(icd,1,1)=='E'){
    'E'
    
  }
  else if (substr(icd,1,1)=='V'){
    substr(icd,1,2)
  }
  else{
    code <- as.character(as.numeric(paste0(substr(icd,1,3),'.',substr(icd,4,4))))
    phecode <- phecode_icd_mappings$phecode[phecode_icd_mappings$icd9==code]
    if (length(phecode) == 0){
      code <- as.character(floor(as.numeric(code)))
      phecode <- phecode_icd_mappings$phecode[phecode_icd_mappings$icd9==code]
    }
    phecode
  }
})

all_phecode_firstHalf <- sort(unique(unlist(pheCode_dat_firstHalf$phecode)))
all_pats_firstHalf <- unique(pheCode_dat_firstHalf$patId)
pat_HU_firstHalf <- pheCode_dat_firstHalf[,c('patId','freq')]
HU_firstHalf <- pat_HU_firstHalf %>%
  group_by(patId) %>%
  transmute(Total=sum(freq))
HU_firstHalf <- as.matrix(HU_firstHalf[!duplicated(HU_firstHalf$patId),])
rownames(HU_firstHalf) <- HU_firstHalf[,'patId']
save(HU_firstHalf,file='HU_firstHalf_byPat.RData')

lens_firstHalf <- sapply(pheCode_dat_firstHalf$phecode,length)
pheCode_dat_firstHalf <- pheCode_dat_firstHalf[lens_firstHalf!=0,c('patId','phecode','freq')]
lens_firstHalf <- sapply(pheCode_dat_firstHalf$phecode,length)
pheCode_dat_firstHalf_1cod <- pheCode_dat_firstHalf[lens_firstHalf==1,]
pheCode_dat_firstHalf_1cod$phecode <- unlist(pheCode_dat_firstHalf_1cod$phecode)
pheCode_dat_firstHalf_2cod <- pheCode_dat_firstHalf[lens_firstHalf==2,]
pheCode_dat_firstHalf_2cod <- data.frame('patId'=rep(pheCode_dat_firstHalf_2cod$patId,each=2),
                               'phecode'=unlist(pheCode_dat_firstHalf_2cod$phecode),
                               'freq'=rep(pheCode_dat_firstHalf_2cod$freq,each=2))
pheCode_dat_firstHalf <- rbind(pheCode_dat_firstHalf_1cod,pheCode_dat_firstHalf_2cod)
pheCode_dat_firstHalf$phecode <- unlist(pheCode_dat_firstHalf$phecode)
save(pheCode_dat_firstHalf,file='pheCode_dat_firstHalf.RData')

pheCode_dat_firstHalf_summed <- pheCode_dat_firstHalf %>%
  group_by(patId,phecode) %>%
  transmute(Total=sum(freq))
pheCode_dat_firstHalf_summed <- pheCode_dat_firstHalf_summed[!duplicated(pheCode_dat_firstHalf_summed[,1:2]),]
pheCode_dat_firstHalf_summed <- as.data.frame(pheCode_dat_firstHalf_summed)
save(pheCode_dat_firstHalf_summed,file='pheCode_dat_firstHalf_summed.RData')


load('HU_firstHalf_byPat.RData')
load('pheCode_dat_firstHalf_summed.RData')
all_phecode_firstHalf <- sort(unique(unlist(pheCode_dat_firstHalf_summed$phecode)))
all_pats_firstHalf <- unique(pheCode_dat_firstHalf_summed$patId)
pheCode_dat_firstHalf_uncompressed <- foreach(i=1:length(all_phecode_firstHalf), .combine=cbind, .packages=c('Matrix','MAP')) %do% {
  print(paste('Phenotype',i,'of',length(all_phecode_firstHalf)))
  phe <- all_phecode_firstHalf[i]
  matches <- which(pheCode_dat_firstHalf_summed$phecode == phe)
  pats <- as.character(pheCode_dat_firstHalf_summed$patId[matches])
  out <- rep(0,length(all_pats_firstHalf)); names(out) <- all_pats_firstHalf
  out[pats] <- pheCode_dat_firstHalf_summed$Total[matches]
  out
}

library(dplyr)
library(doParallel)

stopCluster(clust)
logfile <- "map_prior_predictor.txt"
writeLines(c(""), file(logfile,'w'))
clust <- makeCluster(20, outfile=logfile)
registerDoParallel(clust)

map_prior_firstHalf <- foreach(i=1:length(all_phecode_firstHalf), .combine=cbind, .packages=c('Matrix','MAP')) %do% {
  print(paste('Phenotype',i,'of',length(all_phecode_firstHalf)))
  tryCatch({
    phe <- all_phecode_firstHalf[i]
    matches <- which(pheCode_dat_firstHalf_summed$phecode == phe)
    pats <- as.character(pheCode_dat_firstHalf_summed$patId[matches])
    note <- Matrix(HU_firstHalf[pats,'Total'], sparse=TRUE)
    mat <- Matrix(data=pheCode_dat_firstHalf_summed$Total[matches], sparse=TRUE)
    colnames(mat) <- 'ICD'
    res <- MAP(mat=mat, note=note)
    out <- rep(0,length(all_pats_firstHalf)); names(out) <- all_pats_firstHalf
    out[pats] <- res$scores
    out
  }, error=function(e){
    rep(NA,length(all_pats_firstHalf))
  })
}
save(map_prior_firstHalf,file='map_prior_firstHalf.RData')

nas = which(is.na(map_prior_firstHalf[1,]))
for (i in nas){
  print(i)
  phe <- all_phecode_firstHalf[i]
  matches <- which(pheCode_dat_firstHalf_summed$phecode == phe)
  pats <- as.character(pheCode_dat_firstHalf_summed$patId[matches])
  out <- rep(0,length(all_pats_firstHalf)); names(out) <- all_pats_firstHalf
  out[pats] <- 1
  map_prior_firstHalf[,i] <- out
}

map_prior_firstHalf_compressed <- foreach(i=1:length(all_pats_firstHalf), .combine=rbind) %do% {
  if (i %% 10000 == 0){print(paste(i,'/',length(all_pats_firstHalf)))}
  nonZeros <- which(map_prior_firstHalf[i,]>0)
  cbind(rep(all_pats_firstHalf[i],length(nonZeros)), all_phecode_firstHalf[nonZeros], map_prior_firstHalf[i,nonZeros])
}

observed_phecodes <- sort(unique(map_prior_firstHalf_compressed[,2]))
observed_phecodes <- intersect(all_phecode,observed_phecodes)
pheId_metadata_firstHalf <- data.frame('pheId'=seq(length(observed_phecodes))-1, 'phecodes'=observed_phecodes)
save(pheId_metadata_firstHalf,file='pheId_metadata_firstHalf.RData')

rownames(pheId_metadata_firstHalf) <- observed_phecodes
pheIds <- pheId_metadata_firstHalf[as.character(map_prior_firstHalf_compressed[,2]),'pheId']
map_prior_firstHalf_compressed[,2] <- pheIds

save(map_prior_firstHalf_compressed, file='map_prior_firstHalf_compressed.RData')
write.table(map_prior_firstHalf_compressed, file='map_prior_firstHalf_compressed.txt', sep=' ', row.names=FALSE, col.names=FALSE)



## Interpret David's results
load('metadata.RData')
metadata[[2]] <- metadata[[2]][!duplicated(metadata[[2]][,1]),]
metadata[[3]] <- metadata[[3]][!duplicated(metadata[[3]][,1]),]
metadata[[4]] <- metadata[[4]][!duplicated(metadata[[4]][,1]),]

pheId_metadata_1dec <- get(load('pheId_metadata.RData'))
pheId_metadata_0dec <- get(load('pheId_metadata_integer.RData'))
phi_1dec <- read.csv('combinedData_withactes_noStopWords_JCVB0_nmar_K1327_iter200_phi_normalized.csv', header=FALSE)
phi_0dec <- read.csv('combinedData_sample_nonZero_JCVB0_nmar_K522_iter141_phi_normalized.csv', header=FALSE)

cod_1dec <- c(428.2,428.3,496.1,496.2,496.3,250.1,250.2,295.1,296.1,296.2,411.1,411.2)
cod_0dec <- c(428,496,585,571,250,295,411)

findTopFeatures <- function(phi,metadata,cods){
  codIdx <- c(sapply(cods,function(cod){
    if (cod %in% metadata$phecodes){
      which(metadata$phecodes==cod)
    }
    else{NA}
  }))
  
  topFeatures <- foreach (i=codIdx, .combine=cbind, .packages=c('foreach','doParallel')) %do% {
    if (is.na(i)){rep(NA,50)}
    else{
      topIdx = order(phi[,i],decreasing=TRUE)[1:50]
      foreach(j=1:50, .combine=rbind) %do% {
        typeId = phi[topIdx[j],1]
        pheId = phi[topIdx[j],2]
        if (typeId==1){
          c(as.character(metadata[[1]][metadata[[1]][,'pheId']==pheId, 'desc_en']), phi[topIdx[j],i], 'ICD')
        }
        else if (typeId==2){
          c(as.character(metadata[[2]][metadata[[2]][,'pheId']==pheId, 'brand_name']), phi[topIdx[j],i], 'Rx')
        }
        else{
          if (pheId %in% metadata[[3]][,'pheId']){c(as.character(metadata[[3]][metadata[[3]][,'pheId']==pheId, 'description']), phi[topIdx[j],i], 'ACT')}
          else{c(NA, phi[topIdx[j],i], 'ACT')}
        }
      }
    }
  }
}


topFeatures <- foreach (i=3:ncol(phi), .combine=cbind, .packages=c('foreach','doParallel')) %dopar% {
  print(i)
  topIdx = order(phi[,i],decreasing=TRUE)[1:50]
  foreach(j=1:50, .combine=rbind) %do% {
    typeId = phi[topIdx[j],1]
    pheId = phi[topIdx[j],2]
    if (typeId==1){
      c(as.character(metadata[[1]][metadata[[1]][,'pheId']==pheId, 'desc_en']), phi[topIdx[j],i], 'ICD')
    }
    else if (typeId==2){
      c(as.character(metadata[[2]][metadata[[2]][,'pheId']==pheId, 'brand_name']), phi[topIdx[j],i], 'Rx')
    }
    else{
      if (pheId %in% metadata[[3]][,'pheId']){c(as.character(metadata[[3]][metadata[[3]][,'pheId']==pheId, 'description']), phi[topIdx[j],i], 'ACT')}
      else{c(NA, phi[topIdx[j],i], 'ACT')}
    }
  }
}
topFeatures <- array(topFeatures,dim=c(50,3,2653))
topFeatures <- aperm(topFeatures,c(1,3,2))
rownames(topFeatures) <- 1:50
colnames(topFeatures) <- c(as.character(pheId_metadata$phecodes),'E')
save(topFeatures,file='mixEHR_sureLDA_MTL_topFeatures_SES.RData')

nameMap <- c('#0000FF','#00FF00','#FF0000')
names(nameMap) <- c('ICD','Rx','ACT')
wordClouds <- lapply(codIdx,function(cod){
  names <- sapply(strsplit(topFeatures[,cod,1],','),function(splt){splt[1]})
  scores <- as.numeric(topFeatures[,cod,2])
  scores <- as.integer(round(100*scores/max(scores,na.rm=TRUE)))
  type <- nameMap[topFeatures[,cod,3]]
  data.frame(scores,names,type)
})
for (i in 1:length(cods)){
  write.csv(wordClouds[[(2*i-1)]][wordClouds[[(2*i-1)]][,1]>0,],file=paste0('wordCloud_',cods[i],'_MTL_SEShigh.csv'),row.names=FALSE)
  write.csv(wordClouds[[(2*i)]][wordClouds[[(2*i)]][,1]>0,],file=paste0('wordCloud_',cods[i],'_MTL_SESlow.csv'),row.names=FALSE)
}



# LASSO prediction of 1) 1-year readmission,  2) 1-year mortality, 3) # of drugs,
# and 4) # of ED admissions in followup year among patients with at least 1 year of follow-up data

load('combined_data.RData')
load("map_prior_firstHalf_compressed.RData")
admissions <- read.csv('/home/yueli/data/hospital.csv')
patients <- read.csv('/home/yueli/data/patients.csv')

load('EDadmits_0708_table.RData')
keepIdx <- c(which(patients$month_of_death==''),
             which(patients$month_of_death!='')[as.Date(paste0(patients$month_of_death[patients$month_of_death!=''],'-15')) - as.Date('2007-01-01') >= 0])
alivePats <- patients$id[keepIdx]
patients_alive <- patients[patients$id %in% intersect(alivePats,keepPats0715),]
patients_alive <- patients_alive[order(patients_alive$id),]
patient_table <- data.frame('patId'=patients_alive$id,
                            'dead'=sapply(patients_alive$month_of_death,function(mod){
                              (mod!='') && as.Date(paste0(mod,'-15'))-as.Date('2007-01-01')>=0 &&
                              as.Date(paste0(mod,'-15'))-as.Date('2008-01-01')<0}))
save(patient_table,file='patient_table.RData')

dates = as.Date(combined_data$date)
keepPats0007 <- unique(map_prior_firstHalf_compressed[,1])
keepPats0715 <- unique(combined_data$patId[dates-as.Date('2007-01-01')>=0])
keepPats0815 <- unique(combined_data$patId[dates-as.Date('2008-01-01')>=0])
targetIdx <- (dates - as.Date('2007-01-01') >= 0 &
                dates - as.Date('2008-01-01') < 0)
save(keepPats0715,file='keepPats0715.RData')
save(keepPats0815,file='keepPats0815.RData')
save(targetIdx,file='targetIdx.RData')

load('keepPats0715.RData')
load('keepPats0815.RData')
load('targetIdx.RData')
combined_data_0708 <- combined_data[targetIdx,]
combined_data_0708 <- combined_data_0708[combined_data_0708$patId %in% intersect(intersect(keepPats0007,keepPats0815),alivePats),]
save(combined_data_0708,file='combined_data_0708.RData')

targetAdmissions <- (as.Date(admissions$admit) - as.Date('2007-01-01') >= 0 &
                       as.Date(admissions$admit) - as.Date('2008-01-01') < 0)
admissions_0708 <- admissions[targetAdmissions,]
admissions_0708 <- admissions_0708[admissions_0708$id %in% intersect(keepPats0007,keepPats0815),]

admissions_0708_table <- as.data.frame(rbind(cbind(unique(admissions_0708$id),TRUE),
                                             cbind(setdiff(intersect(keepPats0007,keepPats0815),unique(admissions_0708$id)),FALSE)))
colnames(admissions_0708_table) <- c('patId','admission')
admissions_0708_table <- admissions_0708_table[order(admissions_0708_table$patId),]
save(admissions_0708_table,file='admissions_0708_table.RData')

load('admissions_0708_table.RData')
EDadmits_0708_table <- as.data.frame(table(admissions_0708$id[admissions_0708$ed_admit != '']))
colnames(EDadmits_0708_table) <- c('patId','Count'); EDadmits_0708_table$patId <- as.integer(as.character(EDadmits_0708_table$patId))
noEDpats <- as.data.frame(cbind(setdiff(intersect(keepPats0007,keepPats0815),unique(EDadmits_0708_table$patId)),0))
colnames(noEDpats) <- c('patId','Count'); noEDpats$patId <- as.integer(noEDpats$patId)
EDadmits_0708_table <- rbind(EDadmits_0708_table,noEDpats)
EDadmits_0708_table <- EDadmits_0708_table[order(EDadmits_0708_table$patId),]
save(EDadmits_0708_table,file='EDadmits_0708_table.RData')

drugs <- combined_data_0708[combined_data_0708$typeId==2,]
drugs <- drugs[drugs$patId %in% intersect(keepPats0007,keepPats0815),]
dates <- as.Date(drugs$date)
drugs_0708 <- drugs[dates-as.Date('2007-01-01') >= 0 & dates-as.Date('2008-01-01') < 0,]
drugs_table <- as.data.frame(table(drugs_0708$patId))
colnames(drugs_table) <- c('patId','Count'); drugs_table$patId <- as.integer(as.character(drugs_table$patId))
noDrugspats <- as.data.frame(cbind(setdiff(intersect(keepPats0007,keepPats0815),unique(drugs_table$patId)),0))
colnames(noDrugspats) <- c('patId','Count'); noDrugspats$patId <- as.integer(noDrugspats$patId)
drugs_table <- rbind(drugs_table,noDrugspats)
drugs_table <- drugs_table[order(drugs_table$patId),]
save(drugs_table,file='drugs_table_0708.RData')


# LASSO models
library(glmnet)

metaphe_500 <- read.csv("combinedData_firstHalf_combinedData_firstHalf_JCVB0_nmar_K500_iter200_metaphe.csv",header=FALSE)
pheId_metadata_firstHalf <- read.csv('pheId_metadata_firstHalf.csv')[,-1]
metaphe_IDs <- unlist(read.csv('MTL_mixEHR_sureLDA_IDs_patId.csv',header=FALSE))

patients <- read.csv('/home/yueli/data/patients.csv')
patients <- patients[order(patients$id),c('id','month_of_birth','sex')]
deprivation <- read.csv('/home/yueli/data/patients_deprivation.csv')
deprivation <- deprivation[order(deprivation$id),c(1,3,4)]
load('patient_demographics.RData')
patients_augmented$sexInt <- as.integer(patients_augmented$sex=='F')
admissions <- read.csv('admissions_0708.csv')
drugs <- read.csv('drugs_0708.csv')
EDadmits <- read.csv('EDadmits_0708.csv')
mortality <- read.csv('mortality.csv')
HU <- read.csv('HU.csv')
IDs <- unlist(read.csv('MTL_mixEHR_sureLDA_predictions_withactes_noStopWords_firstHalf_patId.csv',header=FALSE))
theta <- read.table('MTL_mixEHR_sureLDA_predictions_withactes_noStopWords_firstHalf.txt',sep=',',header=FALSE)

mortIDs <- intersect(intersect(intersect(IDs,mortality$patId),HU$patId),patients_augmented$id)
admitIDs <- intersect(intersect(intersect(IDs,admissions$patId),HU$patId),patients_augmented$id)
drugIDs <- intersect(intersect(intersect(IDs,drugs$patId),HU$patId),patients_augmented$id)
EDIDs <- intersect(intersect(intersect(IDs,EDadmits$patId),HU$patId),patients_augmented$id)
patIDs <- intersect(intersect(IDs,HU$patId),patients_augmented$id)

metaphe_mort_matches <- which(metaphe_IDs %in% mortIDs)
metaphe_admit_matches <- which(metaphe_IDs %in% admitIDs)
metaphe_drug_matches <- which(metaphe_IDs %in% drugIDs)
metaphe_ED_matches <- which(metaphe_IDs %in% EDIDs)

theta_admit_matches <- which(IDs %in% admitIDs)
HU_admit_matches <- which(HU$patId %in% admitIDs)
pat_admit_matches <- which(patients_augmented$id %in% admitIDs)
theta_mort_matches <- which(IDs %in% mortIDs)
HU_mort_matches <- which(HU$patId %in% mortIDs)
pat_mort_matches <- which(patients_augmented$id %in% mortIDs)
theta_ED_matches <- which(IDs %in% EDIDs)
HU_ED_matches <- which(HU$patId %in% EDIDs)
pat_ED_matches <- which(patients_augmented$id %in% EDIDs)
theta_drug_matches <- which(IDs %in% drugIDs)
HU_drug_matches <- which(HU$patId %in% drugIDs)
pat_drug_matches <- which(patients_augmented$id %in% drugIDs)

theta_matches <- which(IDs %in% patIDs)
HU_matches <- which(HU$patId %in% patIDs)
pat_matches <- which(patients_augmented$id %in% patIDs)
mort_matches <- which(mortality$patId %in%patIDs)
admit_matches <- which(admissions$patId %in%patIDs)
drug_matches <- which(drugs$patId %in%patIDs)
ED_matches <- which(EDadmits$patId %in% patIDs)

mort_mort_matches <- which(mortality$patId %in% mortIDs)
admit_admit_matches <- which(admissions$patId %in% admitIDs)
ED_ED_matches <- which(EDadmits$patId %in% EDIDs)

X <- cbind(log(theta[theta_mort_matches,]+1), log(HU$Total[HU_mort_matches]+1), patients_augmented[pat_mort_matches,4:7])
X <- Matrix(as.matrix(X),sparse=TRUE)
mortality.logtheta.glmnet <- cv.glmnet(X,mortality$dead[mort_mort_matches],family='binomial',type.measure='auc',nfolds=5)
save(mortality.logtheta.glmnet,file='mortality.logtheta.glmnet2.RData')
mortality.logtheta.glmnet.AUC <- max(mortality.logtheta.glmnet$cvm)
save(mortality.logtheta.glmnet.AUC,file='mortality.logtheta.glmnet2.AUC.RData')
mortality.logtheta.glmnet.coefs <- coef(mortality.logtheta.glmnet)[-1]
save(mortality.logtheta.glmnet.coefs,file='mortality.logtheta.glmnet2.coefs.RData')

X <- cbind(log(theta[theta_admit_matches,]+1), log(HU$Total[HU_admit_matches]+1), patients_augmented[pat_admit_matches,4:7])
X <- Matrix(as.matrix(X),sparse=TRUE)
admissions.logtheta.glmnet <- cv.glmnet(X,admissions$admission[admit_admit_matches],family='binomial',type.measure='auc')
save(admissions.logtheta.glmnet,file='admissions.logtheta.glmnet2.RData')
admissions.logtheta.glmnet.AUC <- max(admissions.logtheta.glmnet$cvm)
save(admissions.logtheta.glmnet.AUC,file='admissions.logtheta.glmnet2.AUC.RData')
admissions.logtheta.glmnet.coefs <- coef(admissions.logtheta.glmnet)[-1]
save(admissions.logtheta.glmnet.coefs,file='admissions.logtheta.glmnet2.coefs.RData')

X <- cbind(log(theta[theta_ED_matches,]+1), log(HU$Total[HU_ED_matches]+1), patients_augmented[pat_ED_matches,4:7])
X <- Matrix(as.matrix(X),sparse=TRUE)
EDadmits.logtheta.glmnet <- cv.glmnet(X,EDadmits$Count[ED_ED_matches]>0,family='binomial',type.measure='auc')
save(EDadmits.logtheta.glmnet,file='EDadmits.logtheta.glmnet2.RData')
EDadmits.logtheta.glmnet.AUC <- max(EDadmits.logtheta.glmnet$cvm)
save(EDadmits.logtheta.glmnet.AUC,file='EDadmits.logtheta.glmnet2.AUC.RData')
EDadmits.logtheta.glmnet.coefs <- coef(EDadmits.logtheta.glmnet)[-1]
save(EDadmits.logtheta.glmnet.coefs,file='EDadmits.logtheta.glmnet2.coefs.RData')


## Same with PheCode counts

X <- cbind(log(pheCode_dat_firstHalf_uncompressed[as.character(mortIDs),]+1), log(HU$Total[HU_mort_matches]+1), patients_augmented[pat_mort_matches,4:7])
X <- Matrix(as.matrix(X),sparse=TRUE)
mortality.logtheta.glmnet <- cv.glmnet(X,mortality$dead[mort_mort_matches],family='binomial',type.measure='auc',nfolds=5)
save(mortality.logtheta.glmnet,file='mortality.logtheta.glmnet.PheCodes.RData')
mortality.logtheta.glmnet.AUC <- max(mortality.logtheta.glmnet$cvm)
save(mortality.logtheta.glmnet.AUC,file='mortality.logtheta.glmnet.PheCodes.AUC.RData')
mortality.logtheta.glmnet.coefs <- coef(mortality.logtheta.glmnet)[-1]
save(mortality.logtheta.glmnet.coefs,file='mortality.logtheta.glmnet.PheCodes.coefs.RData')

X <- cbind(log(pheCode_dat_firstHalf_uncompressed[as.character(admitIDs),]+1), log(HU$Total[HU_admit_matches]+1), patients_augmented[pat_admit_matches,4:7])
X <- Matrix(as.matrix(X),sparse=TRUE)
admissions.logtheta.glmnet <- cv.glmnet(X,admissions$admission[admit_admit_matches],family='binomial',type.measure='auc')
save(admissions.logtheta.glmnet,file='admissions.logtheta.glmnet.PheCodes.RData')
admissions.logtheta.glmnet.AUC <- max(admissions.logtheta.glmnet$cvm)
save(admissions.logtheta.glmnet.AUC,file='admissions.logtheta.glmnet.PheCodes.AUC.RData')
admissions.logtheta.glmnet.coefs <- coef(admissions.logtheta.glmnet)[-1]
save(admissions.logtheta.glmnet.coefs,file='admissions.logtheta.glmnet.PheCodes.coefs.RData')

X <- cbind(log(pheCode_dat_firstHalf_uncompressed[as.character(EDIDs),]+1), log(HU$Total[HU_ED_matches]+1), patients_augmented[pat_ED_matches,4:7])
X <- Matrix(as.matrix(X),sparse=TRUE)
EDadmits.logtheta.glmnet <- cv.glmnet(X,EDadmits$Count[ED_ED_matches]>0,family='binomial',type.measure='auc')
save(EDadmits.logtheta.glmnet,file='EDadmits.logtheta.glmnet.PheCodes.RData')
EDadmits.logtheta.glmnet.AUC <- max(EDadmits.logtheta.glmnet$cvm)
save(EDadmits.logtheta.glmnet.AUC,file='EDadmits.logtheta.glmnet.PheCodes.AUC.RData')
EDadmits.logtheta.glmnet.coefs <- coef(EDadmits.logtheta.glmnet)[-1]
save(EDadmits.logtheta.glmnet.coefs,file='EDadmits.logtheta.glmnet.PheCodes.coefs.RData')


## Same with MAP priors
load('map_prior_firstHalf.RData')
map_prior_firstHalf <- map_prior_firstHalf[,!is.na(map_prior_firstHalf[1,])]

X <- cbind(map_prior_firstHalf[as.character(mortIDs),], log(HU$Total[HU_mort_matches]+1), patients_augmented[pat_mort_matches,4:7])
X <- Matrix(as.matrix(X),sparse=TRUE)
mortality.logtheta.glmnet <- cv.glmnet(X,mortality$dead[mort_mort_matches],family='binomial',type.measure='auc',nfolds=5)
save(mortality.logtheta.glmnet,file='mortality.logtheta.glmnet.MAP.RData')
mortality.logtheta.glmnet.AUC <- max(mortality.logtheta.glmnet$cvm)
save(mortality.logtheta.glmnet.AUC,file='mortality.logtheta.glmnet.MAP.AUC.RData')
mortality.logtheta.glmnet.coefs <- coef(mortality.logtheta.glmnet)[-1]
save(mortality.logtheta.glmnet.coefs,file='mortality.logtheta.glmnet.MAP.coefs.RData')

X <- cbind(map_prior_firstHalf[as.character(admitIDs),], log(HU$Total[HU_admit_matches]+1), patients_augmented[pat_admit_matches,4:7])
X <- Matrix(as.matrix(X),sparse=TRUE)
admissions.logtheta.glmnet <- cv.glmnet(X,admissions$admission[admit_admit_matches],family='binomial',type.measure='auc')
save(admissions.logtheta.glmnet,file='admissions.logtheta.glmnet.MAP.RData')
admissions.logtheta.glmnet.AUC <- max(admissions.logtheta.glmnet$cvm)
save(admissions.logtheta.glmnet.AUC,file='admissions.logtheta.glmnet.MAP.AUC.RData')
admissions.logtheta.glmnet.coefs <- coef(admissions.logtheta.glmnet)[-1]
save(admissions.logtheta.glmnet.coefs,file='admissions.logtheta.glmnet.MAP.coefs.RData')

X <- cbind(map_prior_firstHalf[as.character(EDIDs),], log(HU$Total[HU_ED_matches]+1), patients_augmented[pat_ED_matches,4:7])
X <- Matrix(as.matrix(X),sparse=TRUE)
EDadmits.logtheta.glmnet <- cv.glmnet(X,EDadmits$Count[ED_ED_matches]>0,family='binomial',type.measure='auc')
save(EDadmits.logtheta.glmnet,file='EDadmits.logtheta.glmnet.MAP.RData')
EDadmits.logtheta.glmnet.AUC <- max(EDadmits.logtheta.glmnet$cvm)
save(EDadmits.logtheta.glmnet.AUC,file='EDadmits.logtheta.glmnet.MAP.AUC.RData')
EDadmits.logtheta.glmnet.coefs <- coef(EDadmits.logtheta.glmnet)[-1]
save(EDadmits.logtheta.glmnet.coefs,file='EDadmits.logtheta.glmnet.MAP.coefs.RData')


## Finally, same for mixEHR only metaphenotypes
metaphe_500 <- read.csv("combinedData_firstHalf_combinedData_firstHalf_JCVB0_nmar_K500_iter200_metaphe.csv",header=FALSE)
rownames(metaphe_500) <- unlist(read.csv('MTL_mixEHR_sureLDA_IDs_patId.csv',header=FALSE))

X <- cbind(log(metaphe_500[as.character(mortIDs),]+1), log(HU$Total[HU_mort_matches]+1), patients_augmented[pat_mort_matches,4:7])
X <- Matrix(as.matrix(X),sparse=TRUE)
mortality.logtheta.glmnet <- cv.glmnet(X,mortality$dead[mort_mort_matches],family='binomial',type.measure='auc',nfolds=5)
save(mortality.logtheta.glmnet,file='mortality.logtheta.glmnet.mixEHR500.RData')
mortality.logtheta.glmnet.AUC <- max(mortality.logtheta.glmnet$cvm)
save(mortality.logtheta.glmnet.AUC,file='mortality.logtheta.glmnet.mixEHR500.AUC.RData')
mortality.logtheta.glmnet.coefs <- coef(mortality.logtheta.glmnet)[-1]
save(mortality.logtheta.glmnet.coefs,file='mortality.logtheta.glmnet.mixEHR500.coefs.RData')

X <- cbind(log(metaphe_500[as.character(admitIDs),]+1), log(HU$Total[HU_admit_matches]+1), patients_augmented[pat_admit_matches,4:7])
X <- Matrix(as.matrix(X),sparse=TRUE)
admissions.logtheta.glmnet <- cv.glmnet(X,admissions$admission[admit_admit_matches],family='binomial',type.measure='auc')
save(admissions.logtheta.glmnet,file='admissions.logtheta.glmnet.mixEHR500.RData')
admissions.logtheta.glmnet.AUC <- max(admissions.logtheta.glmnet$cvm)
save(admissions.logtheta.glmnet.AUC,file='admissions.logtheta.glmnet.mixEHR500.AUC.RData')
admissions.logtheta.glmnet.coefs <- coef(admissions.logtheta.glmnet)[-1]
save(admissions.logtheta.glmnet.coefs,file='admissions.logtheta.glmnet.mixEHR500.coefs.RData')

X <- cbind(log(metaphe_500[as.character(EDIDs),]+1), log(HU$Total[HU_ED_matches]+1), patients_augmented[pat_ED_matches,4:7])
X <- Matrix(as.matrix(X),sparse=TRUE)
EDadmits.logtheta.glmnet <- cv.glmnet(X,EDadmits$Count[ED_ED_matches]>0,family='binomial',type.measure='auc')
save(EDadmits.logtheta.glmnet,file='EDadmits.logtheta.glmnet.mixEHR500.RData')
EDadmits.logtheta.glmnet.AUC <- max(EDadmits.logtheta.glmnet$cvm)
save(EDadmits.logtheta.glmnet.AUC,file='EDadmits.logtheta.glmnet.mixEHR500.AUC.RData')
EDadmits.logtheta.glmnet.coefs <- coef(EDadmits.logtheta.glmnet)[-1]
save(EDadmits.logtheta.glmnet.coefs,file='EDadmits.logtheta.glmnet.mixEHR500.coefs.RData')



# Log odds ratio of phenotype != 0 | mortality/admission/drugs/EDadmits
theta_mort <- theta[mort_matches,]>0
mortality_ORs <- apply(theta_mort,2,function(theta_i){
  glm(mortality$dead ~ theta_i, family='binomial')$coefficients[-1]
})
save(mortality_ORs,file='mortality_ORs.RData')

theta_mort <- theta[matches,]>0
admissions_ORs <- apply(theta_mort,2,function(theta_i){
  glm(admissions$admission ~ theta_i, family='binomial')$coefficients[-1]
})
save(mortality_ORs,file='admissions_ORs.RData')

EDadmits_ORs <- apply(theta_mort,2,function(theta_i){
  glm(I(EDadmits$Count>0) ~ theta_i, family='binomial')$coefficients[-1]
})
save(EDadmits_ORs,file='mortality_ORs.RData')

drugs_ORs <- apply(theta_mort,2,function(theta_i){
  glm(I(drugs$Count>0) ~ theta_i, family='binomial')$coefficients[-1]
})
save(drugs_ORs,file='drugs_ORs.RData')


## Interpret results

load('metadata.RData')
pheId_metadata <- read.csv('pheId_metadata.csv')[,-1]
pheId_metadata_firstHalf <- read.csv('pheId_metadata_firstHalf.csv')[,-1]
phecode_icd9_map <- read.csv('phecode_icd9_mappings_fixed102920.csv')
phecode_icd9_map <- phecode_icd9_map[!duplicated(phecode_icd9_map$PheCode),c('PheCode','Phenotype')]
phecode_icd9_map <- phecode_icd9_map[order(phecode_icd9_map$PheCode),]
phecode_icd9_map$Phenotype <- as.character(phecode_icd9_map$Phenotype)
pheId_metadata$Phenotype <- sapply(as.character(pheId_metadata$phecodes),function(cod){
  if (cod=='270.4'){NA}
  else if (substr(cod,1,1)!='V'){phecode_icd9_map$Phenotype[phecode_icd9_map$PheCode==as.numeric(cod)]}
  else{cod}
})
pheId_metadata_firstHalf$Phenotype <- sapply(as.character(pheId_metadata_firstHalf$phecodes),function(cod){
  if (cod=='270.4'){NA}
  else if (substr(cod,1,1)!='V'){phecode_icd9_map$Phenotype[phecode_icd9_map$PheCode==as.numeric(cod)]}
  else{cod}
})
write.csv(pheId_metadata,file='pheId_metadata.csv')
write.csv(pheId_metadata_firstHalf,file='pheId_metadata_firstHalf.csv')

pheId_metadata <- read.csv('pheId_metadata.csv')
pheId_metadata_firstHalf <- read.csv('pheId_metadata_firstHalf.csv')
load('mortality.logtheta.glmnet.AUC.RData')
load('mortality.logtheta.glmnet.coefs.RData')
load('admissions.logtheta.glmnet.AUC.RData')
load('admissions.logtheta.glmnet.coefs.RData')
load('drugs.logtheta.glmnet.coefs.RData')
load('EDadmits.logtheta.glmnet.coefs.RData')

# Mortality AUC: 0.931
# Admissions AUC: 0.773

theta_sds <- apply(theta,2,sd)[1:1326]

mortality.logtheta.glmnet.coefs <- mortality.logtheta.glmnet.coefs[1:1326]*theta_sds
mortality.logtheta.glmnet.coefs <- mortality.logtheta.glmnet.coefs/max(mortality.logtheta.glmnet.coefs)
idx <- order(mortality.logtheta.glmnet.coefs,decreasing=TRUE)[1:50]
mortality_topPredictors <- data.frame('coefs'=round(100*mortality.logtheta.glmnet.coefs[idx]),
                                      'phens'=gsub(',',' ',pheId_metadata_firstHalf$Phenotype[idx]))
write.csv(mortality_topPredictors,file='MTL_mortality_topPredictors_wordcloud_SDscaled.csv',row.names=FALSE)

admissions.logtheta.glmnet.coefs <- admissions.logtheta.glmnet.coefs[1:1326]*theta_sds
admissions.logtheta.glmnet.coefs <- admissions.logtheta.glmnet.coefs/max(admissions.logtheta.glmnet.coefs)
idx <- order(admissions.logtheta.glmnet.coefs,decreasing=TRUE)[1:50]
admissions_topPredictors <- data.frame('coefs'=round(100*admissions.logtheta.glmnet.coefs[idx]),
                                       'phens'=gsub(',',' ',pheId_metadata_firstHalf$Phenotype[idx]))
write.csv(admissions_topPredictors,file='MTL_admissions_topPredictors_wordcloud_SDscaled.csv',row.names=FALSE)

EDadmits.logtheta.glmnet.coefs <- EDadmits.logtheta.glmnet.coefs[1:1326]*theta_sds
EDadmits.logtheta.glmnet.coefs <- EDadmits.logtheta.glmnet.coefs/max(EDadmits.logtheta.glmnet.coefs)
idx <- order(EDadmits.logtheta.glmnet.coefs,decreasing=TRUE)[1:50]
EDadmits_topPredictors <- data.frame('coefs'=round(100*EDadmits.logtheta.glmnet.coefs[idx]),
                                     'phens'=gsub(',',' ',pheId_metadata_firstHalf$Phenotype[idx]))
write.csv(EDadmits_topPredictors,file='MTL_EDadmits_topPredictors_wordcloud_SDscaled.csv',row.names=FALSE)

drugs.logtheta.glmnet.coefs <- drugs.logtheta.glmnet.coefs[1:1326]*theta_sds
drugs.logtheta.glmnet.coefs <- drugs.logtheta.glmnet.coefs/max(drugs.logtheta.glmnet.coefs)
idx <- order(drugs.logtheta.glmnet.coefs,decreasing=TRUE)[1:50]
drugs_topPredictors <- data.frame('coefs'=round(100*drugs.logtheta.glmnet.coefs[idx]),
                                  'phens'=gsub(',',' ',pheId_metadata_firstHalf$Phenotype[idx]))
write.csv(drugs_topPredictors,file='MTL_drugs_topPredictors_wordcloud_SDscaled.csv',row.names=FALSE)


mortality.logtheta.glmnet.coefs <- mortality.logtheta.glmnet.coefs/min(mortality.logtheta.glmnet.coefs)
idx <- order(mortality.logtheta.glmnet.coefs,decreasing=TRUE)[1:50]
idx <- idx[mortality.logtheta.glmnet.coefs[idx]>0]
mortality_lowestPredictors <- data.frame('coefs'=round(100*mortality.logtheta.glmnet.coefs[idx]),
                                      'phens'=gsub(',',' ',pheId_metadata_firstHalf$Phenotype[idx]))
write.csv(mortality_lowestPredictors,file='MTL_mortality_lowestPredictors_wordcloud_SDscaled.csv',row.names=FALSE)

admissions.logtheta.glmnet.coefs <- admissions.logtheta.glmnet.coefs/min(admissions.logtheta.glmnet.coefs)
idx <- order(admissions.logtheta.glmnet.coefs,decreasing=TRUE)[1:50]
idx <- idx[admissions.logtheta.glmnet.coefs[idx]>0]
admissions_lowestPredictors <- data.frame('coefs'=round(100*admissions.logtheta.glmnet.coefs[idx]),
                                       'phens'=gsub(',',' ',pheId_metadata_firstHalf$Phenotype[idx]))
write.csv(admissions_lowestPredictors,file='MTL_admissions_lowestPredictors_wordcloud_SDscaled.csv',row.names=FALSE)

EDadmits.logtheta.glmnet.coefs <- EDadmits.logtheta.glmnet.coefs/min(EDadmits.logtheta.glmnet.coefs)
idx <- order(EDadmits.logtheta.glmnet.coefs,decreasing=TRUE)[1:50]
idx <- idx[EDadmits.logtheta.glmnet.coefs[idx]>0]
EDadmits_lowestPredictors <- data.frame('coefs'=round(100*EDadmits.logtheta.glmnet.coefs[idx]),
                                     'phens'=gsub(',',' ',pheId_metadata_firstHalf$Phenotype[idx]))
write.csv(EDadmits_lowestPredictors,file='MTL_EDadmits_lowestPredictors_wordcloud_SDscaled.csv',row.names=FALSE)

drugs.logtheta.glmnet.coefs <- drugs.logtheta.glmnet.coefs/min(drugs.logtheta.glmnet.coefs)
idx <- order(drugs.logtheta.glmnet.coefs,decreasing=TRUE)[1:50]
idx <- idx[drugs.logtheta.glmnet.coefs[idx]>0]
drugs_lowestPredictors <- data.frame('coefs'=round(100*drugs.logtheta.glmnet.coefs[idx]),
                                  'phens'=gsub(',',' ',pheId_metadata_firstHalf$Phenotype[idx]))
write.csv(drugs_lowestPredictors,file='MTL_drugs_lowestPredictors_wordcloud_SDscaled.csv',row.names=FALSE)

mortality_ORs <- mortality_ORs/max(mortality_ORs)
idx <- order(mortality_ORs,decreasing=TRUE)[1:50]
mortality_OR_topPredictors <- data.frame('coefs'=round(100*mortality_ORs[idx]),
                                         'phens'=gsub(',',' ',pheId_metadata_firstHalf$Phenotype[idx]))
write.csv(mortality_OR_topPredictors,file='mortality_OR_topPredictors.csv',row.names=FALSE)

admissions_ORs <- admissions_ORs/max(admissions_ORs)
idx <- order(admissions_ORs,decreasing=TRUE)[1:50]
admissions_OR_topPredictors <- data.frame('coefs'=round(100*admissions_ORs[idx]),
                                         'phens'=gsub(',',' ',pheId_metadata_firstHalf$Phenotype[idx]))
write.csv(admissions_OR_topPredictors,file='admissions_OR_topPredictors.csv',row.names=FALSE)

EDadmits_ORs <- EDadmits_ORs/max(EDadmits_ORs)
idx <- order(EDadmits_ORs,decreasing=TRUE)[1:50]
EDadmits_OR_topPredictors <- data.frame('coefs'=round(100*EDadmits_ORs[idx]),
                                         'phens'=gsub(',',' ',pheId_metadata_firstHalf$Phenotype[idx]))
write.csv(EDadmits_OR_topPredictors,file='EDadmits_OR_topPredictors.csv',row.names=FALSE)

drugs_ORs <- drugs_ORs/max(drugs_ORs,na.rm=T)
idx <- order(drugs_ORs,decreasing=TRUE)[1:50]
drugs_OR_topPredictors <- data.frame('coefs'=round(100*drugs_ORs[idx]),
                                         'phens'=gsub(',',' ',pheId_metadata_firstHalf$Phenotype[idx]))
write.csv(drugs_OR_topPredictors,file='drugs_OR_topPredictors.csv',row.names=FALSE)



## Compare to alternative predictors
library(glmnet)

metaphe_50 <- read.csv("combinedData_firstHalf_combinedData_firstHalf_JCVB0_nmar_K50_iter200_metaphe.csv",header=FALSE)
metaphe_100 <- read.csv("combinedData_firstHalf_combinedData_firstHalf_JCVB0_nmar_K100_iter200_metaphe.csv",header=FALSE)
metaphe_200 <- read.csv("combinedData_firstHalf_combinedData_firstHalf_JCVB0_nmar_K200_iter200_metaphe.csv",header=FALSE)
metaphe_500 <- read.csv("combinedData_firstHalf_combinedData_firstHalf_JCVB0_nmar_K500_iter200_metaphe.csv",header=FALSE)

pheId_metadata <- read.csv('pheId_metadata.csv')[,-1]
pheId_metadata_firstHalf <- read.csv('pheId_metadata_firstHalf.csv')[,-1]

admissions <- read.csv('admissions_0708.csv')
drugs <- read.csv('drugs_0708.csv')
EDadmits <- read.csv('EDadmits_0708.csv')
mortality <- read.csv('mortality.csv')
HU <- read.csv('HU.csv')
IDs <- unlist(read.csv('MTL_mixEHR_sureLDA_IDs_patId.csv',header=FALSE))

matches <- which(IDs %in% admissions$patId)
HU_matches <- which(HU$patId %in% admissions$patId)
mort_matches <- which(IDs %in% mortality$patId & IDs %in% HU$patId)
HU_mort_matches <- which(HU$patId %in% mortality$patId)

theta <- log(metaphe_50[mort_matches,]+1); theta$HU <- log(HU$Total[HU_mort_matches]+1)
theta <- Matrix(as.matrix(theta),sparse=TRUE)
mortality.logtheta.glmnet <- cv.glmnet(theta,mortality$dead,family='binomial',type.measure='auc')
save(mortality.logtheta.glmnet,file='mortality.logtheta.glmnet.50.RData')
mortality.logtheta.glmnet.AUC <- max(mortality.logtheta.glmnet$cvm)
save(mortality.logtheta.glmnet.AUC,file='mortality.logtheta.glmnet.50.AUC.RData')
mortality.logtheta.glmnet.coefs <- coef(mortality.logtheta.glmnet)[-1]
save(mortality.logtheta.glmnet.coefs,file='mortality.logtheta.glmnet.50.coefs.RData')

theta <- log(metaphe_100[mort_matches,]+1); theta$HU <- log(HU$Total[HU_mort_matches]+1)
theta <- Matrix(as.matrix(theta),sparse=TRUE)
mortality.logtheta.glmnet <- cv.glmnet(theta,mortality$dead,family='binomial',type.measure='auc')
save(mortality.logtheta.glmnet,file='mortality.logtheta.glmnet.100.RData')
mortality.logtheta.glmnet.AUC <- max(mortality.logtheta.glmnet$cvm)
save(mortality.logtheta.glmnet.AUC,file='mortality.logtheta.glmnet.100.AUC.RData')
mortality.logtheta.glmnet.coefs <- coef(mortality.logtheta.glmnet)[-1]
save(mortality.logtheta.glmnet.coefs,file='mortality.logtheta.glmnet.100.coefs.RData')

theta <- log(metaphe_200[mort_matches,]+1); theta$HU <- log(HU$Total[HU_mort_matches]+1)
theta <- Matrix(as.matrix(theta),sparse=TRUE)
mortality.logtheta.glmnet <- cv.glmnet(theta,mortality$dead,family='binomial',type.measure='auc')
save(mortality.logtheta.glmnet,file='mortality.logtheta.glmnet.200.RData')
mortality.logtheta.glmnet.AUC <- max(mortality.logtheta.glmnet$cvm)
save(mortality.logtheta.glmnet.AUC,file='mortality.logtheta.glmnet.200.AUC.RData')
mortality.logtheta.glmnet.coefs <- coef(mortality.logtheta.glmnet)[-1]
save(mortality.logtheta.glmnet.coefs,file='mortality.logtheta.glmnet.200.coefs.RData')

theta <- log(metaphe_500[mort_matches,]+1); theta$HU <- log(HU$Total[HU_mort_matches]+1)
theta <- Matrix(as.matrix(theta),sparse=TRUE)
mortality.logtheta.glmnet <- cv.glmnet(theta,mortality$dead,family='binomial',type.measure='auc')
save(mortality.logtheta.glmnet,file='mortality.logtheta.glmnet.500.RData')
mortality.logtheta.glmnet.AUC <- max(mortality.logtheta.glmnet$cvm)
save(mortality.logtheta.glmnet.AUC,file='mortality.logtheta.glmnet.500.AUC.RData')
mortality.logtheta.glmnet.coefs <- coef(mortality.logtheta.glmnet)[-1]
save(mortality.logtheta.glmnet.coefs,file='mortality.logtheta.glmnet.500.coefs.RData')

theta <- log(metaphe_50[matches,]+1); theta$HU <- log(HU$Total[matches]+1)
theta <- Matrix(as.matrix(theta),sparse=TRUE)
readmission.logtheta.glmnet <- cv.glmnet(theta,admissions$admission,family='binomial',type.measure='auc')
save(readmission.logtheta.glmnet,file='readmission.logtheta.glmnet.50.RData')
readmission.logtheta.glmnet.AUC <- max(readmission.logtheta.glmnet$cvm)
save(readmission.logtheta.glmnet.AUC,file='readmission.logtheta.glmnet.50.AUC.RData')
readmission.logtheta.glmnet.coefs <- coef(readmission.logtheta.glmnet)[-1]
save(readmission.logtheta.glmnet.coefs,file='readmission.logtheta.glmnet.50.coefs.RData')

theta <- log(metaphe_100[matches,]+1); theta$HU <- log(HU$Total[matches]+1)
theta <- Matrix(as.matrix(theta),sparse=TRUE)
readmission.logtheta.glmnet <- cv.glmnet(theta,admissions$admission,family='binomial',type.measure='auc')
save(readmission.logtheta.glmnet,file='readmission.logtheta.glmnet.100.RData')
readmission.logtheta.glmnet.AUC <- max(readmission.logtheta.glmnet$cvm)
save(readmission.logtheta.glmnet.AUC,file='readmission.logtheta.glmnet.100.AUC.RData')
readmission.logtheta.glmnet.coefs <- coef(readmission.logtheta.glmnet)[-1]
save(readmission.logtheta.glmnet.coefs,file='readmission.logtheta.glmnet.100.coefs.RData')

theta <- log(metaphe_200[matches,]+1); theta$HU <- log(HU$Total[matches]+1)
theta <- Matrix(as.matrix(theta),sparse=TRUE)
readmission.logtheta.glmnet <- cv.glmnet(theta,admissions$admission,family='binomial',type.measure='auc')
save(readmission.logtheta.glmnet,file='readmission.logtheta.glmnet.200.RData')
readmission.logtheta.glmnet.AUC <- max(readmission.logtheta.glmnet$cvm)
save(readmission.logtheta.glmnet.AUC,file='readmission.logtheta.glmnet.200.AUC.RData')
readmission.logtheta.glmnet.coefs <- coef(readmission.logtheta.glmnet)[-1]
save(readmission.logtheta.glmnet.coefs,file='readmission.logtheta.glmnet.200.coefs.RData')

theta <- log(metaphe_500[matches,]+1); theta$HU <- log(HU$Total[matches]+1)
theta <- Matrix(as.matrix(theta),sparse=TRUE)
readmission.logtheta.glmnet <- cv.glmnet(theta,admissions$admission,family='binomial',type.measure='auc')
save(readmission.logtheta.glmnet,file='readmission.logtheta.glmnet.500.RData')
readmission.logtheta.glmnet.AUC <- max(readmission.logtheta.glmnet$cvm)
save(readmission.logtheta.glmnet.AUC,file='readmission.logtheta.glmnet.500.AUC.RData')
readmission.logtheta.glmnet.coefs <- coef(readmission.logtheta.glmnet)[-1]
save(readmission.logtheta.glmnet.coefs,file='readmission.logtheta.glmnet.500.coefs.RData')


# PCA-LASSO
library(prcomp)

load('combined_data_firstHalf.RData')
load('combinedData_firstHalf_mixEHRready.RData')
combined_data_firstHalf$typepheId <- 1000000*combined_data_firstHalf$typeId + combined_data_firstHalf$pheId
uniqueCases_firstHalf$typepheId <- 1000000*uniqueCases_firstHalf$typeId + uniqueCases_firstHalf$pheId
uniqueCases_firstHalf$logfreq <- log(uniqueCases_firstHalf$freq+1)
uniqueCases_firstHalf <- uniqueCases_firstHalf[order(uniqueCases_firstHalf$patId),]
# dat_firstHalf <- table(combined_data_firstHalf[,c('patId','typepheId')])
# dat_firstHalf_pcs <- prcomp(log(dat_firstHalf+1),scale.=TRUE)$x

dat_firstHalf <- with(uniqueCases_firstHalf, sparseMatrix(i=as.numeric(as.factor(patId)),
                                                          j=as.numeric(as.factor(typepheId)),
                                                          x=as.numeric(logfreq)))
pca <- prcomp(dat_firstHalf, scale.=TRUE)
pca_projection <- pca$x




## NEW ANALYSES 11/25/2020
library(glmnet)

patients <- read.csv('/home/yueli/data/patients.csv')
deprivation <- read.csv('/home/yueli/data/patients_deprivation.csv')
patients <- patients[order(patients$id),c('id','month_of_birth','sex')]
deprivation <- deprivation[order(deprivation$id),c(1,3,4)]

# patients_augmented <- merge(patients,deprivation,by='id',all=FALSE)
# patients_augmented$age <- as.numeric(as.Date('2007-01-01') - as.Date(paste0(patients_augmented$month_of_birth,'-01')))/365.25
# save(patients_augmented, file='patient_demographics.RData')
load('patient_demographics.RData')

admissions <- read.csv('admissions_0708.csv')
drugs <- read.csv('drugs_0708.csv')
EDadmits <- read.csv('EDadmits_0708.csv')
mortality <- read.csv('mortality.csv')
HU <- read.csv('HU.csv')
IDs <- unlist(read.csv('MTL_mixEHR_sureLDA_predictions_withactes_noStopWords_firstHalf_patId.csv',header=FALSE))
theta <- read.table('MTL_mixEHR_sureLDA_predictions_withactes_noStopWords_firstHalf.txt',sep=',',header=FALSE)
mortIDs <- intersect(intersect(intersect(IDs,mortality$patId),HU$patId),patients_augmented$id)
admitIDs <- intersect(intersect(intersect(IDs,admissions$patId),HU$patId),patients_augmented$id)
drugIDs <- intersect(intersect(intersect(IDs,drugs$patId),HU$patId),patients_augmented$id)
EDIDs <- intersect(intersect(intersect(IDs,EDadmits$patId),HU$patId),patients_augmented$id)
patIDs <- intersect(intersect(IDs,HU$patId),patients_augmented$id)

theta_admit_matches <- which(IDs %in% admitIDs)
HU_admit_matches <- which(HU$patId %in% admitIDs)
pat_admit_matches <- which(patients_augmented$id %in% admitIDs)
theta_mort_matches <- which(IDs %in% mortIDs)
HU_mort_matches <- which(HU$patId %in% mortIDs)
pat_mort_matches <- which(patients_augmented$id %in% mortIDs)
theta_ED_matches <- which(IDs %in% EDIDs)
HU_ED_matches <- which(HU$patId %in% EDIDs)
pat_ED_matches <- which(patients_augmented$id %in% EDIDs)
theta_drug_matches <- which(IDs %in% drugIDs)
HU_drug_matches <- which(HU$patId %in% drugIDs)
pat_drug_matches <- which(patients_augmented$id %in% drugIDs)

theta_matches <- which(IDs %in% patIDs)
HU_matches <- which(HU$patId %in% patIDs)
pat_matches <- which(patients_augmented$id %in% patIDs)
mort_matches <- which(mortality$patId %in%patIDs)
admit_matches <- which(admissions$patId %in%patIDs)
drug_matches <- which(drugs$patId %in%patIDs)
ED_matches <- which(EDadmits$patId %in% patIDs)


# P(socioeconomic status | theta)
X <- as.matrix(cbind(log(theta[theta_matches,]+1), log(HU$Total[HU_matches]+1),
           as.integer(patients_augmented$sex[pat_matches]=='F')))
theta.material_0020 <- cv.glmnet(X,patients_augmented$material_deprivation[pat_matches],
                                 weights=I(patients_augmented$age[pat_matches]>=0 & patients_augmented$age[pat_matches]<20))
theta.material_2040 <- cv.glmnet(X,patients_augmented$material_deprivation[pat_matches],
                                 weights=I(patients_augmented$age[pat_matches]>=20 & patients_augmented$age[pat_matches]<40))
theta.material_4060 <- cv.glmnet(X,patients_augmented$material_deprivation[pat_matches],
                                 weights=I(patients_augmented$age[pat_matches]>=40 & patients_augmented$age[pat_matches]<60))
theta.material_6080 <- cv.glmnet(X,patients_augmented$material_deprivation[pat_matches],
                                 weights=I(patients_augmented$age[pat_matches]>=60 & patients_augmented$age[pat_matches]<80))
theta.material_8000 <- cv.glmnet(X,patients_augmented$material_deprivation[pat_matches],
                                 weights=I(patients_augmented$age[pat_matches]>=80 & patients_augmented$age[pat_matches]<100))
theta.social_0020 <- cv.glmnet(X,patients_augmented$social_deprivation[pat_matches],
                                 weights=I(patients_augmented$age[pat_matches]>=0 & patients_augmented$age[pat_matches]<20))
theta.social_2040 <- cv.glmnet(X,patients_augmented$social_deprivation[pat_matches],
                                 weights=I(patients_augmented$age[pat_matches]>=20 & patients_augmented$age[pat_matches]<40))
theta.social_4060 <- cv.glmnet(X,patients_augmented$social_deprivation[pat_matches],
                                 weights=I(patients_augmented$age[pat_matches]>=40 & patients_augmented$age[pat_matches]<60))
theta.social_6080 <- cv.glmnet(X,patients_augmented$social_deprivation[pat_matches],
                                 weights=I(patients_augmented$age[pat_matches]>=60 & patients_augmented$age[pat_matches]<80))
theta.social_8000 <- cv.glmnet(X,patients_augmented$social_deprivation[pat_matches],
                                 weights=I(patients_augmented$age[pat_matches]>=80))

thetas = list(coef(theta.material_0020),coef(theta.material_2040),coef(theta.material_4060),coef(theta.material_6080),coef(theta.material_8000),
              coef(theta.social_0020),coef(theta.social_2040),coef(theta.social_4060),coef(theta.social_6080),coef(theta.social_8000))
save(thetas,file='theta.SES.predictor.coefficients.RData')


load('theta.SES.predictor.coefficients.RData')
load('metadata.RData')
pheId_metadata_firstHalf <- read.csv('pheId_metadata_firstHalf.csv')[,-1]
theta.ses.predictors <- foreach(theta=thetas) %do% {
  theta <- theta[2:1327]
  theta_max <- theta/max(theta);   theta_min <- theta/min(theta)
  idx_max <- order(theta_max,decreasing=TRUE)[1:50]
  idx_min <- order(theta_min,decreasing=TRUE)[1:50]
  predictors <- data.frame('coefs_max'=round(100*theta_max[idx_max]), 'phens_max'=gsub(',',' ',pheId_metadata_firstHalf$Phenotype[idx_max]),
                           'coefs_min'=round(100*theta_min[idx_min]), 'phens_min'=gsub(',',' ',pheId_metadata_firstHalf$Phenotype[idx_min]))
  predictors
}
for (i in 1:length(theta.ses.predictors)){
  write.csv(theta.ses.predictors[[i]][,1:2],file=paste0('theta_ses_predictor_',i,'_max.csv'),row.names=FALSE)
  write.csv(theta.ses.predictors[[i]][,3:4],file=paste0('theta_ses_predictor_',i,'_min.csv'),row.names=FALSE)
}


# P(mortality | theta, socioeconomic status)
X <- as.matrix(cbind(log(theta[theta_mort_matches,]+1), log(HU$Total[HU_mort_matches]+1),
                     as.integer(patients_augmented$sex[pat_mort_matches]=='F')))
theta.material_low.mortality_0020 <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',
                                               weights=I(patients_augmented$age[pat_mort_matches]>=0 & patients_augmented$age[pat_mort_matches]<20 & patients_augmented$material_deprivation[pat_mort_matches]<0))
theta.material_low.mortality_2040 <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',
                                               weights=I(patients_augmented$age[pat_mort_matches]>=20 & patients_augmented$age[pat_mort_matches]<40 & patients_augmented$material_deprivation[pat_mort_matches]<0))
theta.material_low.mortality_4060 <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',
                                               weights=I(patients_augmented$age[pat_mort_matches]>=40 & patients_augmented$age[pat_mort_matches]<60 & patients_augmented$material_deprivation[pat_mort_matches]<0))
theta.material_low.mortality_6080 <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',
                                               weights=I(patients_augmented$age[pat_mort_matches]>=60 & patients_augmented$age[pat_mort_matches]<80 & patients_augmented$material_deprivation[pat_mort_matches]<0))
theta.material_low.mortality_8000 <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',
                                               weights=I(patients_augmented$age[pat_mort_matches]>=80 & patients_augmented$material_deprivation[pat_mort_matches]<0))
theta.material_high.mortality_0020 <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',
                                               weights=I(patients_augmented$age[pat_mort_matches]>=0 & patients_augmented$age[pat_mort_matches]<20 & patients_augmented$material_deprivation[pat_mort_matches]>0))
theta.material_high.mortality_2040 <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',
                                               weights=I(patients_augmented$age[pat_mort_matches]>=20 & patients_augmented$age[pat_mort_matches]<40 & patients_augmented$material_deprivation[pat_mort_matches]>0))
theta.material_high.mortality_4060 <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',
                                               weights=I(patients_augmented$age[pat_mort_matches]>=40 & patients_augmented$age[pat_mort_matches]<60 & patients_augmented$material_deprivation[pat_mort_matches]>0))
theta.material_high.mortality_6080 <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',
                                               weights=I(patients_augmented$age[pat_mort_matches]>=60 & patients_augmented$age[pat_mort_matches]<80 & patients_augmented$material_deprivation[pat_mort_matches]>0))
theta.material_high.mortality_8000 <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',
                                               weights=I(patients_augmented$age[pat_mort_matches]>=80 & patients_augmented$material_deprivation[pat_mort_matches]>0))
theta.social_low.mortality_0020 <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',
                                               weights=I(patients_augmented$age[pat_mort_matches]>=0 & patients_augmented$age[pat_mort_matches]<20 & patients_augmented$social_deprivation[pat_mort_matches]<0))
theta.social_low.mortality_2040 <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',
                                               weights=I(patients_augmented$age[pat_mort_matches]>=20 & patients_augmented$age[pat_mort_matches]<40 & patients_augmented$social_deprivation[pat_mort_matches]<0))
theta.social_low.mortality_4060 <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',
                                               weights=I(patients_augmented$age[pat_mort_matches]>=40 & patients_augmented$age[pat_mort_matches]<60 & patients_augmented$social_deprivation[pat_mort_matches]<0))
theta.social_low.mortality_6080 <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',
                                               weights=I(patients_augmented$age[pat_mort_matches]>=60 & patients_augmented$age[pat_mort_matches]<80 & patients_augmented$social_deprivation[pat_mort_matches]<0))
theta.social_low.mortality_8000 <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',
                                               weights=I(patients_augmented$age[pat_mort_matches]>=80 & patients_augmented$social_deprivation[pat_mort_matches]<0))
theta.social_high.mortality_0020 <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',
                                                weights=I(patients_augmented$age[pat_mort_matches]>=0 & patients_augmented$age[pat_mort_matches]<20 & patients_augmented$social_deprivation[pat_mort_matches]>0))
theta.social_high.mortality_2040 <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',
                                                weights=I(patients_augmented$age[pat_mort_matches]>=20 & patients_augmented$age[pat_mort_matches]<40 & patients_augmented$social_deprivation[pat_mort_matches]>0))
theta.social_high.mortality_4060 <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',
                                                weights=I(patients_augmented$age[pat_mort_matches]>=40 & patients_augmented$age[pat_mort_matches]<60 & patients_augmented$social_deprivation[pat_mort_matches]>0))
theta.social_high.mortality_6080 <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',
                                                weights=I(patients_augmented$age[pat_mort_matches]>=60 & patients_augmented$age[pat_mort_matches]<80 & patients_augmented$social_deprivation[pat_mort_matches]>0))
theta.social_high.mortality_8000 <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',
                                                weights=I(patients_augmented$age[pat_mort_matches]>=80 & patients_augmented$social_deprivation[pat_mort_matches]>0))

thetas = list(coef(theta.material_low.mortality_0020),coef(theta.material_low.mortality_2040),coef(theta.material_low.mortality_4060),coef(theta.material_low.mortality_6080),coef(theta.material_low.mortality_8000),
              coef(theta.material_high.mortality_0020),coef(theta.material_high.mortality_2040),coef(theta.material_high.mortality_4060),coef(theta.material_high.mortality_6080),coef(theta.material_high.mortality_8000),
              coef(theta.social_low.mortality_0020),coef(theta.social_low.mortality_2040),coef(theta.social_low.mortality_4060),coef(theta.social_low.mortality_6080),coef(theta.social_low.mortality_8000),
              coef(theta.social_high.mortality_0020),coef(theta.social_high.mortality_2040),coef(theta.social_high.mortality_4060),coef(theta.social_high.mortality_6080),coef(theta.social_high.mortality_8000))
save(thetas,file='theta.mortality.predictor.coefficients.RData')


load('theta.mortality.predictor.coefficients.RData')
load('metadata.RData')
pheId_metadata_firstHalf <- read.csv('pheId_metadata_firstHalf.csv')[,-1]
theta.sds <- apply(theta,2,sd)[1:1326]
theta.mortality.predictors <- foreach(theta=thetas) %do% {
  theta <- theta[2:1327]*theta.sds
  theta_max <- theta/max(theta);   theta_min <- theta/min(theta)
  idx_max <- order(theta_max,decreasing=TRUE)[1:50]
  idx_min <- order(theta_min,decreasing=TRUE)[1:50]
  predictors <- data.frame('coefs_max'=round(100*theta_max[idx_max]), 'phens_max'=gsub(',',' ',pheId_metadata_firstHalf$Phenotype[idx_max]),
                           'coefs_min'=round(100*theta_min[idx_min]), 'phens_min'=gsub(',',' ',pheId_metadata_firstHalf$Phenotype[idx_min]))
  predictors
}
for (i in 1:length(theta.mortality.predictors)){
  write.csv(theta.mortality.predictors[[i]][theta.mortality.predictors[[i]][,1]>0,1:2],
            file=paste0('theta_mortality_predictor_',i,'_max_scaled.csv'),row.names=FALSE)
  write.csv(theta.mortality.predictors[[i]][theta.mortality.predictors[[i]][,3]>0,3:4],
            file=paste0('theta_mortality_predictor_',i,'_min_scaled.csv'),row.names=FALSE)
}


# LASSO models adjusted for SES

X <- as.matrix(cbind(log(theta[theta_mort_matches,]+1), log(HU$Total[HU_mort_matches]+1),
                     as.numeric(patients_augmented$sex[pat_mort_matches]=='F')))
ages <- cbind(patients_augmented$age[pat_mort_matches]>=0 & patients_augmented$age[pat_mort_matches]<20,
              patients_augmented$age[pat_mort_matches]>=20 & patients_augmented$age[pat_mort_matches]<40,
              patients_augmented$age[pat_mort_matches]>=40 & patients_augmented$age[pat_mort_matches]<60,
              patients_augmented$age[pat_mort_matches]>=60 & patients_augmented$age[pat_mort_matches]<80,
              patients_augmented$age[pat_mort_matches]>=80)
ses <- as.matrix(patients_augmented[pat_mort_matches,4:5])
Xstar <- cbind(ses,X,X*ses[,1], X*ses[,2])

theta.ses.mortality_4060 <- cv.glmnet(Xstar[ages[,3],],mortality$dead[mort_matches[ages[,3]]],family='binomial',type.measure='auc')
theta.ses.mortality_4060.coefs <- coefficients(theta.ses.mortality_4060)
theta.ses.mortality_4060.auc <- max(theta.ses.mortality_4060$cvm)

theta.ses.mortality_6080 <- cv.glmnet(Xstar[ages[,4],],mortality$dead[mort_matches[ages[,4]]],family='binomial',type.measure='auc')
theta.ses.mortality_6080.coefs <- coefficients(theta.ses.mortality_6080)
theta.ses.mortality_6080.auc <- max(theta.ses.mortality_6080$cvm)

theta.ses.mortality_8000 <- cv.glmnet(Xstar[ages[,5],],mortality$dead[mort_matches[ages[,5]]],family='binomial',type.measure='auc')
theta.ses.mortality_8000.coefs <- coefficients(theta.ses.mortality_8000)
theta.ses.mortality_8000.auc <- max(theta.ses.mortality_8000$cvm)

results <- list(theta.ses.mortality_4060.coefs,theta.ses.mortality_4060.auc,theta.ses.mortality_6080.coefs,
                theta.ses.mortality_6080.auc,theta.ses.mortality_8000.coefs,theta.ses.mortality_8000.auc)
save(results,file='SES_mortality_results_120220.RData')



X <- as.matrix(cbind(log(theta[theta_admit_matches,]+1), log(HU$Total[HU_admit_matches]+1),
                     as.numeric(patients_augmented$sex[pat_admit_matches]=='F')))
ages <- cbind(patients_augmented$age[pat_admit_matches]>=0 & patients_augmented$age[pat_admit_matches]<20,
              patients_augmented$age[pat_admit_matches]>=20 & patients_augmented$age[pat_admit_matches]<40,
              patients_augmented$age[pat_admit_matches]>=40 & patients_augmented$age[pat_admit_matches]<60,
              patients_augmented$age[pat_admit_matches]>=60 & patients_augmented$age[pat_admit_matches]<80,
              patients_augmented$age[pat_admit_matches]>=80)
ses <- as.matrix(patients_augmented[pat_admit_matches,4:5])
Xstar <- cbind(ses,X,X*ses[,1], X*ses[,2])

theta.ses.admission_4060 <- cv.glmnet(Xstar[ages[,3],],admissions$admission[admit_matches[ages[,3]]],family='binomial',type.measure='auc')
theta.ses.admission_4060.coefs <- coefficients(theta.ses.admission_4060)
theta.ses.admission_4060.auc <- max(theta.ses.admission_4060$cvm)

theta.ses.admission_6080 <- cv.glmnet(Xstar[ages[,4],],admissions$admission[admit_matches[ages[,4]]],family='binomial',type.measure='auc')
theta.ses.admission_6080.coefs <- coefficients(theta.ses.admission_6080)
theta.ses.admission_6080.auc <- max(theta.ses.admission_6080$cvm)

theta.ses.admission_8000 <- cv.glmnet(Xstar[ages[,5],],admissions$admission[admit_matches[ages[,5]]],family='binomial',type.measure='auc')
theta.ses.admission_8000.coefs <- coefficients(theta.ses.admission_8000)
theta.ses.admission_8000.auc <- max(theta.ses.admission_8000$cvm)

results <- list(theta.ses.admission_4060.coefs,theta.ses.admission_4060.auc,theta.ses.admission_6080.coefs,
                theta.ses.admission_6080.auc,theta.ses.admission_8000.coefs,theta.ses.admission_8000.auc)
save(results,file='SES_admission_results_120220.RData')



X <- as.matrix(cbind(log(theta[theta_ED_matches,]+1), log(HU$Total[HU_ED_matches]+1),
                     as.numeric(patients_augmented$sex[pat_ED_matches]=='F')))
ages <- cbind(patients_augmented$age[pat_ED_matches]>=0 & patients_augmented$age[pat_ED_matches]<20,
              patients_augmented$age[pat_ED_matches]>=20 & patients_augmented$age[pat_ED_matches]<40,
              patients_augmented$age[pat_ED_matches]>=40 & patients_augmented$age[pat_ED_matches]<60,
              patients_augmented$age[pat_ED_matches]>=60 & patients_augmented$age[pat_ED_matches]<80,
              patients_augmented$age[pat_ED_matches]>=80)
ses <- as.matrix(patients_augmented[pat_ED_matches,4:5])
Xstar <- cbind(ses,X,X*ses[,1], X*ses[,2])

theta.ses.EDadmit_4060 <- cv.glmnet(Xstar[ages[,3],],EDadmits$Count[ED_matches[ages[,3]]],family='poisson')
theta.ses.EDadmit_4060.coefs <- coefficients(theta.ses.EDadmit_4060)

theta.ses.EDadmit_6080 <- cv.glmnet(Xstar[ages[,4],],EDadmits$Count[ED_matches[ages[,4]]],family='poisson')
theta.ses.EDadmit_6080.coefs <- coefficients(theta.ses.EDadmit_6080)

theta.ses.EDadmit_8000 <- cv.glmnet(Xstar[ages[,5],],EDadmits$Count[ED_matches[ages[,5]]],family='poisson')
theta.ses.EDadmit_8000.coefs <- coefficients(theta.ses.EDadmit_8000)

results <- list(theta.ses.EDadmit_4060.coefs,theta.ses.EDadmit_6080.coefs,theta.ses.EDadmit_8000.coefs)
save(results,file='SES_EDadmit_results_120220.RData')


X <- as.matrix(cbind(log(theta[theta_drug_matches,]+1), log(HU$Total[HU_drug_matches]+1),
                     as.numeric(patients_augmented$sex[pat_drug_matches]=='F')))
ages <- cbind(patients_augmented$age[pat_drug_matches]>=0 & patients_augmented$age[pat_drug_matches]<20,
              patients_augmented$age[pat_drug_matches]>=20 & patients_augmented$age[pat_drug_matches]<40,
              patients_augmented$age[pat_drug_matches]>=40 & patients_augmented$age[pat_drug_matches]<60,
              patients_augmented$age[pat_drug_matches]>=60 & patients_augmented$age[pat_drug_matches]<80,
              patients_augmented$age[pat_drug_matches]>=80)
ses <- as.matrix(patients_augmented[pat_drug_matches,4:5])
Xstar <- cbind(ses,X,X*ses[,1], X*ses[,2])

theta.ses.drugs_4060 <- cv.glmnet(Xstar[ages[,3],],drugs$Count[ED_matches[ages[,3]]],family='poisson')
theta.ses.drugs_4060.coefs <- coefficients(theta.ses.drugs_4060)

theta.ses.drugs_6080 <- cv.glmnet(Xstar[ages[,4],],drugs$Count[ED_matches[ages[,4]]],family='poisson')
theta.ses.drugs_6080.coefs <- coefficients(theta.ses.drugs_6080)

theta.ses.drugs_8000 <- cv.glmnet(Xstar[ages[,5],],drugs$Count[ED_matches[ages[,5]]],family='poisson')
theta.ses.drugs_8000.coefs <- coefficients(theta.ses.drugs_8000)

results <- list(theta.ses.drugs_4060.coefs,theta.ses.drugs_6080.coefs,theta.ses.drugs_8000.coefs)
save(results,file='SES_drugs_results_120220.RData')



load('SES_mortality_results_120220.RData')
load('metadata.RData')
pheId_metadata_firstHalf <- read.csv('pheId_metadata_firstHalf.csv')[,-1]
X <- as.matrix(cbind(log(theta[theta_mort_matches,]+1), log(HU$Total[HU_mort_matches]+1),
                     as.numeric(patients_augmented$sex[pat_mort_matches]=='F')))
ages <- cbind(patients_augmented$age[pat_mort_matches]>=0 & patients_augmented$age[pat_mort_matches]<20,
              patients_augmented$age[pat_mort_matches]>=20 & patients_augmented$age[pat_mort_matches]<40,
              patients_augmented$age[pat_mort_matches]>=40 & patients_augmented$age[pat_mort_matches]<60,
              patients_augmented$age[pat_mort_matches]>=60 & patients_augmented$age[pat_mort_matches]<80,
              patients_augmented$age[pat_mort_matches]>=80)
ses <- as.matrix(patients_augmented[pat_mort_matches,4:5])
Xstar <- cbind(X*ses[,1], X*ses[,2])
Xstar.sds <- apply(Xstar,2,sd)
thetas <- list(results[[1]][-c(1:1332)],results[[3]][-c(1:1332)],results[[5]][-c(1:1332)])
theta.mortality.predictors <- foreach(theta=thetas) %do% {
  theta_i <- theta_i*Xstar.sds
  theta_material_max <- theta_i[1:1326]/max(theta_i[1:1326]);   theta_material_min <- theta_i[1:1326]/min(theta_i[1:1326])
  theta_social_max <- theta_i[1330:2655]/max(theta_i[1330:2655]);   theta_social_min <- theta_i[1330:2655]/min(theta_i[1330:2655])
  idx_material_max <- order(theta_material_max,decreasing=TRUE)[1:50]
  idx_material_min <- order(theta_material_min,decreasing=TRUE)[1:50]
  idx_social_max <- order(theta_social_max,decreasing=TRUE)[1:50]
  idx_social_min <- order(theta_social_min,decreasing=TRUE)[1:50]
  predictors <- data.frame('coefs_material_max'=round(100*theta_material_max[idx_material_max]), 'phens_material_max'=gsub(',',' ',pheId_metadata_firstHalf$Phenotype[idx_material_max]),
                           'coefs_material_min'=round(100*theta_material_min[idx_material_min]), 'phens_material_min'=gsub(',',' ',pheId_metadata_firstHalf$Phenotype[idx_material_min]),
                           'coefs_social_max'=round(100*theta_social_max[idx_social_max]), 'phens_social_max'=gsub(',',' ',pheId_metadata_firstHalf$Phenotype[idx_social_max]),
                           'coefs_social_min'=round(100*theta_social_min[idx_social_min]), 'phens_social_min'=gsub(',',' ',pheId_metadata_firstHalf$Phenotype[idx_social_min]))
  predictors
}
for (i in seq(1,length(theta.mortality.predictors))){
  write.csv(theta.mortality.predictors[[i]][theta.mortality.predictors[[i]][,1]>0,1:2],
            file=paste0('material_mortality_predictor_',i,'_max.csv'),row.names=FALSE)
  write.csv(theta.mortality.predictors[[i]][theta.mortality.predictors[[i]][,3]>0,3:4],
            file=paste0('material_mortality_predictor_',i,'_min.csv'),row.names=FALSE)
  write.csv(theta.mortality.predictors[[i]][theta.mortality.predictors[[i]][,5]>0,5:6],
            file=paste0('social_mortality_predictor_',i,'_max.csv'),row.names=FALSE)
  write.csv(theta.mortality.predictors[[i]][theta.mortality.predictors[[i]][,7]>0,7:8],
            file=paste0('social_mortality_predictor_',i,'_min.csv'),row.names=FALSE)
}


load('SES_admission_results_120220.RData')
load('metadata.RData')
pheId_metadata_firstHalf <- read.csv('pheId_metadata_firstHalf.csv')[,-1]
X <- as.matrix(cbind(log(theta[theta_admit_matches,]+1), log(HU$Total[HU_admit_matches]+1),
                     as.numeric(patients_augmented$sex[pat_admit_matches]=='F')))
ages <- cbind(patients_augmented$age[pat_admit_matches]>=0 & patients_augmented$age[pat_admit_matches]<20,
              patients_augmented$age[pat_admit_matches]>=20 & patients_augmented$age[pat_admit_matches]<40,
              patients_augmented$age[pat_admit_matches]>=40 & patients_augmented$age[pat_admit_matches]<60,
              patients_augmented$age[pat_admit_matches]>=60 & patients_augmented$age[pat_admit_matches]<80,
              patients_augmented$age[pat_admit_matches]>=80)
ses <- as.matrix(patients_augmented[pat_admit_matches,4:5])
Xstar <- cbind(X*ses[,1], X*ses[,2])
Xstar.sds <- apply(Xstar,2,sd)
thetas <- list(results[[1]][-c(1:1332)],results[[3]][-c(1:1332)],results[[5]][-c(1:1332)])
theta.admit.predictors <- foreach(theta=thetas) %do% {
  theta_i <- theta*Xstar.sds
  theta_material_max <- theta_i[1:1326]/max(theta_i[1:1326]);   theta_material_min <- theta_i[1:1326]/min(theta_i[1:1326])
  theta_social_max <- theta_i[1330:2655]/max(theta_i[1330:2655]);   theta_social_min <- theta_i[1330:2655]/min(theta_i[1330:2655])
  idx_material_max <- order(theta_material_max,decreasing=TRUE)[1:50]
  idx_material_min <- order(theta_material_min,decreasing=TRUE)[1:50]
  idx_social_max <- order(theta_social_max,decreasing=TRUE)[1:50]
  idx_social_min <- order(theta_social_min,decreasing=TRUE)[1:50]
  predictors <- data.frame('coefs_material_max'=round(100*theta_material_max[idx_material_max]), 'phens_material_max'=gsub(',',' ',pheId_metadata_firstHalf$Phenotype[idx_material_max]),
                           'coefs_material_min'=round(100*theta_material_min[idx_material_min]), 'phens_material_min'=gsub(',',' ',pheId_metadata_firstHalf$Phenotype[idx_material_min]),
                           'coefs_social_max'=round(100*theta_social_max[idx_social_max]), 'phens_social_max'=gsub(',',' ',pheId_metadata_firstHalf$Phenotype[idx_social_max]),
                           'coefs_social_min'=round(100*theta_social_min[idx_social_min]), 'phens_social_min'=gsub(',',' ',pheId_metadata_firstHalf$Phenotype[idx_social_min]))
  predictors
}
for (i in seq(1,length(theta.admit.predictors))){
  write.csv(theta.admit.predictors[[i]][theta.admit.predictors[[i]][,1]>0,1:2],
            file=paste0('material_admit_predictor_',i,'_max.csv'),row.names=FALSE)
  write.csv(theta.admit.predictors[[i]][theta.admit.predictors[[i]][,3]>0,3:4],
            file=paste0('material_admit_predictor_',i,'_min.csv'),row.names=FALSE)
  write.csv(theta.admit.predictors[[i]][theta.admit.predictors[[i]][,5]>0,5:6],
            file=paste0('social_admit_predictor_',i,'_max.csv'),row.names=FALSE)
  write.csv(theta.admit.predictors[[i]][theta.admit.predictors[[i]][,7]>0,7:8],
            file=paste0('social_admit_predictor_',i,'_min.csv'),row.names=FALSE)
}


load('SES_EDadmit_results_120220.RData')
load('metadata.RData')
pheId_metadata_firstHalf <- read.csv('pheId_metadata_firstHalf.csv')[,-1]
X <- as.matrix(cbind(log(theta[theta_ED_matches,]+1), log(HU$Total[HU_ED_matches]+1),
                     as.numeric(patients_augmented$sex[pat_ED_matches]=='F')))
ages <- cbind(patients_augmented$age[pat_ED_matches]>=0 & patients_augmented$age[pat_ED_matches]<20,
              patients_augmented$age[pat_ED_matches]>=20 & patients_augmented$age[pat_ED_matches]<40,
              patients_augmented$age[pat_ED_matches]>=40 & patients_augmented$age[pat_ED_matches]<60,
              patients_augmented$age[pat_ED_matches]>=60 & patients_augmented$age[pat_ED_matches]<80,
              patients_augmented$age[pat_ED_matches]>=80)
ses <- as.matrix(patients_augmented[pat_ED_matches,4:5])
Xstar <- cbind(X*ses[,1], X*ses[,2])
Xstar.sds <- apply(Xstar,2,sd)
thetas <- list(results[[1]][-c(1:1332)],results[[2]][-c(1:1332)],results[[3]][-c(1:1332)])
theta.ED.predictors <- foreach(theta=thetas) %do% {
  theta_i <- theta*Xstar.sds
  theta_material_max <- theta_i[1:1326]/max(theta_i[1:1326]);   theta_material_min <- theta_i[1:1326]/min(theta_i[1:1326])
  theta_social_max <- theta_i[1330:2655]/max(theta_i[1330:2655]);   theta_social_min <- theta_i[1330:2655]/min(theta_i[1330:2655])
  idx_material_max <- order(theta_material_max,decreasing=TRUE)[1:50]
  idx_material_min <- order(theta_material_min,decreasing=TRUE)[1:50]
  idx_social_max <- order(theta_social_max,decreasing=TRUE)[1:50]
  idx_social_min <- order(theta_social_min,decreasing=TRUE)[1:50]
  predictors <- data.frame('coefs_material_max'=round(100*theta_material_max[idx_material_max]), 'phens_material_max'=gsub(',',' ',pheId_metadata_firstHalf$Phenotype[idx_material_max]),
                           'coefs_material_min'=round(100*theta_material_min[idx_material_min]), 'phens_material_min'=gsub(',',' ',pheId_metadata_firstHalf$Phenotype[idx_material_min]),
                           'coefs_social_max'=round(100*theta_social_max[idx_social_max]), 'phens_social_max'=gsub(',',' ',pheId_metadata_firstHalf$Phenotype[idx_social_max]),
                           'coefs_social_min'=round(100*theta_social_min[idx_social_min]), 'phens_social_min'=gsub(',',' ',pheId_metadata_firstHalf$Phenotype[idx_social_min]))
  predictors
}
for (i in seq(1,length(theta.ED.predictors))){
  write.csv(theta.ED.predictors[[i]][theta.ED.predictors[[i]][,1]>0,1:2],
            file=paste0('material_ED_predictor_',i,'_max.csv'),row.names=FALSE)
  write.csv(theta.ED.predictors[[i]][theta.ED.predictors[[i]][,3]>0,3:4],
            file=paste0('material_ED_predictor_',i,'_min.csv'),row.names=FALSE)
  write.csv(theta.ED.predictors[[i]][theta.ED.predictors[[i]][,5]>0,5:6],
            file=paste0('social_ED_predictor_',i,'_max.csv'),row.names=FALSE)
  write.csv(theta.ED.predictors[[i]][theta.ED.predictors[[i]][,7]>0,7:8],
            file=paste0('social_ED_predictor_',i,'_min.csv'),row.names=FALSE)
}


load('SES_drugs_results_120220.RData')
load('metadata.RData')
theta <- read.table('MTL_mixEHR_sureLDA_predictions_withactes_noStopWords_firstHalf.txt',sep=',',header=FALSE)
pheId_metadata_firstHalf <- read.csv('pheId_metadata_firstHalf.csv')[,-1]
X <- as.matrix(cbind(log(theta[theta_drug_matches,]+1), log(HU$Total[HU_drug_matches]+1),
                     as.numeric(patients_augmented$sex[pat_drug_matches]=='F')))
ages <- cbind(patients_augmented$age[pat_drug_matches]>=0 & patients_augmented$age[pat_drug_matches]<20,
              patients_augmented$age[pat_drug_matches]>=20 & patients_augmented$age[pat_drug_matches]<40,
              patients_augmented$age[pat_drug_matches]>=40 & patients_augmented$age[pat_drug_matches]<60,
              patients_augmented$age[pat_drug_matches]>=60 & patients_augmented$age[pat_drug_matches]<80,
              patients_augmented$age[pat_drug_matches]>=80)
ses <- as.matrix(patients_augmented[pat_drug_matches,4:5])
Xstar <- cbind(X*ses[,1], X*ses[,2])
Xstar.sds <- apply(Xstar,2,sd)
thetas <- list(results[[1]][-c(1:1332)],results[[2]][-c(1:1332)],results[[3]][-c(1:1332)])
theta.drug.predictors <- foreach(theta=thetas) %do% {
  theta_i <- theta*Xstar.sds
  theta_material_max <- theta_i[1:1326]/max(theta_i[1:1326]);   theta_material_min <- theta_i[1:1326]/min(theta_i[1:1326])
  theta_social_max <- theta_i[1330:2655]/max(theta_i[1330:2655]);   theta_social_min <- theta_i[1330:2655]/min(theta_i[1330:2655])
  idx_material_max <- order(theta_material_max,decreasing=TRUE)[1:50]
  idx_material_min <- order(theta_material_min,decreasing=TRUE)[1:50]
  idx_social_max <- order(theta_social_max,decreasing=TRUE)[1:50]
  idx_social_min <- order(theta_social_min,decreasing=TRUE)[1:50]
  predictors <- data.frame('coefs_material_max'=round(100*theta_material_max[idx_material_max]), 'phens_material_max'=gsub(',',' ',pheId_metadata_firstHalf$Phenotype[idx_material_max]),
                           'coefs_material_min'=round(100*theta_material_min[idx_material_min]), 'phens_material_min'=gsub(',',' ',pheId_metadata_firstHalf$Phenotype[idx_material_min]),
                           'coefs_social_max'=round(100*theta_social_max[idx_social_max]), 'phens_social_max'=gsub(',',' ',pheId_metadata_firstHalf$Phenotype[idx_social_max]),
                           'coefs_social_min'=round(100*theta_social_min[idx_social_min]), 'phens_social_min'=gsub(',',' ',pheId_metadata_firstHalf$Phenotype[idx_social_min]))
  predictors
}
for (i in seq(1,length(theta.drug.predictors))){
  write.csv(theta.drug.predictors[[i]][theta.drug.predictors[[i]][,1]>0,1:2],
            file=paste0('material_drug_predictor_',i,'_max.csv'),row.names=FALSE)
  write.csv(theta.drug.predictors[[i]][theta.drug.predictors[[i]][,3]>0,3:4],
            file=paste0('material_drug_predictor_',i,'_min.csv'),row.names=FALSE)
  write.csv(theta.drug.predictors[[i]][theta.drug.predictors[[i]][,5]>0,5:6],
            file=paste0('social_drug_predictor_',i,'_max.csv'),row.names=FALSE)
  write.csv(theta.drug.predictors[[i]][theta.drug.predictors[[i]][,7]>0,7:8],
            file=paste0('social_drug_predictor_',i,'_min.csv'),row.names=FALSE)
}





theta.material.mortality_0020 <- glm(X,mortality$dead[mort_matches],family='binomial',
                                     weights=I(patients_augmented$age[pat_mort_matches]>=0 & patients_augmented$age[pat_mort_matches]<20 & patients_augmented$material_deprivation[pat_mort_matches]<0))
theta.material_low.mortality_2040 <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',
                                               weights=I(patients_augmented$age[pat_mort_matches]>=20 & patients_augmented$age[pat_mort_matches]<40 & patients_augmented$material_deprivation[pat_mort_matches]<0))
theta.material_low.mortality_4060 <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',
                                               weights=I(patients_augmented$age[pat_mort_matches]>=40 & patients_augmented$age[pat_mort_matches]<60 & patients_augmented$material_deprivation[pat_mort_matches]<0))
theta.material_low.mortality_6080 <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',
                                               weights=I(patients_augmented$age[pat_mort_matches]>=60 & patients_augmented$age[pat_mort_matches]<80 & patients_augmented$material_deprivation[pat_mort_matches]<0))
theta.material_low.mortality_8000 <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',
                                               weights=I(patients_augmented$age[pat_mort_matches]>=80 & patients_augmented$material_deprivation[pat_mort_matches]<0))
theta.material_high.mortality_0020 <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',
                                                weights=I(patients_augmented$age[pat_mort_matches]>=0 & patients_augmented$age[pat_mort_matches]<20 & patients_augmented$material_deprivation[pat_mort_matches]>0))
theta.material_high.mortality_2040 <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',
                                                weights=I(patients_augmented$age[pat_mort_matches]>=20 & patients_augmented$age[pat_mort_matches]<40 & patients_augmented$material_deprivation[pat_mort_matches]>0))
theta.material_high.mortality_4060 <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',
                                                weights=I(patients_augmented$age[pat_mort_matches]>=40 & patients_augmented$age[pat_mort_matches]<60 & patients_augmented$material_deprivation[pat_mort_matches]>0))
theta.material_high.mortality_6080 <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',
                                                weights=I(patients_augmented$age[pat_mort_matches]>=60 & patients_augmented$age[pat_mort_matches]<80 & patients_augmented$material_deprivation[pat_mort_matches]>0))
theta.material_high.mortality_8000 <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',
                                                weights=I(patients_augmented$age[pat_mort_matches]>=80 & patients_augmented$material_deprivation[pat_mort_matches]>0))
theta.social_low.mortality_0020 <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',
                                             weights=I(patients_augmented$age[pat_mort_matches]>=0 & patients_augmented$age[pat_mort_matches]<20 & patients_augmented$social_deprivation[pat_mort_matches]<0))
theta.social_low.mortality_2040 <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',
                                             weights=I(patients_augmented$age[pat_mort_matches]>=20 & patients_augmented$age[pat_mort_matches]<40 & patients_augmented$social_deprivation[pat_mort_matches]<0))
theta.social_low.mortality_4060 <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',
                                             weights=I(patients_augmented$age[pat_mort_matches]>=40 & patients_augmented$age[pat_mort_matches]<60 & patients_augmented$social_deprivation[pat_mort_matches]<0))
theta.social_low.mortality_6080 <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',
                                             weights=I(patients_augmented$age[pat_mort_matches]>=60 & patients_augmented$age[pat_mort_matches]<80 & patients_augmented$social_deprivation[pat_mort_matches]<0))
theta.social_low.mortality_8000 <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',
                                             weights=I(patients_augmented$age[pat_mort_matches]>=80 & patients_augmented$social_deprivation[pat_mort_matches]<0))
theta.social_high.mortality_0020 <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',
                                              weights=I(patients_augmented$age[pat_mort_matches]>=0 & patients_augmented$age[pat_mort_matches]<20 & patients_augmented$social_deprivation[pat_mort_matches]>0))
theta.social_high.mortality_2040 <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',
                                              weights=I(patients_augmented$age[pat_mort_matches]>=20 & patients_augmented$age[pat_mort_matches]<40 & patients_augmented$social_deprivation[pat_mort_matches]>0))
theta.social_high.mortality_4060 <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',
                                              weights=I(patients_augmented$age[pat_mort_matches]>=40 & patients_augmented$age[pat_mort_matches]<60 & patients_augmented$social_deprivation[pat_mort_matches]>0))
theta.social_high.mortality_6080 <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',
                                              weights=I(patients_augmented$age[pat_mort_matches]>=60 & patients_augmented$age[pat_mort_matches]<80 & patients_augmented$social_deprivation[pat_mort_matches]>0))
theta.social_high.mortality_8000 <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',
                                              weights=I(patients_augmented$age[pat_mort_matches]>=80 & patients_augmented$social_deprivation[pat_mort_matches]>0))



# LASSO models
library(glmnet)

load('patient_demographics.RData')
admissions <- read.csv('admissions_0708.csv')
EDadmits <- read.csv('EDadmits_0708.csv')
mortality <- read.csv('mortality.csv')
HU <- read.csv('HU.csv')
IDs <- unlist(read.csv('MTL_mixEHR_sureLDA_predictions_withactes_noStopWords_firstHalf_patId.csv',header=FALSE))
theta <- read.table('MTL_mixEHR_sureLDA_predictions_withactes_noStopWords_firstHalf.txt',sep=',',header=FALSE)
metaphe_500 <- read.csv("combinedData_firstHalf_combinedData_firstHalf_JCVB0_nmar_K500_iter200_metaphe.csv",header=FALSE)
metapheIDs <- unlist(read.csv('MTL_mixEHR_sureLDA_IDs_patId.csv',header=FALSE))
mortIDs <- intersect(intersect(intersect(IDs,mortality$patId),HU$patId),patients_augmented$id)
admitIDs <- intersect(intersect(intersect(IDs,admissions$patId),HU$patId),patients_augmented$id)
patIDs <- intersect(intersect(IDs,HU$patId),patients_augmented$id)

theta_admit_matches <- which(IDs %in% admitIDs)
HU_admit_matches <- which(HU$patId %in% admitIDs)
pat_admit_matches <- which(patients_augmented$id %in% admitIDs)
theta_mort_matches <- which(IDs %in% mortIDs)
HU_mort_matches <- which(HU$patId %in% mortIDs)
pat_mort_matches <- which(patients_augmented$id %in% mortIDs)
theta_matches <- which(IDs %in% patIDs)
HU_matches <- which(HU$patId %in% patIDs)
pat_matches <- which(patients_augmented$id %in% patIDs)
mort_matches <- which(mortality$patId %in%patIDs)
admit_matches <- which(admissions$patId %in%patIDs)
metaphe_admit_matches <- which(metapheIDs %in% admitIDs)
metaphe_mort_matches <- which(metapheIDs %in% mortIDs)

X <- as.matrix(cbind(log(theta[theta_mort_matches,]+1), log(HU$Total[HU_mort_matches]+1), patients_augmented$age[pat_mort_matches],
                     patients_augmented$material_deprivation[pat_mort_matches], patients_augmented$social_deprivation[pat_mort_matches]))
theta.mortality.mixture <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',type.measure='auc')
theta.mortality.mixture.AUC <- max(theta.mortality.mixture$cvm)

X <- as.matrix(cbind(log(theta[theta_admit_matches,]+1), log(HU$Total[HU_admit_matches]+1), patients_augmented$age[pat_admit_matches],
                     patients_augmented$material_deprivation[pat_admit_matches], patients_augmented$social_deprivation[pat_admit_matches]))
theta.admit.mixture <- cv.glmnet(X,admissions$admission[admit_matches],family='binomial',type.measure='auc')
theta.admit.mixture.AUC <- max(theta.admit.mixture$cvm)

X <- as.matrix(cbind(log(metaphe_500[theta_mort_matches,]+1), log(HU$Total[HU_mort_matches]+1), patients_augmented$age[pat_mort_matches],
                     patients_augmented$material_deprivation[pat_mort_matches], patients_augmented$social_deprivation[pat_mort_matches]))
metaphe_500.mortality.mixture <- cv.glmnet(X,mortality$dead[mort_matches],family='binomial',type.measure='auc')
metaphe_500.mortality.mixture.AUC <- max(metaphe_500.mortality.mixture$cvm)

X <- as.matrix(cbind(log(metaphe_500[theta_admit_matches,]+1), log(HU$Total[HU_admit_matches]+1), patients_augmented$age[pat_admit_matches],
                     patients_augmented$material_deprivation[pat_admit_matches], patients_augmented$social_deprivation[pat_admit_matches]))
metaphe_500.admit.mixture <- cv.glmnet(X,admissions$admission[admit_matches],family='binomial',type.measure='auc')
metaphe_500.admit.mixture.AUC <- max(metaphe_500.admit.mixture$cvm)

results <- list('theta.mortality.mixture'=theta.mortality.mixture.AUC,
                'theta.admit.mixture'=theta.admit.mixture.AUC,
                'metaphe_500.mortality.mixture'=metaphe_500.mortality.mixture.AUC,
                'metaphe_500.admit.mixture'=metaphe_500.admit.mixture.AUC)
save(results,file='theta_prediction_aucs_112720.RData')


# MAP compressed to uncompressed
IDs <- unique(map_prior[,1])
map <- sapply(1:length(IDs),function(i){
  if (i %% 10000 == 0){print(paste(i,'/',length(IDs)))}
  id <- IDs[i]
  matches <- which(map_prior[,1]==id)
  out <- rep(0,1326)
  out[map_prior[matches,2]+1] <- map_prior[matches,3]
  out
})


 
# Train mixEHR-sureLDA with different topics for different socioeconomic clusters

map_prior <- read.table('map_prior_compressed.txt',sep=' ',header=FALSE)
colnames(map_prior) <- c('id','pheid','prob')
map_prior <- merge(map_prior,deprivation[,c('id','cluster')],by='id')
map_prior$pheid <- map_prior$pheid + 1326*(map_prior$cluster-1)
write.table(map_prior[,1:3],file='map_prior_compressed_SES.txt',row.names=FALSE,col.names=FALSE)
  
  

