library(foreach)
library(doParallel)
library(Rcpp)
library(RcppArmadillo)


## Preprocess ICD codes

allCodes <- read.csv('/home/yueli/data/services.csv')
icd_metadata <- as.character(read.csv('/home/yueli/data/icd9_qc.csv')[-c(1:2),1])
icd_metadata <- as.data.frame(cbind(seq(length(icd_metadata)), icd_metadata)); colnames(icd_metadata) <- c('pheId','icd')
allCodes$icd <- as.character(allCodes$icd); icd_metadata$icd <- as.character(icd_metadata$icd)
icds_sorted <- sort(unique(allCodes$icd))[-1]
allCodes$codIdx <- sapply(1:nrow(allCodes),function(i){
  match <- which(icd_metadata$icd==allCodes$icd[i])
  if (length(match)==0){NA}
  else{match}
})
icds_mixehr_surelda <- allCodes[!is.na(allCodes[,'codIdx']),c('id','date','codIdx')]
colnames(icds_mixehr_surelda) <- c('patId','date','pheId')


## Preprocess drug codes

drugs <- read.csv('/home/yueli/data/drugs.csv')
cod_list <- sort(as.numeric(unique(drugs[,4])))
drugs$codIdx <- sapply(1:nrow(drugs),function(i){
  match <- which(cod_list==drugs$din[i])
  if (length(match)==0){NA}
  else{match}
})
drugs_mixehr_surelda <- drugs[,c('id','date','codIdx')]
colnames(drugs_mixehr_surelda) <- c('patId','date','pheId')
drug_metadata <- read.csv('/home/yueli/data/drug_product.csv')[,c('drug_code','drug_identification_number','brand_name')]
drug_metadata$drug_identification_number <- as.numeric(as.character(drug_metadata$drug_identification_number))
drug_metadata$pheId <- sapply(drug_metadata$drug_identification_number,function(id){
  if (id %in% cod_list){which(cod_list==id)}
  else{NA}
})
drug_metadata <- drug_metadata[!is.na(drug_metadata$pheId),c('pheId','drug_code','drug_identification_number','brand_name')]
drug_metadata <- drug_metadata[order(drug_metadata$pheId),]


## Preprocess ACT codes and combine all three data types into combined data framers

actes_metadata <- read.csv('/home/yueli/data/actes.csv')[,2:3]
colnames(actes_metadata) <- c('pheId','description')
actes_mixehr_surelda <- data.frame('patId'=icds$id,'date'=icds$date,'pheId'=icds$act_code,
                                   'typeId'=3,'state'=1)

metadata <- list('ICD9'=icd_metadata, 'Drug'=drug_metadata, 'Act'=actes_metadata)
save(metadata,file='metadata.RData')
icds_mixehr_surelda$typeId <- 1; drugs_mixehr_surelda$typeId <- 2
icds_mixehr_surelda$state <- drugs_mixehr_surelda$state <- 1
combined_data <- rbind(icds_mixehr_surelda, drugs_mixehr_surelda, actes_mixehr_surelda)
combined_data <- combined_data[order(combined_data$pheId),]
combined_data <- combined_data[order(combined_data$typeId),]
combined_data <- combined_data[order(combined_data$patId),]
save(combined_data,file='combined_data.RData')


## Only take data before 1/1/2007

firstHalf <- as.Date(combined_data$date) - as.Date('2007-01-01') < 0
combined_data_firstHalf <- combined_data[firstHalf,]
combined_data_firstHalf <- combined_data_firstHalf[order(combined_data_firstHalf$pheId),]
combined_data_firstHalf <- combined_data_firstHalf[order(combined_data_firstHalf$typeId),]
combined_data_firstHalf <- combined_data_firstHalf[order(combined_data_firstHalf$patId),]
save(combined_data_firstHalf,file='combined_data_firstHalf.RData')


## Compile unique cases (get rid of duplicates)

duplicates <- duplicated(combined_data[,c('patId','typeId','pheId')])
uniqueCases <- combined_data[!duplicates,c('patId','typeId','pheId','state')]

duplicates <- duplicated(combined_data_firstHalf[,c('patId','typeId','pheId')])
uniqueCases_firstHalf <- combined_data_firstHalf[!duplicates,c('patId','typeId','pheId','state')]
uniqueCases_firstHalf <- uniqueCases_firstHalf[order(uniqueCases_firstHalf$pheId),]
uniqueCases_firstHalf <- uniqueCases_firstHalf[order(uniqueCases_firstHalf$typeId),]
uniqueCases_firstHalf <- uniqueCases_firstHalf[order(uniqueCases_firstHalf$patId),]


## Produce fully preprocessed data and metadata files for input into mixEHR-sureLDA
# Uses C++ script 

sourceCpp('create_datamat_david.cpp')

freqs <- findCounts(as.matrix(uniqueCases[,1:3]),as.matrix(combined_data[,c(1,4,3)]))
uniqueCases$freq <- freqs
save(uniqueCases,file='combinedData_processed.RData')
write.table(uniqueCases,file='combinedData_processed.txt',sep=' ',row.names=FALSE,col.names=FALSE)

freqs <- findCounts(as.matrix(uniqueCases_firstHalf[,1:3]),as.matrix(combined_data_firstHalf[,c(1,4,3)]))
uniqueCases_firstHalf$freq <- freqs
save(uniqueCases_firstHalf,file='combinedData_firstHalf_processed.RData')
write.table(uniqueCases_firstHalf,file='combinedData_firstHalf_processed.txt',sep=' ',row.names=FALSE,col.names=FALSE)

combined_metadata <- uniqueCases[,c('typeId','pheId')]
combined_metadata <- combined_metadata[order(combined_metadata$pheId),]
combined_metadata <- combined_metadata[order(combined_metadata$typeId),]
combined_metadata <- combined_metadata[!duplicated(combined_metadata),]
combined_metadata$stateCnt <- 1
write.table(combined_metadata,file='combinedMetaData_processed.txt',sep=' ',row.names=FALSE,col.names=FALSE)

combined_metadata_firstHalf <- uniqueCases_firstHalf[,c('typeId','pheId')]
combined_metadata_firstHalf <- combined_metadata_firstHalf[order(combined_metadata_firstHalf$pheId),]
combined_metadata_firstHalf <- combined_metadata_firstHalf[order(combined_metadata_firstHalf$typeId),]
combined_metadata_firstHalf <- combined_metadata_firstHalf[!duplicated(combined_metadata_firstHalf),]
combined_metadata_firstHalf$stateCnt <- 1
write.table(combined_metadata_firstHalf,file='combinedMetaData_firstHalf_processed.txt',sep=' ',row.names=FALSE,col.names=FALSE)


## MAKE SURE TO REMOVE PATIENTS FROM DATA WITH ALL PRIORS=0; WILL RESULT IN ERRORS


## Filter out very common codes (prevalence >= 0.25)

uniqueCodes <- uniqueCases[!duplicated(uniqueCases[,c('typeId','pheId')]), c('typeId','pheId')]
nPats <- length(unique(uniqueCases$patId))
codeCounts <- table(uniqueCases[,c('typeId','pheId')])
codeFreqs <- codeCounts/nPats

stopWords <- as.data.frame(t(sapply(which(codeFreqs>=0.25),function(x){
  i <- (x-1)%%3 + 1
  j <- ceiling(x/3)
  as.integer(c(rownames(codeFreqs)[i], colnames(codeFreqs)[j]))
})))
colnames(stopWords) <- c('typeId','pheId')

uniqueCases_stopIdx <- which((1e7*uniqueCases$typeId + uniqueCases$pheId) %in% (1e7*stopWords$typeId + stopWords$pheId))
uniqueCases <- uniqueCases[-uniqueCases_stopIdx,]
uniqueCases <- uniqueCases[uniqueCases$patId %in% map_prior_compressed[,1],]
save(uniqueCases,file='combinedData_processed_noStopWords.RData')
write.table(uniqueCases,file='combinedData_processed_noStopWords.txt',sep=' ',row.names=FALSE,col.names=FALSE)

uniqueCases_firstHalf_stopIdx <- which((1000000*uniqueCases_firstHalf$typeId + uniqueCases_firstHalf$pheId) %in% (1000000*stopWords$typeId + stopWords$pheId))
uniqueCases_firstHalf <- uniqueCases_firstHalf[-uniqueCases_firstHalf_stopIdx,]
uniqueCases_firstHalf <- uniqueCases_firstHalf[uniqueCases_firstHalf$patId %in% map_prior_firstHalf_compressed[,1],]
save(uniqueCases_firstHalf,file='combinedData_firstHalf_processed_noStopWords.RData')
write.table(uniqueCases_firstHalf,file='combinedData_firstHalf_processed_noStopWords.txt',sep=' ',row.names=FALSE,col.names=FALSE)

combined_metadata <- uniqueCases[,c('typeId','pheId')]
combined_metadata <- combined_metadata[order(combined_metadata$pheId),]
combined_metadata <- combined_metadata[order(combined_metadata$typeId),]
combined_metadata <- combined_metadata[!duplicated(combined_metadata),]
combined_metadata$stateCnt <- 1
write.table(combined_metadata,file='combinedMetaData_processed_noStopWords.txt',sep=' ',row.names=FALSE,col.names=FALSE)

combined_metadata_firstHalf <- uniqueCases_firstHalf[,c('typeId','pheId')]
combined_metadata_firstHalf <- combined_metadata_firstHalf[order(combined_metadata_firstHalf$pheId),]
combined_metadata_firstHalf <- combined_metadata_firstHalf[order(combined_metadata_firstHalf$typeId),]
combined_metadata_firstHalf <- combined_metadata_firstHalf[!duplicated(combined_metadata_firstHalf),]
combined_metadata_firstHalf$stateCnt <- 1
write.table(combined_metadata_firstHalf,file='combinedMetaData_firstHalf_processed_noStopWords.txt',sep=' ',row.names=FALSE,col.names=FALSE)
