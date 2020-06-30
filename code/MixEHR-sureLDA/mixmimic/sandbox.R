setwd('~/Desktop/mixehr-master/mixmimic')

mimic_trainData = read.table('mimic_trainData.txt',sep=' ')
mimic_pats = unique(mimic_trainData[,1])
n = length(mimic_pats); k = 100
mimic_prior = cbind(mimic_pats,
                    matrix(rbinom(n*k,1,0.2)*rbeta(n*k,1,4),n,k))
write.table(mimic_prior,file="mimic_prior.txt",row.names=FALSE,col.names=FALSE)
