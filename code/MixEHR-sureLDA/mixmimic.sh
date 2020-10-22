#!/bin/bash

ehrdata=./mixmimic/mimic_trainData.txt
ehrmeta=./mixmimic/mimic_meta.txt
ehrprior=./mixmimic/MIMIC_prior_phecodes1d.txt

K=1515
# K=50

niter=100

# mar model (evaluation purpose only)
# ./mixehr --outputIntermediates -f $ehrdata -m $ehrmeta -k $K -i $niter --inferenceMethod JCVB0 --mar

# nmar model
./mixehr --outputIntermediates -f $ehrdata -m $ehrmeta -trp $ehrprior -k $K -i $niter --inferenceMethod JCVB0 --maxcores 10

# stochastic variational inference for massive scale EHR
# bathsize=5000
# stepsize=0.6
# burnin=10
# ./mixehr --outputIntermediates -f $ehrdata -m $ehrmeta -b $batchsize -k $K -i $niter -n SJCVB0 -r $burnin -s $stepsize

