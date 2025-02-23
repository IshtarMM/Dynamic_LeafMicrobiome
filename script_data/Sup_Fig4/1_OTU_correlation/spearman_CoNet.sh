#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=03:00:00
#PBS -N CoNet5_25
cd $PBS_O_WORKDIR

export LIB=$PWD/lib
export CLASSPATH=${CLASSPATH}:${LIB}/CoNet.jar
inputfile=input_otutable.txt

java be.ac.vub.bsb.cooccurrence.cmd.CooccurrenceAnalyser --method ensemble --input $inputfile --max 3 --format gdl --ensemblemethods correl_pearson --matrixtype abundance --thresholdguessing quantile --guessingparam 0.05 --resamplemethod shuffle_rows --multigraph --verbosity FATAL --minetdisc equalfreq --minetmiestimator mi.shrink --kernelwidth 0.25 --edgethreshold 0.05 --nantreatment none --nantreatmentparam 1 --networkmergestrategy union --iterations 100 --correlnonrandp fisherz --multicorr bonferroni --inference mrnet --min 1 --measure2 supp --measure1 conf --minsupport 2 --topbottom --scoremergestrategy scoremethodnum --filter row_minocc --filterparameter 1.0  --output thresholds.txt
java be.ac.vub.bsb.cooccurrence.cmd.CooccurrenceAnalyser --method ensemble --input $inputfile --max 3 --format gdl --ensemblemethods correl_pearson --matrixtype abundance --randroutine none --resamplemethod shuffle_rows --multigraph --verbosity FATAL --minetdisc equalfreq --minetmiestimator mi.shrink --kernelwidth 0.25 --edgethreshold 0.05 --nantreatment none --nantreatmentparam 1 --networkmergestrategy union --iterations 100 --correlnonrandp fisherz --multicorr bonferroni --inference mrnet --min 1 --measure2 supp --measure1 conf --minsupport 2 --scoremergestrategy scoremethodnum --filter row_minocc --filterparameter 1.0 --output all.gdl --ensembleparamfile thresholds.txt
