#Script runs all the steps of the scaden package - deconvolution

##Usage
#input_path_sc_data="/home/ipscAstrocyteNeuron/results/deconv_scaden/input/"
#intermediary_files_path="/home/ipscAstrocyteNeuron/results/deconv_scaden/optimisation/"

#prediction_path="/home/ipscAstrocyteNeuron/results/deconv_scaden/optimisation/predictions/"

#./scaden_deconv.sh $input_path_sc_data $bulk_data_file $intermediary_files_path $prediction_path $var_cutoff $nCells $culture $nSamples

#if the last parameter is not provided, it is set to 1000 - scaden's default

#!/bin/bash

sc_path=$1
bulk_data_file=$2
inter_path=$3
pred_path=$4
var_cutoff=$5
nCells=$6
culture=$7 #the sc samples used to create the training data
nSamples=$8

if [ -z "${nSamples}" ]; then
    echo "nSamples is set to default 1000"
    nSamples=1000 #default
fi


echo "Single cell data : " $sc_path
echo "Bulk data file: " $sc_path
echo "Intermediary files to be saved in : " $inter_path
echo "Prediction output file in : " $pred_path
echo "Variance cutoff : " $var_cutoff
echo "No. of cells : " $nCells
echo "Culture : " $culture
echo "No. of samples : " $nSamples

#Make the directories
newvar_cutoff=`echo $var_cutoff |  sed 's/\./_/'`

#opt_path=$inter_path"/var"$newvar_cutoff"_cells"$nCells"/"
opt_path=$inter_path"/var"$newvar_cutoff"_cells"$nCells"_samples"$nSamples"/"
sim_path=$opt_path"simulated_samples/"
model_path=$opt_path"model/"

#pred_file=$pred_path"/scaden_pred_"$culture"_var"$newvar_cutoff"_cells"$nCells".txt"
pred_file=$pred_path"/scaden_pred_"$culture"_var"$newvar_cutoff"_cells"$nCells"_samples"$nSamples".txt"


echo "Creating directories for intermediate file :"
echo $opt_path
echo $sim_path
echo $model_path
echo $pred_file

mkdir -p $sim_path
mkdir $model_path


 
cd $opt_path

#---Generating training data----------------------------------------
# bulk simulation
scaden simulate --cells $nCells --n_samples $nSamples --data $sc_path --out $sim_path

#---pre-processing of data----------------------------------------
#Create a new file for training which only contains the intersection of genes between the training and the prediction data. #Also the training data will be log2-transformed and scaled to the range [0,1].
scaden process ./data.h5ad $bulk_data_file --processed_path ./processed.h5ad --var_cutoff $var_cutoff

#---Training a Scaden ensemble model.----------------------------------------
scaden train ./processed.h5ad --model_dir $model_path --steps 5000 #default 5000 steps recommended

#---Prediction----------------------------------------
scaden predict --model_dir $model_path --outname $pred_file $bulk_data_file

