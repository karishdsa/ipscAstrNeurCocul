#Running the scaden_deconvolve.sh script

#using combined samples single file with raw counts - scdata
#-----------------------------------------------------------
#with the raw bulk counts

cd ~/ipscAstrocyteNeuron/ipscAstroNeu/bash/

input_path_sc_data="/home/ipscAstrocyteNeuron/results/deconv_scaden/input/sc_allsamples_rawcounts/"
bulk_data_file="/home/ipscAstrocyteNeuron/results/deconv_scaden/input/bulk_data/ipsc_bulk_raw.txt"

intermediary_files_path="/home/ipscAstrocyteNeuron/results/deconv_scaden/allsamples_run/"
prediction_path="/home/ipscAstrocyteNeuron/results/deconv_scaden/allsamples_run/predictions/"

mkdir $prediction_path

opt="bulkRaw_scCombRaw"
var=0.1
nCells=100
culture="comb"

#nCells_list=(100 500 1000 2000 3000 4000)
#nCells_list=(4500 5000)


for nCells in ${nCells_list[@]};
 do
   newvarlog=`echo $var |  sed 's/\./_/'`
  log_file="/home/ipscAstrocyteNeuron/ipscAstroNeu/logs/scaden_opt_"$opt"_var"$newvarlog"_cells"$nCells".log"
  echo "#./scaden_deconv.sh "$input_path_sc_data $bulk_data_file $intermediary_files_path $prediction_path $var $nCells $culture" > "$log_file" 2>&1"
done

#####################################################################################################
#----------------------------------------------------------------------------------------------------

#using combined samples single file with raw counts - scdata - varying the var_cut off too ( above runs were done varying the nCells at var 0.1)
#------------------------------------------------------------

#with the raw bulk counts
#------------------------
#***** vary nCells var cutoff

cd ~/ipscAstrocyteNeuron/ipscAstroNeu/bash/

input_path_sc_data="/home/ipscAstrocyteNeuron/results/deconv_scaden/input/sc_allsamples_rawcounts/"
bulk_data_file="/home/ipscAstrocyteNeuron/results/deconv_scaden/input/bulk_data/ipsc_bulk_raw.txt"

intermediary_files_path="/home/ipscAstrocyteNeuron/results/deconv_scaden/allsamples_run/"

prediction_path="/home/ipscAstrocyteNeuron/results/deconv_scaden/allsamples_run/predictions/"

mkdir $prediction_path

opt="bulkRaw_scCombRaw"

culture="comb"

nCells_list=(3000) #(5000)

var_cutoff_list=(0.2 0.3)

for nCells in ${nCells_list[@]};
 do
  
  for var in ${var_cutoff_list[@]};
   do
    newvarlog=`echo $var |  sed 's/\./_/'`
    log_file="/home/ipscAstrocyteNeuron/ipscAstroNeu/logs/scaden_opt_"$opt"_var"$newvarlog"_cells"$nCells".log"
    echo "#./scaden_deconv.sh "$input_path_sc_data $bulk_data_file $intermediary_files_path $prediction_path $var $nCells $culture" > "$log_file" 2>&1"
  done
done


#******varying the nSamples

input_path_sc_data="/home/ipscAstrocyteNeuron/results/deconv_scaden/input/sc_allsamples_rawcounts/"
bulk_data_file="/home/ipscAstrocyteNeuron/results/deconv_scaden/input/bulk_data/ipsc_bulk_raw.txt"

intermediary_files_path="/home/ipscAstrocyteNeuron/results/deconv_scaden/allsamples_run/"

prediction_path="/home/ipscAstrocyteNeuron/results/deconv_scaden/allsamples_run/predictions_nsamples/"

mkdir $prediction_path

opt="bulkRaw_scCombRaw"

culture="comb"

var=0.1
#nCells_list=(5000)
nCells_list=(3000)
nSamples_list=(8000 10000) #(2000 3000 5000)
for nCells in ${nCells_list[@]};
 do
  
  for nSamples in ${nSamples_list[@]};
   do
    newvarlog=`echo $var |  sed 's/\./_/'`
    log_file="/home/ipscAstrocyteNeuron/ipscAstroNeu/logs/scaden_opt_"$opt"_var"$newvarlog"_cells"$nCells"_samples"$nSamples".log"
    echo "#./scaden_deconv.sh "$input_path_sc_data $bulk_data_file $intermediary_files_path $prediction_path $var $nCells $culture $nSamples" > "$log_file" 2>&1"
  done
done


