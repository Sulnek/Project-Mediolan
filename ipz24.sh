#!/bin/bash

FLAG_FILE="initialization_flagfile"
default_file="eToxPred/tcm600_nr.smi"

# Function to display usage
usage() {
    echo "Usage: $0 [--clean] [--help|-h] [input_file]"
    echo "  --clean       Clean up created files and folders"
    echo "  --help, -h    Display this help message"
    echo "  input_file    Specify the input file for toxicity prediction"
}

# Check if the flag file exists
if [ ! -f "$FLAG_FILE" ]; then
    echo "Running initialization part..."
    # Create a virtual environment with Python 3.7
    python3.8 -m venv venv_etoxpred
    source venv_etoxpred/bin/activate
    pip install -r requirements_etoxpred.txt
    # Download toxity prediction library
    git clone https://github.com/pulimeng/eToxPred.git
    # Get library for ligand generation
    git clone https://github.com/LaMavia/mlf.git
    cd eToxPred
    tar -xzvf etoxpred_best_model.tar.gz
    cd ..

    # Create the flag file to indicate initialization has been done
    touch "$FLAG_FILE"
else
    echo "Initialization already done. Skipping..."
fi

# Check for arguments
if [[ "$1" == "--help" || "$1" == "-h" ]]; then
    usage
    exit 0
elif [[ "$1" == "--clean" ]]; then
    echo "Cleaning up..."
    rm -rf venv_etoxpred eToxPred "$FLAG_FILE"
    echo "Cleanup done."
    exit 0
elif [ $# -eq 1 ] && [ -f "$1" ]; then
    input_file=$1
    mkdir -p eToxPred
else
    input_file=$default_file
fi

echo "Running toxicity prediction..."

source venv_etoxpred/bin/activate
echo "Using file: $input_file"
python lipinski.py "$input_file" -o eToxPred/result_lipinski.csv 
cd eToxPred
python etoxpred_predict.py --datafile result_lipinski.csv --modelfile etoxpred_best_model.joblib --outputfile ../results.csv
cd ..

cut -d ',' -f 2 results.csv > smiles_list.txt

INPUT_FILE="smiles_list.txt"
REC="5zhl_ligandless.pdb"

cd $(dirname $0)

function smiles_to_sdf_file() {
  local SMILES=$1
  local F=$2

  echo "$SMILES" | obabel -i smiles -o sdf -O "$F" --gen2D
}

if test \( ! -f "$REC" \); then
  echo "MISSING PDB[QT]"
  exit 1
fi

if test \( ! -f "$INPUT_FILE" \); then
  echo "MISSING input file"
  exit 1
fi

while IFS= read -r smiles
do
  temp_lig=$(mktemp tmp.XXXXXXXX.sdf)

  smiles_to_sdf_file "$smiles" "$temp_lig"
  GNINA_ARGS=(
    -r "/scr/$REC"
    -l "/scr/$temp_lig"
    # --flex "$REC_FLEX" 
    --autobox_ligand "/scr/5zhl_B_9D0.sdf" # absolutnie było coś innego co innego
    # --flexdist_ligand "$LIG" 
    # --flexdist 3.5
    -o "/scr/whole_docked$smiles.sdf.gz"
    --exhaustiveness 64
  )

  sudo docker run -v $(pwd):/scr gnina/gnina gnina "${GNINA_ARGS[@]}" | tee $smiles.txt

  rm "$temp_lig"
done < "$INPUT_FILE"

source venv_etoxpred/bin/activate
python docking_results.py