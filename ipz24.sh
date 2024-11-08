#!/bin/bash

FLAG_FILE="initialization_flagfile"
default_file="tcm600_nr.smi"

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
    mv "$input_file" "eToxPred/$input_file"
else
    input_file=$default_file
fi

# Rest of the script goes here
echo "Running toxicity prediction..."

source venv_etoxpred/bin/activate
echo "Using file: $input_file"
cd eToxPred
python etoxpred_predict.py --datafile "$input_file" --modelfile etoxpred_best_model.joblib --outputfile ../results.csv
cd ..