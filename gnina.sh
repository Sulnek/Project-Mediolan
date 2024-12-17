#!/usr/bin/env bash

INPUT_FILE=$(realpath "$1")
REC=$(realpath "$2")

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
    -r "$REC"
    -l "$temp_lig"
    # --flex "$REC_FLEX" 
    --autobox_ligand "$temp_lig" # pewnie co innego
    # --flexdist_ligand "$LIG" 
    # --flexdist 3.5
    -o whole_docked.sdf.gz
    --exhaustiveness 128
  )

  gnina "${GNINA_ARGS[@]}"

  rm "$temp_lig"
done < "$INPUT_FILE"




