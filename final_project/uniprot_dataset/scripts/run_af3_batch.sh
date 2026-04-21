#!/bin/bash

AF3_DIR="/media/Data/shanezhu/alphafold3_install/alphafold3"
INPUT_DIR="$AF3_DIR/input/uniprot"
OUTPUT_BASE="/media/Data/shanezhu/proteina_hands_on/final_project/uniprot_dataset/structures/from_af3"

mkdir -p $OUTPUT_BASE

cd $AF3_DIR

for json_file in $INPUT_DIR/*.json; do
    uid=$(basename $json_file .json)
    
    # skip already done 
    if [ -d "$OUTPUT_BASE/$uid" ]; then
        echo "Skipping $uid (already done)"
        continue
    fi
    
    echo "Folding $uid..."
    python run_alphafold.py \
        --json_path=$json_file \
        --model_dir=$AF3_DIR/cache \
        --output_dir=$OUTPUT_BASE/$uid
done

echo "All done!"
