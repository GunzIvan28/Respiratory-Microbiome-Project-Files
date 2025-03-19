#! /usr/bin/env bash

### Selection for positive samples
# Check if correct number of arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_excel_file> <output_tab_file>"
    exit 1
fi

# Assign input and output file from user arguments
INPUT_FILE="$1"
OUTPUT_FILE="$2"

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file '$INPUT_FILE' not found!"
    exit 1
fi

# Convert Excel to CSV (assuming first sheet), then extract first column
TEMP_CSV=$(mktemp)
ssconvert --export-type=Gnumeric_stf:stf_csv "$INPUT_FILE" "$TEMP_CSV" 2>/dev/null

# Extract the first column (Sample IDs) and remove header
TEMP_SAMPLE_IDS=$(mktemp)
cut -d',' -f1 "$TEMP_CSV" | tail -n +2 > "$TEMP_SAMPLE_IDS"

# Ensure the output file exists and has a header
if [ ! -f "$OUTPUT_FILE" ]; then
    head -n 1 "$TEMP_CSV" > "$OUTPUT_FILE"
fi

# Get the list of already selected IDs
TEMP_SELECTED_IDS=$(mktemp)
cut -f1 "$OUTPUT_FILE" | tail -n +2 > "$TEMP_SELECTED_IDS"

# Keep selecting new participants until we reach 50
while [ "$(wc -l < "$TEMP_SELECTED_IDS")" -lt 50 ]; do
    # Select one new ID deterministically
    NEW_ID=$(shuf --random-source=<(yes 42) "$TEMP_SAMPLE_IDS" | grep -vxFf "$TEMP_SELECTED_IDS" | head -n 1)

    # If no new ID is found, break the loop
    if [ -z "$NEW_ID" ]; then
        echo "Warning: Not enough unique samples to reach 50."
        break
    fi

    # Append new ID to the selected list
    echo "$NEW_ID" >> "$TEMP_SELECTED_IDS"
    
    echo "Added participant: $NEW_ID" >> "$OUTPUT_FILE"
done

# Clean up temporary files
rm -f "$TEMP_CSV" "$TEMP_SAMPLE_IDS" "$TEMP_SELECTED_IDS"
sed -i 's/Added participant: //g' "$OUTPUT_FILE" && \
sed -i 's/"Record ID","Study ID","Has the patient received TB treatment for the last episode of TB?","Does the patient have any other chronic conditions?","Xpert Ultra result"//g' "$OUTPUT_FILE" && \
sed -i '/^[[:space:]]*$/d' "$OUTPUT_FILE"

echo "Selection complete. Selected samples saved to '$OUTPUT_FILE'."