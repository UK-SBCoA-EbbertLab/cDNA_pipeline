#!/usr/bin/env python3
import sys
import os
import pandas as pd

def main():
    # 1. Parse the command-line arguments into a list (excluding the script name)
    args = sys.argv[1:]
    if len(args) < 1:
        sys.exit("Usage: {} file1 [file2 ...]".format(sys.argv[0]))
    
    # 2. Open the first file as a dataframe using tab separation.
    try:
        merged_df = pd.read_csv(args[0], sep="\t")
    except Exception as e:
        sys.exit("Error reading {}: {}".format(args[0], e))
    
    # 3. Store the number of rows and columns into variables.
    initial_rows, initial_cols = merged_df.shape

    # 4. Check for a "TXNAME" column and set the merge columns accordingly.
    if "TXNAME" in merged_df.columns:
        merge_columns = ["TXNAME", "GENEID"]
    else:
        merge_columns = "GENEID"
    
    # 5. Loop through the remaining command-line arguments (starting at index 1).
    for filename in args[1:]:
        try:
            # 6. Open the current file as a dataframe using tab separation.
            current_df = pd.read_csv(filename, sep="\t")
        except Exception as e:
            sys.exit("Error reading {}: {}".format(filename, e))
        
        # Save the current number of columns (before merge) for the sanity check.
        prev_num_cols = merged_df.shape[1]
        
        # 7. Inner merge the current file with the merged dataframe on the merge columns.
        merged_df = pd.merge(merged_df, current_df, on=merge_columns, how='inner')
        
        # 8. Sanity checks:
        #    - The number of rows should remain constant.
        if merged_df.shape[0] != initial_rows:
            sys.exit("Error: After merging with '{}', row count changed from {} to {}."
                     .format(filename, initial_rows, merged_df.shape[0]))
        #    - The number of columns should have increased by exactly 1.
        if merged_df.shape[1] != prev_num_cols + 1:
            sys.exit("Error: After merging with '{}', column count did not increase by 1 (from {} to {})."
                     .format(filename, prev_num_cols, merged_df.shape[1]))
    
    # 9. Outside the loop, generate the output file name.
    #     Take the basename of the first file (in case it includes a path).
    base_name = os.path.basename(args[0])
    #     Split the name by '_' and join the last two parts.
    parts = base_name.split('_')
    if len(parts) >= 2:
        output_filename = '_'.join(parts[-2:])
    else:
        # Fallback: if there aren’t enough underscores, use the whole basename.
        output_filename = base_name

    # Save the merged dataframe as a TSV file without the index.
    merged_df.to_csv(output_filename, sep="\t", index=False)
    print("Merged dataframe saved to:", output_filename)

if __name__ == "__main__":
    main()

