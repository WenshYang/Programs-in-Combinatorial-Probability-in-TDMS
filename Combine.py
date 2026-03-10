import os
import pandas as pd
from tqdm import tqdm

def combine_excel_sheets_in_directory(directory_path):
    all_data = {}  # Dictionary to store all data. Key = sheet_name, Value = combined data for that sheet

    # Extract a list of all .xlsx files in the directory
    excel_files = [f for f in os.listdir(directory_path) if f.endswith(".xlsx")]

    # Iterate through each file in the directory with a progress bar
    for filename in tqdm(excel_files, desc="Processing files"):
        filepath = os.path.join(directory_path, filename)
        with pd.ExcelFile(filepath) as xls:
            for sheet_name in xls.sheet_names:
                df = pd.read_excel(xls, sheet_name)

                # Keep ion label, fragment mass, and counts
                if 'b-ion Label' in df.columns:
                    df = df[['b-ion Label', 'b-ion Mass', 'b-ion Count']]
                elif 'y-ion Label' in df.columns:
                    df = df[['y-ion Label', 'y-ion Mass', 'y-ion Count']]
                else:  # internal fragment
                    df = df[['Internal Fragment Mass', 'Internal Fragment Count']]
                
                # Rename columns to make count column unique
                count_column_name = filename + "_" + df.columns[-1]
                df.columns = [*df.columns[:-1], count_column_name]

                # If the sheet name already exists in our dictionary, merge
                if sheet_name in all_data:
                    merge_keys = list(df.columns[:-1])
                    all_data[sheet_name] = all_data[sheet_name].merge(df, on=merge_keys, how='outer')
                else:
                    all_data[sheet_name] = df

    output_filename = os.path.join(directory_path, "combined_simulation.xlsx")

    # Write the combined data to the output file
    with pd.ExcelWriter(output_filename, engine='openpyxl') as writer:
        for sheet_name, data in all_data.items():
            # Sort by the columns (excluding the counts)
            sort_columns = list(data.columns[:-len(excel_files)])
            data = data.sort_values(by=sort_columns)
            data = data.fillna(0)  # Fill NaN values with 0
            data.to_excel(writer, sheet_name=sheet_name, index=False)

if __name__ == "__main__":
    directory_path = os.path.dirname(os.path.realpath(__file__))
    combine_excel_sheets_in_directory(directory_path)
    print(f"Combined data written to combined_simulation.xlsx in the directory {directory_path}")
