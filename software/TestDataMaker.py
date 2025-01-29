import pandas as pd

def process_dataframe(folder_path,file_path):
    """
    Processes a DataFrame from an Excel file:
    1. Creates versions with each column removed (except ignored columns).
    2. Creates a version with column positions swapped.
    3. Saves all versions as new Excel files.

    Parameters:
        file_path (str): Path to the input Excel file.
    """
    # Load the dataframe from the Excel file
    df = pd.read_csv(folder_path + file_path)
    print('read')

    # Columns to ignore
    ignore_columns = ['Batch', 'Sample', 'Percent_PostTrim', '% Uniquely Aligned Reads', 'Percent_Exonic', 'Perc_Aligned_Reads_Overlapping_rRNA']

    # Create versions with each column removed (except the ignored ones)
    for column in df.columns:
        print(column)
        if column not in ignore_columns:
            df_dropped = df.drop(columns=[column])
            df_dropped.to_csv(f'{folder_path}SCRIPTB11dropped_{column}.csv', index=False)

    # Create a version with column positions swapped
    # Example: Swap the first and second columns that are not in the ignore list
    columns_to_swap = [col for col in df.columns if col not in ignore_columns]
    if len(columns_to_swap) >= 2:
        df_swapped = df.copy()
        df_swapped[[columns_to_swap[0], columns_to_swap[1]]] = df_swapped[[columns_to_swap[1], columns_to_swap[0]]]
        df_swapped.to_csv(f'{folder_path}SCRIPTB11swapped_columns.csv', index=False)

    print("Files have been created successfully.")

# Example usage
folderpath = '../data/SCRIPT/'
file_path = 'SCRIPT_B11_QCTable.csv'  # Replace with your actual file path
print('starting')
process_dataframe(folderpath,file_path)
