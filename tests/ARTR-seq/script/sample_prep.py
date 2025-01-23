import os
import pandas as pd

# Define the directory path
sample_dir = "/nemo/lab/ulej/data/STPs/babs/inputs/hanzhong.bai/asf/PM24162/20240807_LH00442_0034_A22MYGKLT3/fastq"

# Create a DataFrame with the specified columns
df = pd.DataFrame(columns=["group", "replicate", "fastq_1", "fastq_2"])

# Loop through samples A9 to A24
for i in range(9, 25):
    sample_name = f"BAI7152A{i}"
    fastq_1 = os.path.join(sample_dir, f"{sample_name}_S{i+338}_L008_R1_001.fastq.gz")
    replicate = 2 - (i % 2)
    df = df.append({"group": "", "replicate": replicate, "fastq_1": fastq_1, "fastq_2": ""}, ignore_index=True)

# Save the DataFrame as a CSV file
output_csv = "sample_sheet.csv"
df.to_csv(output_csv, index=False)

print(f"Sample sheet created: {output_csv}")

