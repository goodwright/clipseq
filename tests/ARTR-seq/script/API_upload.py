import flowbio
import os

client = flowbio.Client()
client.login("hanzhong", "H5HUeHc!5rDaRzs")


# sample dir path
sample_dir = "/nemo/lab/ulej/data/STPs/babs/inputs/hanzhong.bai/asf/PM24162/20240807_LH00442_0034_A22MYGKLT3/fastq"

# Loop through samples A9 to A24
for i in range(12, 25):
    sample_name = f"BAI7152A{i}"
    reads1 = os.path.join(sample_dir, f"{sample_name}_S{i+338}_L008_R1_001.fastq.gz")
    reads2 = os.path.join(sample_dir, f"{sample_name}_S{i+338}_L008_R2_001.fastq.gz")  # Assuming paired-end reads

    # Upload sample
    sample = client.upload_sample(
        sample_name,
        reads1,
        reads2,  # Optional
        progress=True,
        metadata={
            "Category": "CLIP",
            "Organism": "Human",
            "Project": "ARTR-seq-20240924",
            "5' Barcode Sequence": "NNNNNGGG" 
        }
    )
    print(sample)
