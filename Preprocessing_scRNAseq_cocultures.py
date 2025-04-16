#!/usr/bin/python

import os
import pandas as pd
import subprocess
import glob
import shutil

# Paths
input_tsv = "/LNStroma/202503_LNS019_scRNAseq/data/42585_meta.tsv"
output_base = "/LNStroma/202503_LNS019_scRNAseq/202503_CellRanger_logs"
fastq_base = "/LNStroma/202503_LNS019_scRNAseq/data/250305_A01382_0631_AH2MHTDSXF"
feature_ref = "/LNStroma/202503_LNS019_scRNAseq/hashing_3prime.csv"
transcriptome = "/refdata/refdata-gex-GRCh38-2024-A"

#Read metadata file
df = pd.read_csv(input_tsv, sep="\t")

#Define samples
samples = ['S1', 'S2', 'S3']

for sample in samples:
    sample_df = df[(df["SAMPLE_NAME"] == sample) | (df["SAMPLE_NAME"] == f"Hash {sample}")]
    
    for lane in sorted(sample_df["LANE_NO"].unique()):
        lane_df = sample_df[sample_df["LANE_NO"] == lane]
        fastq_folders = {fq.rsplit("_R", 1)[0] for fq in lane_df["FASTQ_FILE"]}
        
        ##### fastq renaming (within the same folder) based on cellranger requirments #####
        # ie SampleName_S1_L001_R1_001.fastq.gz
        # I had to do this one, but not again so I dont run it now 
        for fastq_folder in fastq_folders:
            fastq_folder_full = f"{fastq_base}/{fastq_folder}/fastq"
            fastq_files = glob.glob(os.path.join(fastq_folder_full, "*.fastq.gz"))

            for fastq_file in fastq_files:
                filename = os.path.basename(fastq_file)

                #extract read pair (R1 or R2)
                read_pair = filename[filename.index("_R") : filename.index("_R") + 3]  # e.g., _R1 or _R2
        
                #ensure naming: SampleName_S1_L001_R1_001.fastq.gz
                new_name = f"{sample}_L{lane:03}{read_pair}_001.fastq.gz"
        
                #define renamed path
                renamed_file_path = os.path.join(fastq_folder_full, new_name)
        
                #new_name = f"{sample}_{lane}_L001" + filename[filename.index("_R"): ]
                #renamed_file_path = os.path.join(fastq_base, fastq_folder, 'fastq', new_name)

                #rename the file within the same directory
                #shutil.move(fastq_file, renamed_file_path)
                print(f"Copied: {fastq_file} to {renamed_file_path}")

        # GEX
        #lane_gex = f"{sample}_{lane}_L001"
        lane_gex = list(lane_df[lane_df["SAMPLE_NAME"] == sample]["FASTQ_FILE"])[0].rsplit("_R", 1)[0]
        fastq_path_gex = f"{fastq_base}/{lane_gex}/fastq"

        # hashtags (Feature Barcode)
        #lane_hash = f"{sample}_{lane}_L001"
        lane_hash = list(lane_df[lane_df["SAMPLE_NAME"] == f"Hash {sample}"]["FASTQ_FILE"])[0].rsplit("_R", 1)[0]
        fastq_path_hash = f"{fastq_base}/{lane_hash}/fastq"

        #create output dir
        output_dir = os.path.join(output_base, f"{sample}_{lane}")
        #os.makedirs(output_dir, exist_ok=True)
        
        # Generate multi-config.csv
        multi_config_path = os.path.join(output_dir, "multi-config.csv")

        
        multi_config_content = f"""\
        [gene-expression]
        reference,{transcriptome}
        create-bam,false
        
        [libraries]
        fastq_id,fastqs,feature_types
        {sample}_{lane},{fastq_path_gex},Gene Expression
        {sample}_{lane},{fastq_path_hash},Multiplexing Capture
        
        [feature]
        reference,{feature_ref}

        [samples]
        sample_id,cmo_ids
        C1,HASHTAG1
        C2,HASHTAG2
        C3,HASHTAG3
        C4,HASHTAG4
        C5,HASHTAG5
        C6,HASHTAG6
        C7,HASHTAG7
        """
        
        with open(multi_config_path, "w") as f:
            f.write(multi_config_content)

        bsub_cmd = f"""bsub -q "verylong" -n 32 -R "rusage[mem=200G]" \
        -o {output_dir}/cellranger_%J.log \
        "module load cellranger/8.0.1 && \
        cellranger multi --id={sample}_{lane} \
        --csv={multi_config_path}"
        """
        
        os.system(bsub_cmd)
        print(f"{sample}_{lane} submitted")
