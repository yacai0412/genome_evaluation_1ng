# genome_evaluation_1ng
code for bioinfomatics evaluation of Hi-Fi data


## Pipeline for the Third-generation Sequencing data
Using generate_shell_script.py to generate a.sh file:
```
python3 /public5/home/sch5655/wmx/generate_shell_script.py filenames.tsv PTA_20250811

```

Please read a.sh first.
Put all the .sh .py .pl files into the data_process directory, running a.sh.




## Pipeline for the Next-generation sequencing data (pair-end reads):

Using automatical process pipeline shell code:
```
EQ_42_1_fq_gz=/nfs119/rd1/caiya/wmx/MDA_20251223/20251218_E251218002_U4924_FQ251200955/RawData/EQ-42/E251218002_L01_EQ-42_1.fq.gz
EQ_42_2_fq_gz=/nfs119/rd1/caiya/wmx/MDA_20251223/20251218_E251218002_U4924_FQ251200955/RawData/EQ-42/E251218002_L01_EQ-42_2.fq.gz

NEB_phi_1_fq_gz=/nfs119/rd1/caiya/wmx/MDA_20251223/20251218_E251218002_U4924_FQ251200955/RawData/NEB-phi/E251218002_L01_NEB-phi_1.fq.gz
NEB_phi_2_fq_gz=/nfs119/rd1/caiya/wmx/MDA_20251223/20251218_E251218002_U4924_FQ251200955/RawData/NEB-phi/E251218002_L01_NEB-phi_2.fq.gz



bash /rd/caiya/wmx/run_shortread_pipeline.sh \
    --workdir /rd/caiya/wmx/testdata/dataprocess \
    --sample \
        --fq1-short EQ-42_1.fq.gz --fq1-abs /nfs119/rd1/caiya/wmx/MDA_20251223/20251218_E251218002_U4924_FQ251200955/RawData/EQ-42/E251218002_L01_EQ-42_1.fq.gz \
        --fq2-short EQ-42_2.fq.gz --fq2-abs /nfs119/rd1/caiya/wmx/MDA_20251223/20251218_E251218002_U4924_FQ251200955/RawData/EQ-42/E251218002_L01_EQ-42_2.fq.gz \
        --label EQ-42 \
    --sample \
        --fq1-short NEB-phi_1.fq.gz --fq1-abs /nfs119/rd1/caiya/wmx/MDA_20251223/20251218_E251218002_U4924_FQ251200955/RawData/NEB-phi/E251218002_L01_NEB-phi_1.fq.gz \
        --fq2-short NEB-phi_2.fq.gz --fq2-abs /nfs119/rd1/caiya/wmx/MDA_20251223/20251218_E251218002_U4924_FQ251200955/RawData/NEB-phi/E251218002_L01_NEB-phi_2.fq.gz \
        --label NEB-phi 

```

