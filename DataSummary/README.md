# Data summary

We have sequenced a lot of different things, and this folder has the summary counts of our sequencing.

# [MGI Counts Summary](MGI_Sequence_Counts.tsv.gz)

This tab-separated file is a summary of the MGI sequence reads. It has the following columns:

Column name | Meaning | Sum
--- | --- | ---
Sample | The sample identifier. Note that we have counted R1 and R2 files separately! | 127 R1 and R2 samples
Raw sequence | Number of reads from the sequencer | 1,178,307,208 reads
After QC/QA | Numbers of reads after `fastp` clean up | 1,176,881,906 reads
Mapped to the human genome | Number of reads that mapped to the human genome using `minimap2` | 672,903,626 reads
Did not map to the human genome | Number of reads that did _not_ map to the human genome | 515,964,090 reads
Had a match to UniRef50 | Number of reads that mapped to the UniRef50 protein database using MMSeqs2 | 429,542,979 reads


# [MinION Sample Counts](MinION_Sample_Counts.tsv)

This three-column table has sequence ID, # reads, #bp.

