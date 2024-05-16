# All the CF Metagenomes in the SRA

We get all the metagenomes that have CF data from the SRA archive, and use `mash` to compare them to our analysis.

First, get a list of all the [BioProjects at NCBI](https://www.ncbi.nlm.nih.gov/bioproject/?term=cystic+fibrosis) related to `cystic fibrosis`. Save the accessions as a file `cf_bioprojects.txt` for those (at the time of writing, there are 964 bioprojects. You could probably manually curate them but I didn't.

Upload those to [Google BigQuery](https://console.google.com) and use the query shown to get a list of the SRA accessions for those bioprojects that are not amplicon sequences:

```sql
create temp table AMPLICON(acc STRING) as select acc as amplicon from `nih-sra-datastore.sra.metadata` where assay_type = 'AMPLICON' or libraryselection = 'PCR';
create temp table BIOPROJ(bioproject STRING) as SELECT string_field_0 FROM `sra-searches.cystic_fibrosis.bioproject_accs` WHERE string_field_0 IS NOT NULL;
select * from `nih-sra-datastore.sra.metadata` where acc not in (select acc from AMPLICON) and bioproject in (select bioproject from BIOPROJ) and (librarysource = "METAGENOMIC" or librarysource = 'METATRANSCRIPTOMIC' or organism like "%microbiom%" OR organism like "%metagenom%");
```

This provides 4,008 results which I save to Google Drive as JSONL (newline delimited), and download the file here:

```bash
conda activate rclone
rclone copy GoogleDrive:bq-results-20230331-064211-1680245088432.json .
mv bq-results-20230331-064211-1680245088432.json bigquery.json
gzip bigquery.json
conda deactivate
```

Now we have a results file called `bigquery.json.gz` we can parse.

For bash, I prefer to use [jq](https://cameronnokes.com/blog/jq-cheatsheet/) to parse `json`, and so we can extract all the accession numbers:

```bash
gunzip -c bigquery.json.gz  | jq '.acc' | sed -e 's/"//g' > sra_ids.txt
```

(Of course, if that is all we want, we could just select the accessions in the SQL and download those as text. This way we have the SRA metadata too!)

Now we can use [fasterq-dump](https://edwards.flinders.edu.au/fastq-dump/) to download the data, using either of these two approaches:

- `fastqdump.sh`: read a file called `cf_metagenomes.txt` and download all the data in `fasta` format from SRA using [fasterq-dump](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump). We use fasta format because we don't do anything with the quality scores.
- `fastqdump.slurm`: download the data in `cf_metagenomes.txt` in parallel. This works, but is not kind to your internet connection!

After downloading the data, we had to use `split_fasta.slurm` because when we download the data from SRA using `--fasta-unsorted` we get the R1 and R2 reads in one file. This splits the fasta files into `R1` and `R2` files, in parallel. It is probably better to use `--fasta` and get the data separated.

I used `bigquery_json.py` to parse the lat lon and geolocation data from the SRA. Then I only downloaded the samples with a lat lon since that's all we need (for now)

> Note that these scripts are based on [atavide-lite](https://github.com/linsalrob/atavide_lite) and we used that for a lot of the data analysis. However, we also ran some of this on different clusters, including [NCI](https://www.nci.org.au/) which uses PBS and (Pawsey)[https://pawsey.org.au/] which uses slurm. For the PBS scripts, we generally optimise by writing to the temporary `$PBS_JOBFS` filesystem and then copying the data to the final location. 

Other scripts included here:

- `fastp.pbs`: Runs [fastp](https://github.com/OpenGene/fastp) using the PBS queue on NCI. 
- `fastp.sh`: an array job that reads a file called `R1_reads.txt` and processes the reads in that file
- `mash.sh`: calculate sketches for all the read data


