# Functional summary for the MinION data

All data comes from `*.unmapped.fasta.functions` files in [CFTeam:CF_Hackathon_2023/Minion_read_based_annotations](https://flinders.sharepoint.com/:f:/r/sites/CysticFibrosisResearch2/Shared%20Documents/CF_Hackathon_2023/Minion_read_based_annotations?csf=1&web=1&e=xep7WU) but they have been passed through the [overlap filter](https://github.com/linsalrob/CFHackathon2023/blob/main/parse_m8/non_overlapping_mmseqs_hits.py) using a flags `-c -e 1e-20`: 

 ``` 
 python CFHackathon2023/parse_m8/non_overlapping_mmseqs_hits.py -c -f $FILE -e 1e-20 
 ``` 

The processed files are `*.unmapped.fasta.functions.eVal20.nooverlap.gz` in [CFTeam:CF_Hackathon_2023/Minion_read_based_annotations](https://flinders.sharepoint.com/:f:/r/sites/CysticFibrosisResearch2/Shared%20Documents/CF_Hackathon_2023/Minion_read_based_annotations?csf=1&web=1&e=xep7WU)
