# gnavigator
usage: gnavigator.py [-h] [-p PREFIX] [-d DB_DIR] [-n DB_NAME] [-t THREADS]
                     [-m GENETIC_MAP]
                     cDNA genome

Assess assembly quality and completeness using cDNA sequences

positional arguments:
  cDNA                  FASTA file of cDNA sequences to align to assembly
  genome                FASTA file of genome assembly to assess

optional arguments:
  -h, --help            show this help message and exit
  -p PREFIX, --prefix PREFIX
                        Prefix to use for intermediate and output files
                        [gnavigator]
  -d DB_DIR, --db_dir DB_DIR
                        Path to directory containing prebuilt GMAP index
                        [optional]
  -n DB_NAME, --db_name DB_NAME
                        Name of prebuilt GMAP index [optional]
  -t THREADS, --threads THREADS
                        Number of threads for GMAP alignment [1]
  -m GENETIC_MAP, --genetic_map GENETIC_MAP
                        Genetic map file as tsv with LG:cDNA pairs [optional]
