# Lep fusion split finder 

A script to assign ancestral linkage units and/or identify fusion/fission events in Lepidopteran chromosomes based using a reference karyotype and BUSCO genes as markers.

### Running the script
`fusion_split_finder.py` takes the `full_table.tsv` output file for two species, along with an optional prefix (specified with -f, default "fsf"), e.g.:

```
python3 fusion_split_finder.py -q test_data/Aglais_io_full_table.tsv -r test_data/Melitaea_cinxia_full_table.tsv -f Agalis_io
```

### The output 

This will write three files:

`Aglais_chromosome_assignments.tsv`: a summary of the assignments for each scaffold in the query genome

`Aglais_fused_chromosomes.tsv`: a list of all fused chromosomes and their putative origins (note: any empty file is created when there are no fusion events inferred)

`Aglais_split_chromosomes.tsv`: a list of all split chromosomes and their putative origins (note: any empty file is created when there are no fission events inferred)


## Full usage 

```
usage: fusion_split_finder.py [-h] -r REFERENCE_TABLE -q QUERY_TABLE [-p MIN_PROPORTION] [-m REPORT_PROPORTION] [-f PREFIX]

optional arguments:
  -h, --help            show this help message and exit
  -r REFERENCE_TABLE, --reference_table REFERENCE_TABLE
                        full_table.tsv file for reference species
  -q QUERY_TABLE, --query_table QUERY_TABLE
                        full_table.tsv for query species
  -p MIN_PROPORTION, --min_proportion MIN_PROPORTION
                        Minimum proportion of BUSCO genes used to identify ancestral chromosomes
  -m REPORT_PROPORTION, --report_proportion REPORT_PROPORTION
                        Minimum proportion of BUSCOs required report for each fused/split chromosome
  -f PREFIX, --prefix PREFIX
                        Prefix for all output files
```
