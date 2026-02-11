# MycoCosm Genome Downloader

Download genome FASTA files and GFF3 annotation files from all JGI MycoCosm portals.

This project is a fork of https://github.com/WesterdijkInstitute/MycoCosm_genome_downloader by [Jorge Navarro](https://github.com/jorgecnavarrom). It has been updated from the original (v1.3.0) to function with `ete4`.


## Requirements and Installation

* A JGI account
* Python 3.14 or newer
* `ete4` and `libxml2`

Clone this git repository, add it to the PATH, change into it, and install the dependencies. 

```
python -m venv venv
source venv/bin/activate
pip install ete4 lxml
```

**Optional but recommended:** save your JGI login credentials as a text file (`credentials.txt`) in the git repository to automatically log in for commands that require login credentials. The username should be on the first line and the password on the second with no other text. 


## Usage 

### Update the taxonomy database

As part of the output of this program, a "taxonomy" file is created, which includes the lineage for each species (this lineage is also used to create the folder structure for the output). For this, `ete4` is used. Use the following command to update the database.

```
python mycocosm_genome_downloader.py --update
```

### Obtain data files

It is optional but recommended to specify an output folder when downloading files. Otherwise, they will be downloaded into the directory `output` in the git repository. 

1. Obtain the [list of genomes](https://mycocosm.jgi.doe.gov/ext-api/mycocosm/catalog/download-group?flt=&seq=all&pub=all&grp=fungi&srt=released&ord=asc) using either the default output directory or specifying your own.

```
python mycocosm_genome_downloader.py --getgenomelist
python mycocosm_genome_downloader.py --getgenomelist -o $out 
```

The genomes list is a comma-separated file which includes all current projects ("**portals**") in MycoCosm (`MycoCosm_Genome_list.csv`). **Note:** it uses ISO-8859-15 encoding.

2. Obtain a list of files to download using either the default output directory or specifying your own and optionally including a credentials files for easier log in.

```
python mycocosm_genome_downloader.py --getxml
python mycocosm_genome_downloader.py --getxml -o $out
python mycocosm_genome_downloader.py --getxml -o $out -j $credentials 
```

The XML file is named `MycoCosm_data.xml`. It takes about an hour to download and is around 70 MB in size. After downloading, the XML file is re-formatted to be human-readable. If there is already a file called `MycoCosm_data.xml` in the output folder, it will be updated with any missing data.


### Download data

Once you have the genomes and avialable file lists, use them to download the data. It took about four hours to download 2864 portals. The compressed files sum up to about 40 GB.

```
python mycocosm_genome_downloader.py --csv [path to MycoCosm_Genome_list.csv] --xml [path to MycoCosm_data.xml] -o $out 
```

This will create a directory structure that follows the fungal taxonomy tree at MycoCosm [main page](https://mycocosm.jgi.doe.gov/mycocosm/home).

**Optional:** use parameter `--simulate` to create the directory structure without downloading any files.

### Re-use data

If you already have downloaded data, you can simply use your local copy of the files instead of downloading them. This makes it easier to update with new genomes. Use the option `--getprevious`, which will scan a base folder with previous results and create a file called `previously_downloaded_files.tsv`. Use this file with option `--previous` to skip files already downloaded.


## Output

Each 'leaf' folder within the folder structure represents a genome from MycoCosm and contains the genome assembly FASTA and gene annotation GFF3 files (compressed as `.gz`). Additionally, three files are found inside the base output folder:
* `JGI_taxonomy.tsv`. Columns: `Short name` (Portal), `Accession` (same as short name in this case), `TaxId` (from NCBI), `Name` (project name), `Path` (relative path from output folder to Portal), `Assembly file`, `GFF file`, `lineage` (comma-separated)
* `JGI_download_list_[date].txt`: All the files that were downloaded
* `List_gene_gff_filenames.txt`: All gff files for each portal

## Known limitations and notes

* The incorrect GFF might have been chosen
* Some old GFF files can't be processed correctly by other tools (e.g. `bcbio-gff`)
* Doesn't handle dikaryon genomes ("primary/secondary alleles")

## More details

### Target download folders

These are the folders that are scanned for files, as they appear in each Portal's download section:

* Assemblies: "Genome Assembly (unmasked)" (previously called "Assembled scaffolds (unmasked)")
* Gene annotations: "Filtered Models ('best')"

**Note:** Some genomes seem to be only available in the "masked" version (e.g. `Pyrtr1`, `Altbr1`, etc.)

### Choosing the right annotation file

There are usually a few files to choose from in the "Filtered Models ('best')" folder. The criteria to choose a single file is implemented in the following order:
* In some cases, the choosing was done manually and incorporated to the `hardcoded_gff_files.tsv` file. If the Portal has an entry here, use the indicated file. This is a work in progress!
* Skip some specific files: `Aciri1_meta_GeneCatalog_genes_20111216.gff.gz`, `Exoaq1_GeneCatalog_20160901.gff3.gz`, `Exoaq1_GeneCatalog_20160828.gff3.gz`, `Fonpe1_GeneCatalog_20160901.gff3.gz`, `Copmic2_FM1_removed_alleles.gff.gz`
* Skip files with `proteins` or `secondary_alleles` in their name
* Skip files with the following extensions: `gtf.gz`, `tgz`, and any other that is *not* `gz`
* Prefer `gff3` files
* If we have two `gff3` (or two `gff`) files, prefer the most recent

### Ignored assembly files / portals

* Files with `MitoAssembly`, `MitoScaffolds`, `PrimaryAssemblyScaffolds` or `SecondaryAssemblyScaffolds` in their name
* Excluded assemblies (hardcoded). Ignore these filenames as they're not related to assemblies, are old versions or are assemblies of meta-samples: `1034997.Tuber_borchii_Tbo3840.standard.main.scaffolds.fasta.gz`, `Spofi1.draft.mito.scaffolds.fasta.gz`, `Patat1.draft.mito.scaffolds.fasta.gz`, `PleosPC9_1_Assembly_scaffolds.fasta.gz`, `Neuhi1_PlasmidAssemblyScaffolds.fasta.gz`, `CocheC5_1_assembly_scaffolds.fasta.gz`, `Alternaria_brassicicola_masked_assembly.fasta.gz`, `Aciri1_meta_AssemblyScaffolds.fasta.gz`, `Rhoto_IFO0880_2_AssemblyScaffolds.fasta.gz`
* Ignored portals (hardcoded). Metaprojects or old versions: `Rhoto_IFO0880_2`, `Aciri1_meta`, `Pospl1`
