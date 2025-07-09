# VMR-BIOINFORMATICS-PIPELINE-NF

A nextflow pipeline for bioinformatics analysis from the GRDI-AMR2 data.




## Overview

The nextflow script `download_databases.nf` is used to download the databases
used by RGI and MOB-suite (if needed).

The nextflow script `pipeline.nf` runs starAMR, RGI, Mobsuite and Abricate. Both RGI and MOB-suite output their results in the form of tabular
files. These results are merged, associating each detected AMR determinant with
MOB-suite's results. This enables the user to identify if the AMR determinants
are present on the chromosome, or on a plasmid. If multiple samples are passed
to the program, the results will be collated into a single table.  StarAMR will identify MLST and Resfinder genes.
Sistr and Ectyper can be used depending on the species isolated in the project to perform the serotyping.
Abricate is configured to search into the VFDB database to retrieve virulence genes .Once detected a Klebsiella species, Keborate will run 
automatically with the parameters for either Pneumoniae or Oxytocal complexes.
On top of that if Escherichia coli is identified, ectyper and virulence finder will automatically be ran.


## Installation

Dependencies:

* NextFlow

* If using docker/singularity profile: 
    - Docker
    - Singularity

## Usage

Clone this pipeline into a working directory.
### Slurm compatibility

The pipeline is compatible with slurm, that is needed when used at PHAC/Waffles server. Just add the "-c path/to/vmr-bioinformatics-pipeline-nf/nextflow_slurm.config"" parameter.


### Download databases

Both MOB-suite and RGI require databases to be downloaded in order to run the
tools. MOB-RGI-NF provides a pipeline to download the databases if the user does
not already have them. You may specify `--download_mobDB` to download the
MOB-suite databases, `--download_card_json` to download the CARD json file, or
`--download_all` to download both.

```bash
nextflow run download_databases.nf \
    -profile [docker|singularity] \
    --download_all
```

By default, this command will download both databases into a `databases` folder
in nextflow's `work` directory. This keeps the databases out of the way and
prevents them from being deleted when nextflow runs are cleaned during analyses.

However, the user can specify their own locations instead. In this case, make
sure to also specify the location of these databases when running the analysis
pipeline.

```bash
nextflow run download_databases.nf \
    -profile [conda|docker|singularity] \
    --download_all \
    --download_mobDB "/path/to/mobDB/directory" \
    --download_card_json "/path/to/card/directory" \
    --download_vfDB "path/to/virulencefinder/directory"
```

By default, VMR-BIOINFORMATICS-PIPELINE-NF will not overwrite the databases if
it detects that they are already installed at the given path; set `--overwrite`
to overwrite the existing databases if needed.

### Detect AMR determinants and plasmids

#### Basic Usage

Use the script `pipeline.nf` to run mobSuite and RGI. The input should be
assembled contigs. Use a bash glob to designate multiple sequences to be run by
the pipeline. Note that the:wq file names will be used to designate sample ID.

Make sure to select a run profile, depending on the container
you want to use. Currently conda (using mamba) and docker are supported. The
pipeline will handle downloading and initializing the containers. Therefore, the
very first run will take some time as the proper containers are downloaded.

```bash
nextflow run pipeline.nf \
    -profile [conda|docker|singularity] \
    --species "species" \  default "Escherichia coli"
    --contigs "dir/to/sequences/*.fasta"
    --sistr [true|false] \ default: false
    --mobDB PATH to the DIRECTORY of the MOB-suite databases
    --card_json PATH to CARD's card.json file.
    --vfinder [listeria|s.aureus_exoenzyme|s.aureus_hostimm|s.aureus_toxin| \
               stx|virulence_ecoli|virulence_ent|virulence_entfm]
    
```

This command assumes that the databases have first been downloaded into the
default directory by the `download_databases.nf` pipeline. Otherwise, specify
the location of the databases. Use the option `--mobDB` to specify the path to
the parent directory of the MOB-suite databases, and `--cardDB` to specify the
path to the card.json file



The results will be collected in the directory `results` in the project
directory by default. This can be changed with the option `--outDir`.

### Tools Included in the Pipeline

- **StarAMR**: Antimicrobial resistance detection using ResFinder, PlasmidFinder, MLST.
- **MOB-suite**: Plasmid typing and reconstruction.
- **RGI**: AMR gene detection using CARD database.
- **ECTyper**: E. coli serotyping.
- **SISTR**: Salmonella serotyping.
- **Kleborate**: Klebsiella virulence and resistance profiling.
- **Prokka**: Genome annotation.
- **Refseq_masher**: Species identification through similarities.
- **VirulenceFinder**: Detection of virulence genes.
- **Abricate**: Screening contigs against the VFDB database.
- **ICEberg**: Identification of integrative and conjugative elements (ICEs).
- **IntegronFinder**: Identification of integrons in genomes.
- **IslandPath-DIMOB**: Prediction of genomic islands using sequence composition.
- **DIGIS**: Focused detection of **insertion sequences (IS elements)** using a curated IS database.



=======
# vmr-bioinformatics-pipeline
>>>>>>> 8ba31a4a7f49c970d128063b5c48fff70a063d68
