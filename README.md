# VMR-BIOINFORMATICS-PIPELINE-NF

A nextflow pipeline for bioinformatics analysis from the GRDI-AMR2 data.




## Overview

The nextflow script `download_databases.nf` is used to download the databases
used by RGI, MOB-suite, ICEberg, PHASTER, and BacMet2 (if needed).

The nextflow script `pipeline.nf` runs StarAMR, RGI, MOB-suite, Abricate, ICEberg,
PHASTER, IntegronFinder, IslandPath, DIGIS, and Prokka. Both RGI and MOB-suite output
their results in the form of tabular files. These results are merged, associating each
detected AMR determinant with MOB-suite's results. This enables the user to identify if
the AMR determinants are present on the chromosome or on a plasmid. If multiple samples
are passed to the program, the results will be collated into a single table.

StarAMR identifies MLST and ResFinder genes. Sistr and ECTyper can be used depending on
the species to perform serotyping. Once a Klebsiella species is detected, Kleborate runs
automatically with parameters for either the Pneumoniae or Oxytoca complexes. If
Escherichia coli is identified, ECTyper and VirulenceFinder run automatically.

Abricate is configured to search the VFDB database for virulence genes and the BacMet2
(confirmed and predicted) databases for biocide/metal resistance genes.

PHASTER (PHAge Search Tool Enhanced Release) protein database is searched via DIAMOND
blastp to identify phage-related proteins in annotated genes.


## Installation

Dependencies:

* NextFlow

* If using docker/singularity profile:
    - Docker
    - Singularity

## Usage

Clone this pipeline into a working directory.

### Slurm compatibility

The pipeline is compatible with SLURM, which is needed when used at the PHAC/Waffles
server. Add `-c path/to/vmr-pipeline/nextflow_slurm.config` to your nextflow command.


### Download databases

Several tools require databases to be downloaded before running. `download_databases.nf`
provides a pipeline to download all of them. Use `--download_all` to download everything,
or individual flags to download specific databases.

```bash
nextflow run download_databases.nf \
    -profile [docker|singularity] \
    --download_all
```

Individual download flags:

| Flag | Database |
|------|----------|
| `--download_mobDB` | MOB-suite databases |
| `--download_card_json` | CARD json file (for RGI) |
| `--download_vfDB` | VirulenceFinder database |
| `--download_prokkaDB` | Prokka HMM databases (TIGRFAM + PGAP) |
| `--download_iceberg` | ICEberg nucleotide and protein databases |
| `--download_phaster` | PHASTER prophage/phage protein database |
| `--download_bacmet2` | BacMet2 confirmed and predicted databases (for Abricate) |

By default, databases are downloaded into a `databases` folder in Nextflow's `work`
directory. You can specify custom locations:

```bash
nextflow run download_databases.nf \
    -profile [conda|docker|singularity] \
    --download_all \
    --mobDB            "/path/to/mobDB" \
    --card_json        "/path/to/card" \
    --virulencefinderDB "/path/to/virulencefinder_db" \
    --prokkaDB         "/path/to/prokka_db" \
    --iceberg          "/path/to/iceberg" \
    --phaster          "/path/to/phaster" \
    --bacmet2DB        "/path/to/bacmet2_abricate_db"
```

Set `--overwrite` to force re-download if databases are already present at the given path.


### Detect AMR determinants and plasmids

#### Basic Usage

Use `pipeline.nf` to run the full analysis. The input should be assembled contigs.
Use a bash glob to specify multiple sequences. File names are used as sample IDs.

Select a run profile depending on the container system available. The first run will
take time as containers are downloaded.

```bash
nextflow run pipeline.nf \
    -profile [conda|docker|singularity] \
    --contigs "dir/to/sequences/*.fasta"
```

#### Pipeline Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--contigs` | Path to input contigs (glob supported) | *(required)* |
| `--outDir` | Output directory | `results/` |
| `--num_threads` | Threads per process | `1` |
| `--species` | Species for StarAMR samplesheet | `"Escherichia coli"` |
| `--plasmids_only` | Run RGI only on plasmid contigs | `false` |
| `--ectyper` | Force ECTyper for all samples (auto-runs for *E. coli*) | `false` |
| `--sistr` | Enable SISTR serotyping | `false` |
| `--vfinder` | VirulenceFinder database(s) to use (see below) | `""` |
| `--mobDB` | Path to MOB-suite database directory | workDir default |
| `--card_json` | Path to CARD card.json file | workDir default |
| `--virulencefinderDB` | Path to VirulenceFinder database directory | workDir default |
| `--prokkaDB` | Path to Prokka HMM database directory | workDir default |
| `--icebergDB` | Path to ICEberg database directory | workDir default |
| `--phasterDB` | Path to PHASTER database directory | workDir default |
| `--bacmet2DB` | Path to BacMet2 abricate database directory | workDir default |
| `--blast_MGEs_minid` | Minimum % identity for MGE BLAST searches (ICEberg, PHASTER) | `35` |
| `--blast_MGEs_mincov` | Minimum % coverage for MGE BLAST searches (ICEberg, PHASTER) | `80` |

#### Using pre-computed MOB-suite results (Mikrokondo integration)

If MOB-suite was already run externally (e.g. via the Mikrokondo pipeline), you can
skip the internal `mob_recon` step and point the pipeline to the existing results:

```bash
nextflow run pipeline.nf \
    -profile [conda|docker|singularity] \
    --contigs "dir/to/sequences/*.fasta" \
    --mob_recon_dir "/path/to/mikrokondo/results" \
    --mob_suffix "contig_report.mob.recon.annotation.filled.txt"
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--mob_recon_dir` | Path to directory containing per-sample MOB-suite output folders | `""` (runs mob_recon internally) |
| `--mob_suffix` | Filename suffix of the MOB-suite contig report inside each sample folder | `"contig_report.mob.recon.annotation.filled.txt"` |

The pipeline expects the directory structure `<mob_recon_dir>/<sample_name>/*<mob_suffix>`.

#### VirulenceFinder database options

Pass one or more database names (space-separated, quoted) to `--vfinder`:

```
listeria, s.aureus_exoenzyme, s.aureus_hostimm, s.aureus_toxin,
stx, virulence_ecoli, virulence_ent, virulence_entfm
```

The results will be collected in the `results/` directory by default.
Change this with `--outDir`.


### Tools Included in the Pipeline

- **StarAMR**: Antimicrobial resistance detection using ResFinder, PlasmidFinder, and MLST.
- **MOB-suite**: Plasmid typing and reconstruction.
- **RGI**: AMR gene detection using the CARD database.
- **ECTyper**: *E. coli* serotyping (auto-triggered when *E. coli* is detected).
- **SISTR**: *Salmonella* serotyping (auto-triggered when *Salmonella* is detected).
- **Kleborate**: *Klebsiella* virulence and resistance profiling (auto-triggered for KPSC/KOSC).
- **Prokka**: Genome annotation (generates GFF, GBK, protein, gene, and genome FASTA).
- **Refseq_masher**: Species identification used to drive species-specific tool activation.
- **VirulenceFinder**: Detection of virulence genes (auto-triggered for *E. coli*).
- **Abricate (VFDB)**: Screening contigs against the VFDB virulence gene database.
- **Abricate (BacMet2)**: Screening contigs against BacMet2 confirmed and predicted
  biocide/metal resistance gene databases.
- **ICEberg**: Identification of integrative and conjugative elements (ICEs) via both
  protein (DIAMOND blastp) and nucleotide (blastn) searches.
- **PHASTER**: Detection of phage-related proteins in annotated genes via DIAMOND blastp
  against the PHAST protein database.
- **IntegronFinder**: Identification of integrons in genomes.
- **IslandPath-DIMOB**: Prediction of genomic islands using sequence composition and mobility genes.
- **DIGIS**: Detection of insertion sequences (IS elements) using a curated IS database.
