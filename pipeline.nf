// Pipeline input parameters
params.contigs = ""

// Database Locations
params.mobDB = "$workDir/databases/mobDB"
params.card_json = "$workDir/databases/card.json"
params.virulencefinderDB = "$workDir/databases/virulencefinder_db"

// Full RGI or on plasmids only?
params.plasmids_only = false

// Output for results
params.outDir = "$launchDir/results"

// Process parameters
params.num_threads = 1

// Parameters for StarAMR
params.species = "Escherichia coli"
// Parameters to activate ectyper
params.ectyper = false
// Parameters to activate sistr
params.sistr = false
// Parameters for vfinder database 
params.vfinder = ""

// Help
params.help = false

def helpMessage() { 

    log.info """
        Usage: 
        nextflow run pipeline.nf \\
            -profile [conda|docker] \\
            --contigs "path/to/contigs/*.fasta"

        Options:
          --contigs         PATH to contigs to run the pipeline on. Multiple 
                            files can be specified using bash globs.
          --mobDB           PATH to the DIRECTORY of the MOB-suite databases
          --card_json       PATH to CARD's card.json file.
          --plasmids_only   Run RGI only on contigs identified as plasmids, thus
                            ignoring any chromosomal AMR determinants.
          --outDir          DIRECTORY to link results to.
          --num_threads     Number of threads to use in downstream processes, 
                            per sample.
	  --species	    Species for starAMR samplesheet.
          --ectyper         if used it will  serotype E.coly .Default=false 
          --sistr           if used it will serotype Salmonella.If used make sure to 
                            use parameter 'species'. Default=false
          --vfinder         Use one or more of vfinder available databases: listeria, 
                            s.aureus_exoenzyme, s.aureus_hostimm, s.aureus_toxin, 
                            stx, virulence_ecoli, virulence_ent, virulence_entfm.
    """.stripIndent(true)
}

if (params.help) { 
    helpMessage() 
    exit 0
}


// Log
log.info """\
    Mob-Suite RGI pipeline
    ===================================
    Profile     : ${workflow.profile}
    Contigs     : ${params.contigs}
    Species     : ${params.species}
    Sistr       : ${params.sistr}
    Ectyper     : ${params.ectyper}
    mobDB       : ${params.mobDB}
    vfDB        : ${params.virulencefinderDB}
    card.json   : ${params.card_json}
    outDir      : ${params.outDir} 
   
    PlasmidsOnly? ${params.plasmids_only}
    """
   .stripIndent(true)

// Note: what happens if these files are not generated?, e.g., with
// optional flag?


process runStarAMR {
    label "STARAMR"
    publishDir "${params.outDir}/$sample/STARAMR"

    input:
    tuple val(sample), path(contigs)
    

    output:
    tuple val(sample), path("out/results.xlsx")          , emit: results_xlsx
    tuple val(sample), path("out/detailed_summary.tsv")  , emit: summary_tsv
    tuple val(sample), path("out/resfinder.tsv")         , emit: resfinder_tsv
    tuple val(sample), path("out/plasmidfinder.tsv")     , emit: plasmidfinder_tsv
    tuple val(sample), path("out/mlst.tsv")              , emit: mlst_tsv
    tuple val(sample), path("out/settings.txt")          , emit: settings_txt
    tuple val(sample), path("out/pointfinder.tsv")       , emit: pointfinder_tsv, optional: true 

    script:
    """
    if [ -d "out" ]; then
	    rm -r "out"
            echo "Directory 'out' removed."
    fi
    staramr search -o out  $contigs  
    """
}
process run_ectyper {
    label "ECTYPER"
    publishDir "${params.outDir}/$sample/ECTYPER"

    input:
    tuple val(sample), path(contigs)
    output:
    tuple val(sample), path("out/output.tsv"), emit: ectyper_tsv    
    script:
    """
    ectyper -i $contigs -o out
    """
}
process run_sistr {
    label "SISTR"
    publishDir "${params.outDir}/$sample/SISTR"

    input:
    tuple val(sample), path(contigs)
    output:
    tuple val(sample), path("sistr.tab"), emit: sistr_tsv
    script:
    """
    sistr -i $contigs "$params.species" -o sistr -f tab
    """
}
process run_abricate {
    label "ABRICATE"
    publishDir "${params.outDir}/$sample/ABRICATE"

    input:
    tuple val(sample), path(contigs)
    output:
    tuple val(sample), path("amr.vfdb.results.tsv"), emit: vfdb_tsv
    script:
    """
    abricate $contigs --db vfdb > amr.vfdb.results.tsv
    """


}

process run_virulencefinder {
    label "VFINDER"
    publishDir "${params.outDir}/$sample/VIRULENCEFINDER"

    input:
    tuple val(sample), path(contigs)
    output:
    tuple val(sample), path("out/"), emit: vf_tsv
    
    script:
    """
    mkdir -p out
    virulencefinder.py -i "$contigs" -o "out" -d "$params.vfinder" -p $params.virulencefinderDB
    """


}
process load_RGI_database { 
    label "RGI" 

    input:
    path card_json

    output:
    path 'localDB', type: "dir", emit: out

    script:
    """

    rgi load --local --card_json $card_json 

    """

}

process run_RGI { 
    label "RGI"
    publishDir "${params.outDir}/$sample/RGI"
    cpus params.num_threads

    input:
    tuple val(sample), path(contigs)
    path local_DB

    output:
    tuple val(sample), path("rgi_results.txt"), emit: table
    tuple val(sample), path("rgi_results.json"), emit: json

    script: 
    """
    echo "plasmids only? $params.plasmids_only"
    rgi main \
        --local \
        --num_threads ${task.cpus} \
        --input_sequence $contigs \
        --output_file "rgi_results"

    """ 

    stub: 
    """
    touch rgi_results.txt
    touch rgi_results.json
    """
}

process concatenate_plasmid_seqs {

    input:
    tuple val(sample), path(plasmid_contigs)

    output:
    tuple val(sample), path("*.fasta"), emit: contigs

    script:
    """
    cat $plasmid_contigs > ${sample}.fasta

    """
}


process run_mobSuite {
    label "MOB"
    publishDir "${params.outDir}/$sample"
    cpus params.num_threads

    input: 
    tuple val(sample), path(contigs)

    output: 
    tuple val(sample), path("mobSuite/contig_report.txt"), 
      emit: contig_table
    tuple val(sample), path("mobSuite/plasmid*.fasta"), 
      emit: plasmid_fastas, optional: true
    tuple val(sample), path("mobSuite/mobtyper_results.txt"), 
      emit: typer, optional: true
    tuple val(sample), path("mobSuite/mge.report.txt"), 
      emit: mge, optional: true

    script:
    """

    mob_recon \
      --infile $contigs \
      --outdir mobSuite \
      --num_threads 1 \
      --database_directory $params.mobDB \
      --force 

    """
    stub: 
    """
    mkdir -p mobSuite
    touch mobSuite/contig_report.txt
    touch mobSuite/chromosome.fasta
    """
}

process merge_tables {
    label "RGI"
    publishDir "${params.outDir}/$sample/Merge"
    cache false

    input:
    tuple val(sample), path(tables)

    output:
    path('merged_tables.csv'), emit: out

    script: 
    """
    python $projectDir/bin/merge.py ${tables[0]} ${tables[1]}

    """
}
process json_generator {
    label "RGI"
    publishDir "${params.outDir}", mode: 'copy'
    cache 'none'

    // This process runs after all samples are processed, taking only the directory as an argument
    input:
    path tables

    output:
    path("results_summary.json"), emit: json_out
    script:
    """
    python $projectDir/bin/json_generator.py ${params.outDir}
    """
}

process create_report {
    label "RGI"
    publishDir "${params.outDir}"
    cache false

    input: 
    path(table)

    output: 
    path('report.html'), emit: out

    script: 
    """
    #!/usr/bin/env Rscript 

    rmarkdown::render(
        input = "${projectDir}/bin/generate_report.rmd",
        params = list(data = "${table}"),
        output_file = "report.html",
        knit_root_dir = getwd(),
        output_dir = getwd()
        )
    """
}

workflow {
    def ECTYPER_RESULTS = Channel.empty()
    def SISTR_RESULTS = Channel.empty()
    def VFINDER_RESULTS = Channel.empty()    
    // Get the Contigs into a channel
    CONTIGS = Channel
                .fromPath(params.contigs)
                .map { file -> tuple(file.baseName, file) }

    if ( params.ectyper == true) {
       ECTYPER_RESULTS = run_ectyper(CONTIGS)
    }
       
    
    if (params.sistr == true) {
       SISTR_RESULTS = run_sistr(CONTIGS) 
    }
    if (params.vfinder){
    VFINDER_RESULTS = run_virulencefinder(CONTIGS)
    }
    
   // Run mob_recon on the contigs.
    MOB_RESULTS = run_mobSuite(CONTIGS)
    // Run star_amr
    STARAMR_RESULTS=runStarAMR(CONTIGS)
    // Run Abricate
    ABRICATE_RESULTS=run_abricate(CONTIGS)	    
    // Get the CARD Json
    JSON = Channel.fromPath(params.card_json)
    // Load RGI database locally.
    LOCAL_DB = load_RGI_database(JSON)

    // Run RGI
    // Does the user want to run RGI on plasmids only?
    if ( params.plasmids_only ){
        // Merge plasmid seqs into a single file
        PLASMID_CONTIGS = concatenate_plasmid_seqs(MOB_RESULTS.plasmid_fastas)
        // Run RGI on the plasmid contigs only
        RGI_RESULTS = run_RGI(  PLASMID_CONTIGS.contigs,
                                LOCAL_DB.out.collect()  )
    } else { 
        // Otherwise, just run RGI on the full sample set
        RGI_RESULTS = run_RGI(CONTIGS, LOCAL_DB.out.collect())
    }

    // Create channel with tables to be combined
    TABLES = RGI_RESULTS.table 
                .concat(MOB_RESULTS.contig_table)
                .groupTuple(size: 2)

    // Merge the tables using an included script
    MERGE_TAB = merge_tables(TABLES)

    // This operation concatenates the CSV files, but leaves just one header at 
    // the top
 
    CAT_TAB = MERGE_TAB.out
                .collectFile(keepHeader: true, 
                             skip: 1, 
                             name: 'Mob_rgi_contig_results.csv', 
                             storeDir: params.outDir )
    TOTAL_JSON = json_generator(CAT_TAB)
   
    // Create report
//    create_report(CAT_TAB)
 }

