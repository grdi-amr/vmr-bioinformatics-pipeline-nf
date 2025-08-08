// Pipeline input parameters
params.contigs = ""

// Database Locations
params.mobDB = "$workDir/databases/mobDB"
params.card_json = "$workDir/databases/card.json"
params.virulencefinderDB = "$workDir/databases/virulencefinder_db"
params.icebergDB = "$workDir/databases/iceberg"
params.prokkaDB = "$workDir/databases/prokka_db"


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
// Parameters to activate kleborate
//params.kleborate = false

// Help
params.help = false

//def staticKleborate = params.kleborate

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
    stub:
    """
    touch out/results.xlsx \
      out/detailed_summary.tsv \
      out/resfinder.tsv \
      out/plasmidfinder.tsv \
      out/mlst.tsv \
      out/settings.txt \
      out/pointfinder.tsv
    """
}
process run_ectyper {
    label "ECTYPER"
    publishDir "${params.outDir}/$sample/ECTYPER"

    input:
    tuple val(sample), path(contigs), val(flag)
    output:
    tuple val(sample), path("out/output.tsv"), emit: ectyper_tsv    
    script:
    """
    ectyper -i $contigs -o out
    """
    stub:
    """
    mkdir -p out
    echo -e "Isolate\tSerotype" > out/output.tsv
    echo -e "${sample}\tO157:H7" >> out/output.tsv
    """
}
process run_kleborate {
    label "KLEBORATE"
    publishDir "${params.outDir}/$sample/KLEBORATE"

    input:
    tuple val(sample), path(contigs), val(flag)
    
    output:
    tuple val(sample), path("out/*"), emit: keblorate_txt
    script:
    """
    kleborate -a $contigs -o out -p "$flag"
    """
    stub:
    """
    mkdir -p out
    echo "This is a stub result for sample: $sample" > out/${sample}_kleborate.txt
    """
}
process run_prokka{
    label "PROKKA"
    publishDir "${params.outDir}/$sample/PROKKA"
    cpus params.num_threads
    
    input:
    tuple val(sample), path(contigs)
    output:
    tuple val(sample), path("annotation/${sample}.gff"), emit: gff
    tuple val(sample), path("annotation/${sample}.gbk"), emit: gbk
    tuple val(sample), path("annotation/${sample}.fna"), emit: genome
    tuple val(sample), path("annotation/${sample}.faa"), emit: proteins
    tuple val(sample), path("annotation/${sample}.ffn"), emit: genes
    script:
    """
    #!/usr/bin/env bash
     
    dbs_dir=\$(prokka --listdb 2>&1 >/dev/null |  grep "databases in" | cut -f 4 -d ":" | tr -d " ") ;
    cp -r \$dbs_dir prokka_db
    cp ${params.prokkaDB}/PGAP_NCBI.hmm prokka_db/hmm
    ( cd  prokka_db/hmm/ ; for i in *.hmm ; do hmmpress -f \$i ; done )
    awk '{ if (\$0 ~ />/) print substr(\$0,1,21) ; else print \$0 }' $contigs > cleaned_header.fasta
    prokka \\
        --dbdir prokka_db \\
        --outdir annotation \\
        --cpus $task.cpus \\
        --mincontiglen 200 \\
        --prefix ${sample} \\
        --genus '' \\
        --species '' \\
        --strain \"${sample}\" \\
        cleaned_header.fasta
    rm -r prokka_db
    """
    stub:
    """
    mkdir -p annotation
    echo "##gff-stub" > annotation/${sample}.gff
    echo "##gbk-stub" > annotation/${sample}.gbk
    echo ">genome-sequence-stub" > annotation/${sample}.fna
    echo ">protein-sequence-stub" > annotation/${sample}.faa
    echo ">gene-sequence-stub" > annotation/${sample}.ffn
    """







}
process run_refseq_masher{
    label "REFSEQ_MASHER"
    publishDir "${params.outDir}/$sample/REFSEQ_MASHER"
    input:
    tuple val(sample), path(genome)

    output:
    tuple val(sample), path("refseq_masher_results.txt"), emit: results
    tuple val(sample), path("kleborate_flag.txt"), emit: flag

    script:
    """

    # Run tool
    refseq_masher \\
      -vv matches \\
      --top-n-results 10 \\
      --output-type tab \\
      $genome > refseq_masher_results.txt
    top_hit=\$(tail -n +2 "refseq_masher_results.txt" | head -n 1 | cut -f2)
    if echo "\$top_hit" | grep -E "Klebsiella pneumoniae|Klebsiella quasipneumoniae|Klebsiella variicola"; then
        echo "kpsc" > kleborate_flag.txt
    elif echo "\$top_hit" | grep -E "Klebsiella oxytoca|Klebsiella michiganensis|Klebsiella grimontii|Klebsiella pasteurii|Klebsiella huaxiensis"; then
         echo "kosc" > kleborate_flag.txt
    elif echo "\$top_hit" | grep -E "Escherichia coli"; then
         echo "virulence_ecoli" > kleborate_flag.txt
    elif echo "\$top_hit" | grep -E "Salmonella enterica"; then
         echo "Salmonella enterica" > kleborate_flag.txt
    elif echo "\$top_hit" | grep -E "Salmonella bongori"; then
         echo "Salmonella bongori" > kleborate_flag.txt
    else
        echo "none" > kleborate_flag.txt
    fi

    """
    stub:
    """
    # Create dummy refseq_masher_results.txt with header + one line
    if [[ "$(basename $genome)" == SRR26893262* ]]; then
        echo -e "rank\\tname\\tscore" > refseq_masher_results.txt
        echo -e "1\\tSalmonella enterica\\t100" >> refseq_masher_results.txt
        echo "Salmonella enterica" > kleborate_flag.txt
    elif [[ "\$filename" == SRR32114460* ]]; then
        echo -e "rank\\tname\\tscore" > refseq_masher_results.txt
        echo -e "1\\tKlebsiella pneumoniae\\t100" >> refseq_masher_results.txt
        echo "kpsc" > kleborate_flag.txt
    elif [[ "\$filename" == SRR32428968* ]]; then
        echo -e "rank\\tname\\tscore" > refseq_masher_results.txt
        echo -e "1\\tEscherichia coli\\t100" >> refseq_masher_results.txt
        echo "virulence_ecoli" > kleborate_flag.txt
    elif [[ "\$filename" == SRR24127505* ]]; then
        echo -e "rank\\tname\\tscore" > refseq_masher_results.txt
        echo -e "1\\tKlebsiella oxytoca\\t100" >> refseq_masher_results.txt
        echo "kosc" > kleborate_flag.txt
    else
        echo -e "rank\\tname\\tscore" > refseq_masher_results.txt
        echo -e "1\\tEnterobacter\\t100" >> refseq_masher_results.txt
        echo "none" > kleborate_flag.txt
    """
    

}

process run_sistr {
    label "SISTR"
    publishDir "${params.outDir}/$sample/SISTR"

    input:
    tuple val(sample), path(contigs), val(flag)
    output:
    tuple val(sample), path("sistr.tab"), emit: sistr_tsv
    script:
    """
    sistr -i $contigs "$flag" -o sistr -f tab
    """
    stub:
    """
    echo -e "genome\tserovar\tantigen\tcgmlst_ST\n\$sample\tEnteritidis\t1,9,12:g,m:-\t1234" > sistr.tab
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
    stub:
    """
    echo -e "SEQUENCE\tGENE\t%IDENTITY\t%COVERAGE\tDATABASE" > amr.vfdb.results.tsv
    echo -e "$sample\tmock_gene\t99.8\t100.0\tvfdb" >> amr.vfdb.results.tsv
    """


}

process run_virulencefinder {
    label "VFINDER"
    publishDir "${params.outDir}/$sample/VIRULENCEFINDER"

    input:
    tuple val(sample), path(contigs), val(flag)
    output:
    tuple val(sample), path("out/"), emit: vf_tsv
    
    script:
    """
    mkdir -p out
    virulencefinder.py -i "$contigs" -o "out" -d "$flag" -p $params.virulencefinderDB
    """
    stub:
    """
    mkdir -p out
    echo -e "Isolate ID\tGene\t%Identity\t%Coverage\tDatabase" > out/mock_virulencefinder.tsv
    echo -e "$sample\tmock_gene\t99.5\t98.7\t$flag" >> out/mock_virulencefinder.tsv
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
    stub:
    """
    mkdir -p localDB
    echo "stubbed" > localDB/dummy.txt
    """

}

process run_iceberg {
    label = [ 'misc', 'process_low' ]
    publishDir "${params.outDir}/$sample/ICEBERG"
    cpus params.num_threads

    input:
    tuple val(sample), file(genes_aa)
    tuple val(sample), file (genome)
    output:
    tuple val(sample), path("${sample}_iceberg_blastp_onGenes.summary.txt") , emit: genes_summary
    tuple val(sample), path("${sample}_iceberg_blastp_onGenes.txt")         , emit: results
    tuple val(sample), path("${sample}_iceberg_blastn_onGenome.summary.txt"), emit: genome_summary
    tuple val(sample), path('*.txt'), emit: all

    script:
    """
    # ICEberg is a protein and nucleotide dabatase
    # In protein are the genes found inside ICEs
    # In nucleotide are the full-length ICEs
     python $projectDir/bin/run_blasts.py \\
        blastp \\
        --query $genes_aa \\
        --db ${params.icebergDB}/diamond.dmnd \\
        --minid 85 \\
        --mincov 85 \\
        --threads $task.cpus \\
        --out ${sample}_iceberg_blastp_onGenes.txt --2way | \\
    sed -e 's/GENE/ICEBERG_ID/g' > ${sample}_iceberg_blastp_onGenes.summary.txt ;
   
    ## Checking for full-length ICEs
    ### The blast db was throwing errors
    makeblastdb \\
        -dbtype nucl \\
        -in ${params.icebergDB}/sequences \\
        -out sequences ;
     python $projectDir/bin/run_blasts.py \\
        blastn \\
        --query $genome \\
        --db sequences \\
        --minid 0 \\
        --mincov 0 \\
        --threads ${task.cpus} \\
        --out ${sample}_iceberg_blastn_onGenome.txt | \\
    sed -e 's/GENE/ICEBERG_ID/g' > ${sample}_iceberg_blastn_onGenome.summary.txt ;
    """
    stub:
    """
    ouch ${sample}_iceberg_blastp_onGenes.txt
    echo -e "ICEBERG_ID\\tidentity\\tcoverage" > ${sample}_iceberg_blastp_onGenes.summary.txt

    touch ${sample}_iceberg_blastn_onGenome.txt
    echo -e "ICEBERG_ID\\tidentity\\tcoverage" > ${sample}_iceberg_blastn_onGenome.summary.txt
    """    


}
process run_integron_finder{
    label "INTEGRON_FINDER"
    publishDir "${params.outDir}/$sample/INTEGRON_FINDER"
    cpus params.num_threads

    input:
    tuple val(sample), file(genome)

    output:
    tuple val(sample), path("*")                      , emit: all
    tuple val(sample), path("${sample}_integrons.gbk"), emit: gbk, optional: true
    path("integronfinder_version.txt")

    script:
    """
    # Get version
    integron_finder --version > integronfinder_version.txt ;

    # run tool
    integron_finder \\
        --local-max \\
        --func-annot \\
        --pdf \\
        --gbk \\
        --cpu $task.cpus \\
        $genome
    
    # move results
    mv Results_Integron_Finder_${sample}/* . ;
    rm -rf Results_Integron_Finder_${sample} ;
    
    # convert to gff if available
    for gbk in \$(ls *.gbk) ; do
        cat \$gbk >> ${sample}_integrons.gbk ;
    done
    """
    stub:
    """
    mkdir -p Results_Integron_Finder_${sample}
    touch Results_Integron_Finder_${sample}/dummy_output.txt
    touch ${sample}_integrons.gbk
    echo 'IntegronFinder vX.X.X' > integronfinder_version.txt
    """   




}
process run_island_path{
    label "ISLANDPATH"
    publishDir "${params.outDir}/$sample/ISLAND_PATH"
    cpus params.num_threads

    input:
    tuple val(sample), file("annotation.gbk")

    output:
    tuple val(sample), path("${sample}_predicted_GIs.bed"), emit: results
    script:
    """
    # Split genbank files
    splitgenbank.pl annotation.gbk && rm annotation.gbk ;

    # Run islandpath in each
    touch ${sample}_predicted_GIs.bed ;
    for file in \$(ls *.gbk); do \
      touch \${file%%.gbk}_GIs.txt ;
      ( sed '/CDS.*::.*0/d' \$file | grep -q "CDS" ) && \\
          islandpath \\
          \$file \\
          \${file%%.gbk}_GIs.txt 2> dimob.err ;
          name="\${file%%.gbk}" ;
          awk -v contig=\$name 'BEGIN { FS = "\\t"; OFS="\\t" } { print contig,\$2,\$3 }' \${file%%.gbk}_GIs.txt >> ${sample}_predicted_GIs.bed ;
        done
    """
    stub:
    """
    touch ${sample}_predicted_GIs.bed
    echo -e "${sample}\t0\t1000" > ${sample}_predicted_GIs.bed
    """


}
process run_digis{
    label = [ 'misc', 'process_low' ]
    publishDir "${params.outDir}/$sample/DIGIS"
    cpus params.num_threads
    
    input:
    tuple val(sample), path(genome), path(genbank)
     
    output:
    tuple val(sample), path("digIS")                      , emit: all
    tuple val(sample), path("digIS/results/${sample}.gff"), emit: gff
    tuple val(sample), path("${sample}_IS.gff"), path("digIS/results/fastas/${sample}_IS.fa"), path("digIS/results/fastas/${sample}_IS.faa"), emit: gff_and_sequences

    script:
    """
    # run digIS
    conda run -n digIS python3 \$(which digIS_search.py) -i $genome -g $genbank -o digIS

    # parse digIS to get nucleotide and aminoacide
    # also put ids in uppercase
    # required for annotation merging and sqldb

   ## dir for fastas
   mkdir -p digIS/results/fastas ;

   ## save info in gff
   sed \\
     -e 's/id=/ID=/g' \\
     digIS/results/${sample}.gff > ${sample}_IS.gff ;

   ## get nucl sequences
   gff-toolbox \\
     convert \\
     -i ${sample}_IS.gff  \\
     -f fasta-nt \\
     --fasta $genome \\
     --fasta_features transposable_element > digIS/results/fastas/${sample}_IS.fa  ;
  
   ## get prot sequences
   gff-toolbox \\
     convert \\
     -i ${sample}_IS.gff  \\
     -f fasta-aa \\
     --fasta $genome \\
     --fasta_features transposable_element > digIS/results/fastas/${sample}_IS.faa ;
   """
   stub:
   """
   mkdir -p digIS/results/fastas
   touch digIS/results/${sample}.gff
   cp digIS/results/${sample}.gff ${sample}_IS.gff
   touch digIS/results/fastas/${sample}_IS.fa
   touch digIS/results/fastas/${sample}_IS.faa
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
    def KLEBORATE_RESULTS = Channel.empty()
    def VFINDER_RESULTS = Channel.empty()
    //def staticKleborate = params.kleborate    
    // Get the Contigs into a channel
    CONTIGS = Channel
                .fromPath(params.contigs)
                .map { file -> tuple(file.baseName, file) }
    //Run_Prokka
    PROKKA_RESULTS = run_prokka(CONTIGS)
    MASHER_RESULTS = run_refseq_masher (PROKKA_RESULTS.genome)
    
    //check_kleborate_species(MASHER_RESULTS.results)
    //    .set { kleborate_flag_ch }

    //kleborate_flag_ch.view { println "kleborate_flag_ch: $it" }
    //KLEBORATE_RUN = SPECIE.kleborate_flag_ch.filter { it == 'kpsc' || it == 'kosc' }
    // Run kleborate only if the flag says so
    run_kleborate_flag = MASHER_RESULTS.flag
    .map { sample, flag_file -> 
        def flag = flag_file.text.trim()
        //def valid = (flag == 'kpsc' || flag == 'kosc') ? flag : false
        //tuple(sample, valid)
        (flag == 'kpsc' || flag == 'kosc') ? tuple(sample, flag) : null
    }
    .filter { it != null }
    kleborate_input = PROKKA_RESULTS.genome
    .join(run_kleborate_flag)
    kleborate_input.view()
    KLEBORATE_RESULTS = run_kleborate(kleborate_input)
    run_ectyper_flag = MASHER_RESULTS.flag
    .map { sample, flag_file ->
        def flag = flag_file.text.trim()
        (flag == 'virulence_ecoli') ? tuple(sample, flag) : null
    }
    .filter { it != null }
    ectyper_input = PROKKA_RESULTS.genome
    .join(run_ectyper_flag)
    ECTYPER_RESULTS = run_ectyper(ectyper_input)
    run_salmonella_flag = MASHER_RESULTS.flag
    .map { sample, flag_file ->
        def flag = flag_file.text.trim()
        (flag == 'Salmonella enterica' || flag == 'Salmonella bongori') ? tuple(sample, flag) : null
    }
    .filter { it != null }
    salmonella_input = PROKKA_RESULTS.genome
    .join(run_salmonella_flag)
    SISTR_RESULTS = run_sistr(salmonella_input)
    run_vfinder_flag = MASHER_RESULTS.flag
    .map { sample, flag_file ->
        def flag = flag_file.text.trim()
        (flag == 'virulence_ecoli' || flag == 'virulence_ent') ? tuple(sample, flag) : null
    }
    .filter { it != null }
    vfinder_input = PROKKA_RESULTS.genome
    .join(run_vfinder_flag)
   vfinder_input.view() 
   VFINDER_RESULTS = run_virulencefinder(vfinder_input) 
   if ( params.ectyper == true) {
       ECTYPER_RESULTS = run_ectyper(PROKKA_RESULTS.genome)
    }
       
    if (params.test) {
    log.info "Test mode enabled. Exiting pipeline early for CI testing."
    System.exit(0)
    }    
    //if (params.sistr == true) {
    //   SISTR_RESULTS = run_sistr(PROKKA_RESULTS.genome) 
    //}
    //if (params.kleborate == 'kpsc' || params.kleborate == 'kosc'){
    //   KLEBORATE_RESULTS = run_kleborate(PROKKA_RESULTS.genome)
    //}
    if (params.vfinder){
       VFINDER_RESULTS = run_virulencefinder(PROKKA_RESULTS.genome)
    }
    //Run Iceberg
    ICEBERG_RESULTS = run_iceberg(PROKKA_RESULTS.proteins, PROKKA_RESULTS.genome)
    INTEGRON_FINDER_RESULTS =run_integron_finder(PROKKA_RESULTS.genome)
    ISLAND_PATH_RESULTS = run_island_path(PROKKA_RESULTS.gbk)
    DIGIS_RESULTS = run_digis(PROKKA_RESULTS.genome.join(PROKKA_RESULTS.gbk))
   // Run mob_recon on the contigs.
    MOB_RESULTS = run_mobSuite(PROKKA_RESULTS.genome)
    // Run star_amr
    STARAMR_RESULTS=runStarAMR(PROKKA_RESULTS.genome)
    // Run Abricate
    ABRICATE_RESULTS=run_abricate(PROKKA_RESULTS.genome)	    
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
        RGI_RESULTS = run_RGI(PROKKA_RESULTS.genome, LOCAL_DB.out.collect())
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
    //TOTAL_JSON = json_generator(CAT_TAB)
   
    // Create report
//    create_report(CAT_TAB)
} 

