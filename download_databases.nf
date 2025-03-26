// Parameters
// Database Locations
params.mobDB = "$workDir/databases/mobDB"
params.card_json = "$workDir/databases"
params.virulencefinderDB = "$workDir/databases/virlencefinder_db"
params.prokkaDB = "$workDir/databases/prokka_db"

// Database download controls
params.download_mobDB = false
params.download_card_json = false
params.download_vfDB = false
params.download_prokkaDB = false
params.download_all = false
params.overwrite = false

// Set a flag to donwload the mob database based on input paramters
if ( params.download_all | params.download_mobDB ){
    DL_mob = true 
} else { 
    DL_mob = false
}

// Set a flag to donwload the card json based on input paramters
if ( params.download_all | params.download_card_json ){
    DL_card = true 
} else { 
    DL_card = false
}
// Set a flag to donwload the virulencefinder db based on input paramters
if ( params.download_all | params.download_vfDB ){
    DL_vfdb = true
} else {
    DL_vfdb = false
}
if ( params.download_all | params.download_prokkaDB ){
    DL_prokka = true
} else {
    DL_prokka = false
}

// Log message
log.info """ 

    Download Databases
    -------------------
    MOB-suite databases: ${DL_mob}
    CARD json: ${DL_card}
    Virulencefinder database: ${DL_vfdb}
    Prokka database: ${DL_prokka}
""".stripIndent(true)


process download_CARD_json {
    label "RGI"
    publishDir params.card_json, mode: 'move', overwrite: true

    input: 

    output:
    stdout emit: stdout
    path 'card.json', emit: json

    script:
    """
    echo "Downloading CARD json"
    wget https://card.mcmaster.ca/latest/data
    tar -xvf data ./card.json
    """
    stub: 
    """
    echo Data! > ./card.json
    """

}
process download_virulencefinder_database {

    publishDir params.card_json, mode: 'move', overwrite: true
    input:

    output:
    path "virulencefinder_db/*"
       
    """
    git clone https://bitbucket.org/genomicepidemiology/virulencefinder_db/
    """
}

process download_MOB_database {
    label "MOB"
    publishDir params.mobDB, mode: 'move', overwrite: true

    input:

    output: 
    path "*"

    script:
    """ 
    mob_init --database_directory '.'
    """
    stub:
    """
    touch status.txt
    """
}

process download_prokka_database {
    label = ['db_download','process_low']
    publishDir params.prokkaDB, mode: 'move', overwrite: true

    input:

    output:
    file "*"

    script:
    """
    # download prokka additional database
    wget --tries=10 https://ftp.ncbi.nlm.nih.gov/hmm/TIGRFAMs/release_15.0/TIGRFAMs_15.0_HMM.LIB.gz
    gzip -d TIGRFAMs_15.0_HMM.LIB.gz
    mv TIGRFAMs_15.0_HMM.LIB TIGRFAMs_15.0.hmm
    wget --tries=10 https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.LIB -O PGAP_NCBI.hmm
    """
    stub:
    """
    touch status.txt
    """
}

workflow download_CARD_WF {
        
    // If the JSON does not exist, simply download it: 
    if ( !file(params.card_json).exists() ){
        println "Downloading CARD json file, please wait"
        download_CARD_json()
    // But, if it DOES exist:
    } else {
        // Has overwrite flag been set?
        // If YES: download it
        if ( params.overwrite ){
            println "WARNING: file at $params.card_json detected, overwriting!"
            download_CARD_json()
        // But if NO: skip
        } else {
            println(
                "WARNING: file at $params.card_json detected, " + 
                "but overwrite set to $params.overwrite, skipping")
        }
    }
}


workflow download_MOB_WF {
    
    mobDB_dir = file(params.mobDB)
    // If the MOB database does not exist, simply download it: 
    if ( !mobDB_dir.exists() ){
        println "Downloading MOB-suite databases, please wait"
        download_MOB_database()
    // But, if it DOES exist:
    } else {
        // Is the directory empty?
        if ( mobDB_dir.isEmpty() ){
        // If EMPTY: download it:
            println "Downloading MOB-suite databases, please wait"
            download_MOB_database()
        // But if the directory is not empty:
        } else { 
            // Has overwrite flag been set?
            // If YES: download it
            if ( params.overwrite ){
                println (
                    "WARNING: Contents of $params.mobDB detected, " +
                    "overwriting!")
                download_MOB_database()
            // But if NO: skip
            } else {
                println(
                    "WARNING: Contents of $params.mobDB detected " +
                    "but overwrite set to $params.overwrite, skipping" )
            }       
        }
    }

}

workflow download_PROKKA_WF{
    prokkaDB_dir = file (params.prokkaDB)
    if ( !prokkaDB_dir.exists() ){
        println "Downloading Prokka database, please wait"
        download_prokka_database()
    // But, if it DOES exist:
    } else {
        // Is the directory empty?
        if ( prokkaDB_dir.isEmpty() ){
        // If EMPTY: download it:
            println "Downloading  virulence finder database, please wait"
            download_prokka_database()
        // But if the directory is not empty:
        } else {
            // Has overwrite flag been set?
            // If YES: download it
            if ( params.overwrite ){
                println (
                    "WARNING: Contents of $params.prokkaDB detected, " +
                    "overwriting!")
                download_prokka_database()
            // But if NO: skip
            } else {
                println(
                    "WARNING: Contents of $params.prokkaDB detected " +
                    "but overwrite set to $params.overwrite, skipping" )
            }
        }
    }


}

workflow download_VF_WF {

    vfDB_dir = file(params.virulencefinderDB)
    // If the Virulencefinder database does not exist, simply download it:
    if ( !vfDB_dir.exists() ){
        println "Downloading Virulencefinder database, please wait"
        download_virulencefinder_database()
    // But, if it DOES exist:
    } else {
        // Is the directory empty?
        if ( vfDB_dir.isEmpty() ){
        // If EMPTY: download it:
            println "Downloading  virulence finder database, please wait"
            download_virulencefinder_database()
        // But if the directory is not empty:
        } else {
            // Has overwrite flag been set?
            // If YES: download it
            if ( params.overwrite ){
                println (
                    "WARNING: Contents of $params.virulencefinderDB detected, " +
                    "overwriting!")
                download_virulencefinder_database()
            // But if NO: skip
            } else {
                println(
                    "WARNING: Contents of $params.virulencefinderDB detected " +
                    "but overwrite set to $params.overwrite, skipping" )
            }
        }
    }

}

workflow {

    // Has user specified that the CARD database be downloaded?
    if ( DL_card ){
        download_CARD_WF()
    }
    // Has use specifed that the MOB database be downloaded?
    if ( DL_mob ){
        download_MOB_WF()
    }
    // Has use specifed that the MOB database be downloaded?
    if ( DL_vfdb ){
        download_VF_WF()
    }
    if ( DL_prokka){
        download_PROKKA_WF()
    }
}



