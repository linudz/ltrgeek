#! /usr/bin/env nextflow

/*

         _   _______ _____   ____   ____                  _    
        | | |__   __|  __ \ / /  \ /  \ \                | |   
        | |    | |  | |__) | | () | () | | __ _  ___  ___| | __
        | |    | |  |  _  /| |\__/ \__/| |/ _` |/ _ \/ _ \ |/ /
        | |____| |  | | \ \| |         | | (_| |  __/  __/   < 
        |______|_|  |_|  \_\ |         | |\__, |\___|\___|_|\_\
                            \_\       /_/  __/ |               
                                          |___/                

        LTR GEEK CORE PIPELINE SCRIPT.
        This is the core pipeline script for LTRgeek workflow.

        [SUMMARY]

        Part 1: The LTR Harvest launch and workaraund
        ----------------------------------------------------
            1.1: Launch suffixerator. Creates an index that is needed for LTRharvest run.
            1.2: Launch LTRharvest with the parameters defined in the config file.

        Part 2: The Filtering
        ----------------------------------------------------

        Part 3: The caracterization of the main features (coding regions, classification, etc)
        Part 4: The output processing

*/


// Part 1: The LTR harvest launch and workaraund

process harvestlaunch {

    publishDir 'harvest_process_publishdir'
    when:
    params.run_harvest = 'yes'

    """
    # 1.1: Launch suffixerator. Creates an index that is needed for LTRharvest run

    cp $params.in .
    gt suffixerator -db $params.in -indexname $params.in -tis -suf -lcp -des -ssp -sds -dna

    # 1.2: Launch LTRharvest with the parameters defined in the config file

    echo "gt ltrharvest \
    -index $params.in \
    -seed $harvest_seed  \
    -minlenltr $harvest_minlenltr \
    -maxlenltr $harvest_maxlenltr \
    -mindistltr $harvest_mindistltr \
    -maxdistltr $harvest_maxdistltr \
    -similar $harvest_similar \
    -mintsd $harvest_mintsd \
    -maxtsd $harvest_maxtsd \
    -vic $harvest_vic \
    -overlaps $harvest_overlaps \
    -xdrop $harvest_xdrop \
    -mat $harvest_mat \
    -mis $harvest_mis \
    -ins $harvest_ins \
    -del $harvest_del "

    """
}
