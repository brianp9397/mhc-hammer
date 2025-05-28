#!/usr/bin/env nextflow
/*
========================================================================================
    mhc_hammer
========================================================================================
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

//WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { MHC_HAMMER } from './workflows/mhc_hammer'
include { convert_alleles_to_mhc_hammer_format; mhc_hammer_solutions; mhc_hammer_inventory } from './mhc_hammer_prep.nf'

//
// WORKFLOW: Run main nf-mhc_hammer analysis pipeline
//
workflow NF_MHC_HAMMER {
    convert_alleles_to_mhc_hammer_format( alleles, hla_ref )
    mhc_hammer_solutions( samples.map{ [it[0], it[1], it[2], it[5]] }.join(purity, by:[0,1,2]) )
    mhc_hammer_inventory( samples.join(convert_alleles_to_mhc_hammer_format.out.mhc_hammer_alleles, by:[0,1,2]).join(mhc_hammer_solutions.out.formatted_solutions, by:[0,1,2]) )
    MHC_HAMMER ( mhc_hammer_inventory.out.inventory, params.mhc_coords, params.kmer_file, params.mhc_gtf, params.mhc_fasta, params.mhc_transcriptome_fasta, params.codon_table, params.contigs_placeholder )
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
//
workflow {
    NF_MHC_HAMMER ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
