/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/
nextflow.enable.dsl=2


/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// SUBWORKFLOWS
//

include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { PREPROCESSING; PREPROCESSING_WITHOUT_BAM_SUBSETTING } from '../subworkflows/local/preprocessing_subworkflow'
include { DNA_ANALYSIS; DETECT_MUTATIONS } from '../subworkflows/local/mhc_hammer_dna_subworkflow'
include { GENERATE_DNA_NOVOALIGN_BAMS } from '../subworkflows/local/dna_novoalign_alignment_subworkflow'
include { GENERATE_RNA_NOVOALIGN_BAMS } from '../subworkflows/local/rna_novoalign_alignment_subworkflow'
include { GENERATE_RNA_STAR_BAMS } from '../subworkflows/local/rna_star_alignment_subworkflow'
include { ALT_SPLICING } from '../subworkflows/local/alt_splicing_subworkflow'
include { RNA_ANALYSIS as NOVOALIGN_RNA_ANALYSIS;
            RNA_ANALYSIS as STAR_RNA_ANALYSIS } from '../subworkflows/local/mhc_hammer_rna_subworkflow'

//
// MODULES 
//

include { HLAHD; HLAHD_LOCAL } from '../modules/local/run_hlahd'
include { GENERATE_REFERENCES } from '../modules/local/generate_references'
include { CREATE_TUMOUR_NORMAL_SPLICING_COHORT_TABLE; CREATE_LIBRARY_SIZE_TABLE;
        CREATE_MOSDEPTH_COHORT_TABLE; CREATE_MHC_HAMMER_TABLE; 
        CREATE_KNOWN_SPLICING_COHORT_TABLE; CREATE_NOVEL_SPLICING_COHORT_TABLE } from '../modules/local/make_overview_table'
include { filterByMetadataField; filterByMetadataFieldTwo } from '../lib/core_functions'

/*
========================================================================================
    IMPORT NF-CORE MODULE
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

// import custom functions
include { check_hlahd_results; check_tumour_results; 
        file_collector as cohort_known_file_collector;
        file_collector as cohort_novel_file_collector;
        file_collector as cohort_mosdepth_file_collector;
        file_collector as cohort_tumour_normal_file_collector;
        file_collector as library_size_file_collector;
        file_collector as allele_table_file_collector;
        file_collector as hlahd_gl_samples_file_collector } from '../lib/core_functions'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/


workflow MHC_HAMMER {
    take:
    ch_input
    mhc_coords_ch
    mhc_kmer_ch
    mhc_gtf_ch
    mhc_fasta_ch
    mhc_transcriptome_fasta_ch
    codon_table_ch
    contigs_file_ch

    main:
    ch_versions = Channel.empty()

    // Define empty channels here so the RNA analysis can be run independently of the DNA analysis
    cohort_mosdepth_input_ch = Channel.empty()
    rna_all_snps_ch = Channel.empty()
    rna_bam_read_count_ch = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //

    INPUT_CHECK ( ch_input, mhc_fasta_ch )


    // Add software used in INPUT_CHECK
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
 
    // These channels will be used to prevent pipeline stalling
    tumour_wxs_sample_count_ch = INPUT_CHECK.out.tumour_wxs_sample_count 

    rna_sample_count_ch = INPUT_CHECK.out.rna_sample_count
    
    germline_sample_count_ch = INPUT_CHECK.out.germline_sample_count

    normal_rna_count = INPUT_CHECK.out.normal_rna_sample_count

    gl_samples = INPUT_CHECK.out.gl_samples

    purity_ploidy_ch = INPUT_CHECK.out.purity_ploidy

    total_tumour_count = INPUT_CHECK.out.total_tumour_sample_count

    //
    // SUBWORKFLOW: Run Preprocessing steps - gererates bams / fqs containing putative HLA aligned reads, runs flagstat on input bams
    // 

    PREPROCESSING ( 
        INPUT_CHECK.out.bams_ch, 
        germline_sample_count_ch,
        mhc_coords_ch,
        contigs_file_ch,
        mhc_kmer_ch
    )

    // update downstream input channels
    preprocessing_versions = PREPROCESSING.out.versions
    hlahd_input_ch = PREPROCESSING.out.hlahd_input
    wxs_fqs_ch = PREPROCESSING.out.wxs_fqs
    dna_library_size_ch = PREPROCESSING.out.dna_library_size
    rna_fqs_ch = PREPROCESSING.out.rna_fqs
    rna_bams_ch = PREPROCESSING.out.rna_bams
    rna_library_size_ch = PREPROCESSING.out.rna_library_size


    // Add software used in PREPROCESSING
    ch_versions = ch_versions.mix(preprocessing_versions)

    //
    // MODULE: Run HLAHD on WXS normal sample
    //

    // define outputs
    //hlahd_genotype_ch = INPUT_CHECK.out.hla_alleles
    successful_hlahd_samples = INPUT_CHECK.out.hla_alleles

    // Now we check the output of HLAHD to see if it has failed for any samples
    // Failed samples will be excluded from the rest of the pipeline
    // Failed samples will have no typed HLA genes in their HLAHD output

    // Define branch criteria
    /*
    def criteria = branchCriteria {
        failed_hlahd_samples: check_hlahd_results(it)
        successful_hlahd_samples: !check_hlahd_results(it)
    }
    
    // Apply the branch criteria to the HLAHD genotype output
    hlahd_samples_by_status = hlahd_genotype_ch
        .branch(criteria)
    
    // output a warning if any samples failed HLAHD
    hlahd_samples_by_status.failed_hlahd_samples.view { 
    patient_id, file -> 
    log.warn("HLAHD failed to produce HLA allele estimates for sample ${patient_id}.\n" +
         "As a result, MHCHammer will not be executed for this sample.\n" +
         "The issue may be due to overly strict subset_bam process parameters. Consider allowing more reads to be included.")
    }
    */
    //hlahd_samples_by_status.successful_hlahd_samples
    // MODULE: generates sample specific reference files for the rest of the pipeline
    GENERATE_REFERENCES ( 
        successful_hlahd_samples, 
        mhc_fasta_ch, 
        mhc_gtf_ch, 
        mhc_transcriptome_fasta_ch 
    )

    ch_versions = ch_versions.mix(GENERATE_REFERENCES.out.versions.first())

    // SUBWORKFLOW: align wxs fqs to personalised reference and generate allele specific bams
    GENERATE_DNA_NOVOALIGN_BAMS (
        wxs_fqs_ch,
        GENERATE_REFERENCES.out.genome_reference,
        GENERATE_REFERENCES.out.mosdepth_reference
    )
    
    // ch_versions = ch_versions.mix(GENERATE_DNA_NOVOALIGN_BAMS.out.versions)

    // SUBWORKFLOW: detect HLA allele copy number and allelic imbalance
    
    DNA_ANALYSIS (
        dna_library_size_ch,
        GENERATE_REFERENCES.out.genome_snp_positions,
        GENERATE_REFERENCES.out.patient_gtf,
        purity_ploidy_ch,
        GENERATE_DNA_NOVOALIGN_BAMS.out.hla_allele_bams_ch,
        GENERATE_DNA_NOVOALIGN_BAMS.out.passed_heterozygous_hla_genes_ch
    )
    
    ch_versions = ch_versions.mix(DNA_ANALYSIS.out.versions)

    // SUBWORKFLOW: detect HLA allele mutations
    DETECT_MUTATIONS (
        GENERATE_DNA_NOVOALIGN_BAMS.out.hla_allele_bams_ch,
        GENERATE_DNA_NOVOALIGN_BAMS.out.passed_heterozygous_hla_alleles_ch,
        GENERATE_REFERENCES.out.mutation_calling_reference,
        tumour_wxs_sample_count_ch,
        germline_sample_count_ch,
        gl_samples,
        INPUT_CHECK.out.matched_inventory
    )

    ch_versions = ch_versions.mix(DETECT_MUTATIONS.out.versions)


    // SUBWORKFLOW: align rna fqs to personalised reference and generate allele specific bams
               
    GENERATE_RNA_NOVOALIGN_BAMS (
        rna_fqs_ch,
        GENERATE_REFERENCES.out.transcriptome_reference,
        GENERATE_REFERENCES.out.mosdepth_exons_reference
        )
        
    ch_versions = ch_versions.mix(GENERATE_RNA_NOVOALIGN_BAMS.out.versions)

    //
    // SUBWORKFLOW: Detect HLA allele expression and allelic imbalance in all samples
    //              Detect HLA allelic repression in samples with a match normal
    //


    NOVOALIGN_RNA_ANALYSIS (
        rna_library_size_ch,
        GENERATE_REFERENCES.out.transcriptome_snp_positions,
        GENERATE_REFERENCES.out.transcriptome_allele_tables,
        GENERATE_RNA_NOVOALIGN_BAMS.out.hla_allele_bams_ch,
        GENERATE_RNA_NOVOALIGN_BAMS.out.passed_heterozygous_hla_genes_ch,
        GENERATE_RNA_NOVOALIGN_BAMS.out.passed_heterozygous_hla_alleles_ch,
        GENERATE_RNA_NOVOALIGN_BAMS.out.hla_bam_read_count_ch,
        INPUT_CHECK.out.checked_inventory,
        "transcriptome",
        GENERATE_REFERENCES.out.star_input,
        "novoalign"
    )

    ch_versions = ch_versions.mix(NOVOALIGN_RNA_ANALYSIS.out.versions)

    // overwrite rna_all_snps_ch with the output from NOVOALIGN_RNA_ANALYSIS when the rnaseq analysis is run
    rna_all_snps_ch = NOVOALIGN_RNA_ANALYSIS.out.rna_all_snps_output_ch

    GENERATE_RNA_STAR_BAMS (
        GENERATE_REFERENCES.out.mosdepth_reference,
        rna_bams_ch,
        GENERATE_REFERENCES.out.star_input,
        mhc_kmer_ch
    )

    STAR_RNA_ANALYSIS (
        rna_library_size_ch,
        GENERATE_REFERENCES.out.genome_snp_positions,
        GENERATE_REFERENCES.out.genome_allele_tables,
        GENERATE_RNA_STAR_BAMS.out.hla_allele_bams_ch,
        GENERATE_RNA_STAR_BAMS.out.passed_heterozygous_hla_genes_ch,
        GENERATE_RNA_STAR_BAMS.out.passed_heterozygous_hla_alleles_ch,
        GENERATE_RNA_STAR_BAMS.out.hla_bam_read_count_ch,
        INPUT_CHECK.out.checked_inventory,
        "genome",
        GENERATE_REFERENCES.out.star_input,
        "star"
    )

    // updates the mosdepth cohort channel only if RNA samples are being run
    cohort_mosdepth_input_ch = cohort_mosdepth_input_ch
        .mix(GENERATE_RNA_STAR_BAMS.out.mosdepth_ch)
        .mix(GENERATE_RNA_NOVOALIGN_BAMS.out.mosdepth_ch)
    // overwrite rna_bam_read_count_ch with the output from GENERATE_RNA_NOVOALIGN_BAMS when the rnaseq analysis is run
    rna_bam_read_count_ch = GENERATE_RNA_NOVOALIGN_BAMS.out.hla_bam_read_count_ch


        // SUBWORKFLOW: Generate STAR aligned bams to detect exon skipping and intron retention

        ALT_SPLICING ( 
            GENERATE_RNA_STAR_BAMS.out.splice_junction_ch,
            GENERATE_REFERENCES.out.star_input,
            codon_table_ch
        )

        inventory_table_format = INPUT_CHECK.out.matched_inventory
            .map { meta, inventory -> tuple(meta.patient_id, inventory) }


        // Collate known splicing cohort table
        cohort_known_splicing_input_ch = ALT_SPLICING.out.known_sjs_ch
        cohort_known_splicing_input_ch
            .map { [it[0], it[1].getName()] }
            .groupTuple()
            .set{ cohort_known_splicing_input_ch_file_collector }


        cohort_known_file_collector( cohort_known_splicing_input_ch_file_collector, "CREATE_KNOWN_SPLICING_COHORT_TABLE" )
        cohort_known_file_collector.out.csvs.set{ cohort_known_splicing_input_file }

        //cohort_known_splicing_input_file.view()
        
        CREATE_KNOWN_SPLICING_COHORT_TABLE(
            cohort_known_splicing_input_ch.join(INPUT_CHECK.out.matched_inventory).join(cohort_known_splicing_input_file)
        )    

        // Collate novel splicing cohort table
        cohort_novel_splicing_input_ch = ALT_SPLICING.out.novel_sjs_ch
        cohort_novel_splicing_input_ch
            .map { [it[0], it[1].getName()] }
            .groupTuple()
            .set{ cohort_novel_splicing_input_ch_file_collector }


       cohort_novel_file_collector( cohort_novel_splicing_input_ch_file_collector, "CREATE_NOVEL_SPLICING_COHORT_TABLE" )
       cohort_novel_file_collector.out.csvs.set{ cohort_novel_splicing_input_file }

        CREATE_NOVEL_SPLICING_COHORT_TABLE(
            cohort_novel_splicing_input_ch.join(INPUT_CHECK.out.checked_inventory).join(cohort_novel_splicing_input_file)
        )           

        // Collate tumour/normal enrichment splicing cohort table
        // ALT_SPLICING.out.tumour_normal_sjs_ch.first().view()
        cohort_tumour_normal_splicing_input_ch = ALT_SPLICING.out.tumour_normal_sjs_ch
        cohort_tumour_normal_splicing_input_ch
            .map { [it[0], it[1].getName()] }
            .groupTuple()
            .set{ cohort_tumour_normal_splicing_input_ch_file_collector }

        cohort_tumour_normal_file_collector( cohort_tumour_normal_splicing_input_ch_file_collector, "CREATE_TUMOUR_NORMAL_SPLICING_COHORT_TABLE" )
        cohort_tumour_normal_file_collector.out.csvs.set{ cohort_tumour_normal_splicing_input_file }

        CREATE_TUMOUR_NORMAL_SPLICING_COHORT_TABLE(
            cohort_tumour_normal_splicing_input_ch.join(INPUT_CHECK.out.checked_inventory).join(cohort_tumour_normal_splicing_input_file)
        )           

        ch_versions = ch_versions.mix(ALT_SPLICING.out.versions)


    // Collate mosdepth output
   cohort_mosdepth_input_ch = cohort_mosdepth_input_ch
            .mix(GENERATE_DNA_NOVOALIGN_BAMS.out.mosdepth_ch)
            .groupTuple()
            .map{ patient_id, meta, csvs ->
                sorted_csvs = csvs.sort()
                unique_csvs = sorted_csvs.unique { it.toString().tokenize('/').last() }
                tuple(patient_id, unique_csvs)
                }

    cohort_mosdepth_input_ch
            .transpose()
            .map { [it[0], it[1].getName()] }
            .groupTuple()
            .set{ cohort_mosdepth_input_file_collector }


    cohort_mosdepth_file_collector( cohort_mosdepth_input_file_collector, "CREATE_MOSDEPTH_COHORT_TABLE" )
    cohort_mosdepth_file_collector.out.csvs.set{ cohort_mosdepth_input_text_file }

    CREATE_MOSDEPTH_COHORT_TABLE(
        cohort_mosdepth_input_ch.join(inventory_table_format).join(cohort_mosdepth_input_text_file)
    )

    // collate library sizes table - only if flagstat is run
    library_size_input_ch = PREPROCESSING.out.unmapped_count_ch
                .mix(PREPROCESSING.out.mapped_count_ch)
                .groupTuple()
                .map{ patient_id, csvs ->
                    sorted_csvs = csvs.sort()
                    unique_csvs = sorted_csvs.unique { it.toString().tokenize('/').last() }
                    tuple(patient_id, unique_csvs)
                    }

        // Create text file with input csvs for final table
        library_size_input_ch
            .transpose()
            .map { [it[0], it[1].getName()] }
            .groupTuple()
            .set { library_size_input_file_collector }

        library_size_file_collector( library_size_input_file_collector, "CREATE_LIBRARY_SIZE_TABLE" )
        library_size_file_collector.out.csvs.set{ library_size_text_file }


        CREATE_LIBRARY_SIZE_TABLE(
            library_size_input_ch.join(INPUT_CHECK.out.checked_inventory).join(library_size_text_file)
        )


    
    // create cohort mutations table
        

    // collate overview table
    genome_allele_tables_ch = GENERATE_REFERENCES.out.genome_allele_tables
            //.map{patient_id, genome_allele_tables -> tuple(genome_allele_tables)}

    transcriptome_allele_tables_ch = GENERATE_REFERENCES.out.transcriptome_allele_tables
            //.map{patient_id, transcriptome_allele_tables -> tuple(transcriptome_allele_tables)}
    
    dna_library_size_for_output = dna_library_size_ch
            .map{patient_id, sample_name, library_size_path -> tuple(patient_id, library_size_path)}
    rna_library_size_for_output = rna_library_size_ch
            .map{patient_id, sample_name, library_size_path -> tuple(patient_id, library_size_path)}            

    hla_bam_read_count_ch_input = GENERATE_DNA_NOVOALIGN_BAMS.out.hla_bam_read_count_ch
            .map{meta, hla_bam_read_count -> tuple(meta.patient_id, hla_bam_read_count)}
    rna_bam_read_count_ch_input = rna_bam_read_count_ch
            .map{meta, rna_bam_read_count -> tuple(meta.patient_id, rna_bam_read_count)}

    allele_table_input_ch = rna_bam_read_count_ch_input
            .mix(hla_bam_read_count_ch_input)
            .mix(DNA_ANALYSIS.out.dna_all_snps_output_ch)
            .mix(rna_all_snps_ch)
            .mix(genome_allele_tables_ch)
            .mix(transcriptome_allele_tables_ch)
            .mix(rna_library_size_for_output)
            .mix(dna_library_size_for_output)
            .groupTuple()
            .map{ patient_id, csvs ->
                sorted_csvs = csvs.sort()
                unique_csvs = sorted_csvs.unique { it.toString().tokenize('/').last() }
                tuple(patient_id, unique_csvs)
                }
    
    // Create text file with input csvs for final table
    allele_table_input_ch
            .transpose()
            .map { [it[0], it[1].getName()] }
            .groupTuple()
            .set { allele_table_input_file_collector }
    
    allele_table_file_collector( allele_table_input_file_collector, "CREATE_MHC_HAMMER_TABLE_alleles" )
    allele_table_file_collector.out.csvs.set{ allele_input_text_file }

    // We also want a text file that has the germline samples that were used for hlahd
    hlahd_gl_samples_ch = hlahd_input_ch.map{meta, fastqs -> tuple(meta.patient_id, meta.sample_id)}
                .groupTuple()
    
    hlahd_gl_samples_ch
        .set { hlahd_gl_samples_input_file_collector }

    hlahd_gl_samples_file_collector( hlahd_gl_samples_input_file_collector, "CREATE_MHC_HAMMER_TABLE_hlahd" )
    hlahd_gl_samples_file_collector.out.csvs.set{ hlahd_germline_samples_text_file }

    CREATE_MHC_HAMMER_TABLE(
        allele_table_input_ch.join(inventory_table_format).join(allele_input_text_file).join(hlahd_germline_samples_text_file),
        params.min_depth
    )

    

    // CUSTOM_DUMPSOFTWAREVERSIONS (
    //     ch_versions.unique().collectFile(name: 'collated_versions.yml')
    // )


}


/*
========================================================================================
    COMPLETION SUMMARY
========================================================================================
*/


/*
========================================================================================
    THE END
========================================================================================
*/
