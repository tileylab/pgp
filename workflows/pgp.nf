/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_pgp_pipeline'

// Modules added from nf-core
include { BBMAP_BBDUK } from '../modules/nf-core/bbmap/bbduk/main'
include { BWA_INDEX } from '../modules/nf-core/bwa/index/main'
include { BWA_MEM } from '../modules/nf-core/bwa/mem/main' 
include { FASTP } from '../modules/local/fastp/main'

include { PICARD_MARKDUPLICATES } from '../modules/nf-core/picard/markduplicates/main'

include { GATK4_CREATESEQUENCEDICTIONARY } from '../modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_GENOMICSDBIMPORT } from '../modules/local/gatk4/genomicsdbimport/main'
include { GATK4_GENOTYPEGVCFS } from '../modules/local/gatk4/genotypegvcfs/main'
include { GATK4_HAPLOTYPECALLER } from '../modules/local/gatk4/haplotypecaller/main'
include { GATK4_SELECTVARIANTS } from '../modules/nf-core/gatk4/selectvariants/main'
include { GATK4_SELECTVARIANTS as SELECT_ALL } from '../modules/local/gatk4/selectvariants/main'
include { GATK4_SELECTVARIANTS as SELECT_PASS } from '../modules/local/gatk4/selectvariants/main'
include { GATK4_SELECTVARIANTS as SELECT_BIPASS } from '../modules/local/gatk4/selectvariants/main'
include { GATK4_VARIANTFILTRATION } from '../modules/nf-core/gatk4/variantfiltration/main'

include { SAMTOOLS_FAIDX } from '../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_INDEX } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_PASSDUPS } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_VIEW } from '../modules/nf-core/samtools/view/main'

// VCFtools iterative filtering — single-task processes chained with aliases
include { VCFTOOLS_RECODE as VCF_MINDP        } from '../modules/local/vcftools/recode/main'
include { VCFTOOLS_RECODE as VCF_SITE_1       } from '../modules/local/vcftools/recode/main'
include { VCFTOOLS_RECODE as VCF_SITE_2       } from '../modules/local/vcftools/recode/main'
include { VCFTOOLS_RECODE as VCF_SITE_3       } from '../modules/local/vcftools/recode/main'
include { VCFTOOLS_RECODE as VCF_FINAL        } from '../modules/local/vcftools/recode/main'
include { VCFTOOLS_RECODE as VCF_RMFILTERED   } from '../modules/local/vcftools/recode/main'
include { VCFTOOLS_RECODE as VCF_BIALLELIC    } from '../modules/local/vcftools/recode/main'
include { VCFTOOLS_RECODE as VCF_MAC          } from '../modules/local/vcftools/recode/main'
include { VCFTOOLS_RECODE as VCF_THIN         } from '../modules/local/vcftools/recode/main'
include { VCFTOOLS_RECODE as VCF_INGROUP      } from '../modules/local/vcftools/recode/main'
include { VCFTOOLS_MISSING_INDV as VCF_RMINDV_1 } from '../modules/local/vcftools/missingindv/main'
include { VCFTOOLS_MISSING_INDV as VCF_RMINDV_2 } from '../modules/local/vcftools/missingindv/main'
include { VCFTOOLS_MISSING_INDV as VCF_RMINDV_3 } from '../modules/local/vcftools/missingindv/main'
include { VCFTOOLS_DEPTH_OUTLIER              } from '../modules/local/vcftools/depthoutlier/main'
include { VCFTOOLS_SUMMARY                    } from '../modules/local/vcftools/summary/main'
include { BGZIP_TABIX } from '../modules/local/bgzip_tabix/main'

// Polyploid genotyping (EBG + updog) with an EBG-vs-updog concordance comparison
include { POLYPLOID_GENOTYPING } from '../subworkflows/local/polyploid_genotyping/main'

// Subworkflows added from nf-core
include { BAM_STATS_SAMTOOLS } from '../subworkflows/nf-core/bam_stats_samtools/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PGP {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()
    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

        //
    // Module: Run fastp
    //

    adapter_fasta = []
    discard_trimmed_pass = false
    save_trimmed_fail = false
    save_merged = false

    FASTP (
        ch_samplesheet,
        adapter_fasta,
        discard_trimmed_pass,
        save_trimmed_fail,
        save_merged
    )
    ch_versions = ch_versions.mix(FASTP.out.versions_fastp.first())
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it[1]})
    // checking that output of FASTP process is a tuple mapping ids to reads
    //FASTP.out.reads.view()
    
    //
    // Module: Run bbmap_bbduk
    //

    fastp_reads = FASTP.out.reads
    contamination_file = channel.fromPath("${projectDir}/assets/phix174_ill.ref.fa", checkIfExists: true)
    BBMAP_BBDUK (
        fastp_reads,
        contamination_file.first()
    )
    ch_versions = ch_versions.mix(BBMAP_BBDUK.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(BBMAP_BBDUK.out.log.collect{it[1]})

    //
    // Module: Run bwa_index
    // Index the reference genome for bwa
    
    //first need to create a channel with meta map connecting reference name with reference genom fasta
    ch_reference = channel.fromPath(params.reference)
        .map { reference -> tuple(reference.simpleName, reference) }
    //ch_reference.view()
   
    BWA_INDEX (
       ch_reference 
    )
    ch_versions = ch_versions.mix(BWA_INDEX.out.versions_bwa.first())

    //
    // Module: Run bwa_mem
    //

    bbduk_reads = BBMAP_BBDUK.out.reads
    indexed_reference = BWA_INDEX.out.index.first()
    //indexed_reference.view()
    sort_bam = true
    BWA_MEM(
        bbduk_reads,
        indexed_reference,
        ch_reference.first(),
        sort_bam
    )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

/*
Report some statistics on the reads mapping
*/
    //
    // Module: Run samtools_index
    // some downstream functions like filtering with samtools_view require the indexed bams
    bwa_bams = BWA_MEM.out.bam
    SAMTOOLS_INDEX (
        bwa_bams
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions_samtools.first())
    //SAMTOOLS_INDEX.out.bai.view()

    //create a tuple with the ids, bams, and indecies for samtools
    ch_bam_bai = SAMTOOLS_INDEX.out.bai.join(bwa_bams).map{meta, bai, bam -> [meta, bam, bai]}
    //ch_bam_bai.view()

    // Create indexed reference genome (needed for SAMTOOLS_VIEW and downstream GATK)
    SAMTOOLS_FAIDX (
        ch_reference,
        false
    )
    reference_faidx = SAMTOOLS_FAIDX.out.fai.first()
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions_samtools.first())
    ch_reference_fasta_fai = ch_reference.join(SAMTOOLS_FAIDX.out.fai).map { meta, fasta, fai -> [ meta, fasta, fai ] }
    ch_reference_fai = ch_reference.join(SAMTOOLS_FAIDX.out.fai).map { meta, fasta, fai -> [ meta, fai ] }

    //
    // Subworkflow: bam_stats_samtools
    // collects the stats, flagstats, and idxstats all in one go

    // consider splitting the existing tuple and accessing the file that way
    //ch_fasta_file = Channel.fromPath(params.reference)
    BAM_STATS_SAMTOOLS (
        ch_bam_bai,
        ch_reference.first()
    )
    ch_multiqc_files = ch_multiqc_files.mix(BAM_STATS_SAMTOOLS.out.stats.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(BAM_STATS_SAMTOOLS.out.flagstat.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(BAM_STATS_SAMTOOLS.out.idxstats.collect{it[1]})

/*
Do some light filtering prior to genotying
GATK should filter out some frass by default such as improperly paired reads and the quality scores will be baked into the genotype quality
But I think some filtering here could allow modularity if we want to brach the process to different genotyping models
Clearing out some junk should ease the computation a little too, which does cost on clod services
*/
    //
    // Module: Run samtools_view
    // Filter bams to include only properly paired mapping q20 reads
    SAMTOOLS_VIEW (
        ch_bam_bai,
        ch_reference_fasta_fai.first(),
        [],
        []
    )
    ch_q20reads = SAMTOOLS_VIEW.out.bam
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions_samtools.first())

/*
 BEGIN: GATK Genotyping
    Create Referece Sequence Dictionary
    Mark Duplicates
    Genotype at individual level with haplotypecaller in gvcf mode
    Create DB
    Joint Genotyping
    Add Annotations
NOTE: Not bothering with splitting on intervals since the microbial genomes are typically very small. Therefore some of the modules have been moved to local and modified.
*/

    // Create Sequence Dictionary from reference genome
    GATK4_CREATESEQUENCEDICTIONARY (
        ch_reference
    )
    reference_dictionary = GATK4_CREATESEQUENCEDICTIONARY.out.dict.first()
    ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions_gatk4.first())
    ch_reference_dict = ch_reference.join(GATK4_CREATESEQUENCEDICTIONARY.out.dict).map { meta, fasta, dict -> [ meta, dict ] }
    //reference_dictionary.view()

    //reference_faidx.view()

    if (params.skip_markduplicates) {
        log.info 'Skipping duplicate marking: --skip_markduplicates is enabled. SAMTOOLS_VIEW-filtered BAMs will be indexed and passed directly to HaplotypeCaller.'
        // Duplicate marking can be inappropriate for targeted/restriction-based libraries.
        SAMTOOLS_INDEX_PASSDUPS (
            ch_q20reads
        )
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX_PASSDUPS.out.versions_samtools.first())
        ch_bam_bai_marked = SAMTOOLS_INDEX_PASSDUPS.out.bai.join(ch_q20reads).map{meta, bai, bam -> [meta, bam, bai]}
    } else {
        log.info 'Duplicate marking is enabled. Running PICARD_MARKDUPLICATES before HaplotypeCaller.'
        // Mark duplicates for standard random-sheared library prep.
        PICARD_MARKDUPLICATES (
            ch_q20reads,
            ch_reference.first(),
            reference_faidx.first()
        )
        ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions_picard.first())
        ch_bam_bai_marked = PICARD_MARKDUPLICATES.out.bam.join(PICARD_MARKDUPLICATES.out.bai).map{meta, bam, bai -> [meta, bam, bai]}
    }
    //ch_bam_bai_marked.view()

    // Genotype at individual level with gatk
    GATK4_HAPLOTYPECALLER (
        ch_bam_bai_marked,
        ch_reference.first(),
        reference_faidx,
        reference_dictionary
    )
    ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions_gatk4.first())
    //GATK4_HAPLOTYPECALLER.out.vcf.view()
    
    
    // created channel of lists for genomicsbdimport based on sarek here https://github.com/nf-core/sarek/blob/5fe5cdff171e3baed603b3990cab7f7fd3fcb992/subworkflows/local/bam_variant_calling_haplotypecaller/main.nf#L11
    
    //first need to create a channel with meta map connecting reference name with reference genom fasta
    ch_intervals = channel.fromPath(params.intervals)
    //ch_intervals.view()
    
    gvcfs = GATK4_HAPLOTYPECALLER.out.vcf.join(GATK4_HAPLOTYPECALLER.out.tbi, failOnMismatch: true).map{
        meta, gvcf, tbi -> [ meta, gvcf, tbi, 0]
    }
    //gvcfs.view()
    gendb_input = gvcfs.map{
       meta, gvcf, tbi, interval -> [
        [id:'joint_genotyping'], gvcf, tbi, interval
       ] 
    }.groupTuple(by:3)
    .map{
        meta_list, gvcf, tbi, interval -> [
            meta_list[0], gvcf, tbi, interval
        ]
    }
    //gendb_input.view()

    // Create Genomics DataBase for GATK4 joint genotyping
    
    GATK4_GENOMICSDBIMPORT (
        gendb_input,
        ch_intervals
    )
    ch_versions = ch_versions.mix(GATK4_GENOMICSDBIMPORT.out.versions_gatk4.first())
    genotype_input = GATK4_GENOMICSDBIMPORT.out.genomicsdb.map{ meta, genomicsdb -> [ meta, genomicsdb, [], [], [] ] }
    //genotype_input.view()

    GATK4_GENOTYPEGVCFS (
       genotype_input,
       ch_reference.first(),
       reference_faidx,
       reference_dictionary
    )
    ch_versions = ch_versions.mix(GATK4_GENOTYPEGVCFS.out.versions_gatk4.first())
    multisample_vcf_tbi = GATK4_GENOTYPEGVCFS.out.vcf.join(GATK4_GENOTYPEGVCFS.out.tbi).map{
        meta, vcf, tbi -> [meta, vcf, tbi]
    }
    //multisample_vcf_tbi.view()

    /*
        BEGIN: Filtering of the multisample VCF
        At this point the joint vcf is ready but has a lot of false positives and includes structural variants. It might be desirable to take this and start from the unfiltered results. These will be avilable as an output but I still try to address the best-practices-like filtering based on GATK annotations followed by some reasonable hard filtering with VCFtools.
        The nf-core GATK4_SELECTVARIANTS module was modified to get rid of an unnecessary intervals option. I trust that the desired intervals have been handled upfront.
    */

    SELECT_ALL (
        multisample_vcf_tbi
    )
    ch_versions = ch_versions.mix(SELECT_ALL.out.versions_gatk4.first())
    selectall_vcf_tbi = SELECT_ALL.out.vcf.join(SELECT_ALL.out.tbi).map{
        meta, vcf, tbi -> [meta, vcf, tbi]
    }
    //selectall_vcf_tbi.view()

    GATK4_VARIANTFILTRATION (
        selectall_vcf_tbi,
        ch_reference,
        ch_reference_fai,
        ch_reference_dict
    )
    ch_versions = ch_versions.mix(GATK4_VARIANTFILTRATION.out.versions_gatk4.first())
    annotated_vcf_tbi = GATK4_VARIANTFILTRATION.out.vcf.join(GATK4_VARIANTFILTRATION.out.tbi).map{
        meta, vcf, tbi -> [meta, vcf, tbi]
    }
    //annotated_vcf_tbi.view()

    SELECT_PASS (
        annotated_vcf_tbi
    )
    ch_versions = ch_versions.mix(SELECT_PASS.out.versions_gatk4.first())
    pass_vcf_tbi = SELECT_PASS.out.vcf.join(SELECT_PASS.out.tbi).map{
        meta, vcf, tbi -> [meta, vcf, tbi]
    }
    //pass_vcf_tbi.view()

    SELECT_BIPASS (
        pass_vcf_tbi
    )
    ch_versions = ch_versions.mix(SELECT_BIPASS.out.versions_gatk4.first())
    bipass_vcf_tbi = SELECT_BIPASS.out.vcf.join(SELECT_BIPASS.out.tbi).map{
        meta, vcf, tbi -> [meta, vcf, tbi]
    }
    //bipass_vcf_tbi.view()

    /*
        Iterative VCFtools filtering of the joint VCF, decomposed into single-task
        processes chained with aliases (mirrors the SELECT_* pattern). Operates on
        the .snps.gatkfilters joint VCF (SELECT_PASS). Thresholds are tuned via the
        vcf_* params; outgroups are flagged in the samplesheet (outgroup column).
        Each stage's VCF is collected for the per-stage summary.
    */
    // Outgroup sample IDs (empty channel when no sample is flagged -> ingroup step is skipped)
    ch_outgroups = ch_samplesheet
        .filter { meta, reads -> (meta.outgroup as Integer) == 1 }
        .map    { meta, reads -> meta.id }
        .collectFile(name: 'outgroups.txt', newLine: true)

    ch_stages = channel.empty()

    // Seed: mean-depth floor
    VCF_MINDP ( pass_vcf_tbi.map { meta, vcf, tbi -> [ meta, vcf ] }, [] )
    ch_stages = ch_stages.mix(VCF_MINDP.out.vcf.map { meta, vcf -> [ meta, '01_min_meanDP', vcf ] })

    // Iterative site/individual missingness rounds
    VCF_SITE_1   ( VCF_MINDP.out.vcf, [] )
    ch_stages = ch_stages.mix(VCF_SITE_1.out.vcf.map { meta, vcf -> [ meta, '02_site_1', vcf ] })
    VCF_RMINDV_1 ( VCF_SITE_1.out.vcf )
    ch_stages = ch_stages.mix(VCF_RMINDV_1.out.vcf.map { meta, vcf -> [ meta, '03_indiv_1', vcf ] })

    VCF_SITE_2   ( VCF_RMINDV_1.out.vcf, [] )
    ch_stages = ch_stages.mix(VCF_SITE_2.out.vcf.map { meta, vcf -> [ meta, '04_site_2', vcf ] })
    VCF_RMINDV_2 ( VCF_SITE_2.out.vcf )
    ch_stages = ch_stages.mix(VCF_RMINDV_2.out.vcf.map { meta, vcf -> [ meta, '05_indiv_2', vcf ] })

    VCF_SITE_3   ( VCF_RMINDV_2.out.vcf, [] )
    ch_stages = ch_stages.mix(VCF_SITE_3.out.vcf.map { meta, vcf -> [ meta, '06_site_3', vcf ] })
    VCF_RMINDV_3 ( VCF_SITE_3.out.vcf )
    ch_stages = ch_stages.mix(VCF_RMINDV_3.out.vcf.map { meta, vcf -> [ meta, '07_indiv_3', vcf ] })

    // Depth-outlier cut
    VCFTOOLS_DEPTH_OUTLIER ( VCF_RMINDV_3.out.vcf )
    ch_stages = ch_stages.mix(VCFTOOLS_DEPTH_OUTLIER.out.vcf.map { meta, vcf -> [ meta, '08_max_meanDP', vcf ] })

    // Optional final site-missingness pass
    ch_pre_clean = VCFTOOLS_DEPTH_OUTLIER.out.vcf
    if (params.vcf_max_missing_final != null) {
        VCF_FINAL ( VCFTOOLS_DEPTH_OUTLIER.out.vcf, [] )
        ch_stages = ch_stages.mix(VCF_FINAL.out.vcf.map { meta, vcf -> [ meta, '09_final_max_missing', vcf ] })
        ch_pre_clean = VCF_FINAL.out.vcf
    }

    // Cleaned VCF + analysis-ready subsets
    VCF_RMFILTERED ( ch_pre_clean, [] )
    ch_stages = ch_stages.mix(VCF_RMFILTERED.out.vcf.map { meta, vcf -> [ meta, '10_filtered', vcf ] })

    VCF_BIALLELIC ( VCF_RMFILTERED.out.vcf, [] )
    ch_stages = ch_stages.mix(VCF_BIALLELIC.out.vcf.map { meta, vcf -> [ meta, '11_biallelic', vcf ] })
    VCF_MAC       ( VCF_BIALLELIC.out.vcf, [] )
    ch_stages = ch_stages.mix(VCF_MAC.out.vcf.map { meta, vcf -> [ meta, '12_biallelic_mac', vcf ] })
    VCF_THIN      ( VCF_MAC.out.vcf, [] )
    ch_stages = ch_stages.mix(VCF_THIN.out.vcf.map { meta, vcf -> [ meta, '13_biallelic_mac_thin', vcf ] })

    // Ingroup-only product: channel-driven, runs only when >= 1 outgroup is flagged
    ch_ingroup_in = VCF_RMFILTERED.out.vcf.combine(ch_outgroups)
    VCF_INGROUP (
        ch_ingroup_in.map { meta, vcf, og -> [ meta, vcf ] },
        ch_ingroup_in.map { meta, vcf, og -> og }
    )
    ch_stages = ch_stages.mix(VCF_INGROUP.out.vcf.map { meta, vcf -> [ meta, '14_ingroup', vcf ] })

    // Per-stage SNP/individual summary + MultiQC content (00_input = the raw joint VCF)
    ch_stages = ch_stages.mix(pass_vcf_tbi.map { meta, vcf, tbi -> [ meta, '00_input', vcf ] })
    ch_summary = ch_stages
        .toSortedList { a, b -> a[1] <=> b[1] }
        .map { rows -> tuple(rows[0][0], rows.collect { it[1] }, rows.collect { it[2] }) }
    VCFTOOLS_SUMMARY ( ch_summary )
    ch_multiqc_files = ch_multiqc_files.mix(VCFTOOLS_SUMMARY.out.mqc.flatten())

    // Compress + index the final products
    ch_products = VCF_RMFILTERED.out.vcf
        .mix(VCF_BIALLELIC.out.vcf, VCF_MAC.out.vcf, VCF_THIN.out.vcf, VCF_INGROUP.out.vcf)
        .groupTuple()
    BGZIP_TABIX ( ch_products )

    /*
        Polyploid genotyping (EBG + updog) on the candidate biallelic SNP set
        (ingroup, pre-MAC/thin). Samples are grouped by samplesheet ploidy.
        Disabled with --no_polyploids (runs only the standard diploid front end).
    */
    if (!params.no_polyploids) {
        // Candidate VCF: ingroup-biallelic when outgroups are flagged, else biallelic
        ch_candidate = VCF_INGROUP.out.vcf.concat(VCF_BIALLELIC.out.vcf).first()

        // Ingroup samples (outgroup == 0): ploidy map, distinct ploidies, BAMs
        ch_ingroup_sheet = ch_samplesheet.filter { meta, reads -> (meta.outgroup as Integer) == 0 }
        ch_ploidy_map = ch_ingroup_sheet
            .map { meta, reads -> "${meta.id}\t${meta.ploidy}" }
            .collectFile(name: 'sample_ploidy.tsv', newLine: true)
        ch_ploidies = ch_ingroup_sheet.map { meta, reads -> (meta.ploidy as Integer) }.unique()

        ch_ingroup_bb   = ch_bam_bai_marked.filter { meta, bam, bai -> (meta.outgroup as Integer) == 0 }
        ch_ingroup_bams = ch_ingroup_bb.map { meta, bam, bai -> bam }.collect()
        ch_ingroup_bais = ch_ingroup_bb.map { meta, bam, bai -> bai }.collect()

        POLYPLOID_GENOTYPING (
            ch_candidate,
            ch_ingroup_bams,
            ch_ingroup_bais,
            ch_reference.first(),
            reference_faidx.map { meta, fai -> fai },
            ch_ploidy_map,
            ch_ploidies
        )
    }

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_versions_files = ch_versions.filter { it instanceof Path }
    softwareVersionsToYAML(ch_versions_files.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'pgp_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
