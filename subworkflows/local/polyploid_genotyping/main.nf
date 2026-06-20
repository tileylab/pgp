//
// Polyploid genotyping of the candidate biallelic SNP set with EBG and updog,
// followed by an EBG-vs-updog concordance comparison.
//
// Samples are grouped by their samplesheet ploidy; EBG (diseq) and updog are run
// once per ploidy group. A single per-locus sequencing-error vector (from an
// mpileup over all ingroup BAMs) is shared across groups.
//

include { GENOTYPE_EXTRACT } from '../../../modules/local/genotype/extract/main'
include { SAMTOOLS_MPILEUP } from '../../../modules/local/samtools/mpileup/main'
include { GENOTYPE_ERROR   } from '../../../modules/local/genotype/error/main'
include { EBG_DISEQ        } from '../../../modules/local/ebg/diseq/main'
include { UPDOG            } from '../../../modules/local/updog/main'
include { GENOTYPE_COMPARE } from '../../../modules/local/genotype/compare/main'

workflow POLYPLOID_GENOTYPING {

    take:
    ch_candidate     // channel: [ val(meta), path(vcf) ]   ingroup biallelic VCF (single)
    ch_ingroup_bams  // channel: [ [bam, ...] ]             collected ingroup BAMs
    ch_ingroup_bais  // channel: [ [bai, ...] ]             collected ingroup BAM indexes
    ch_reference     // channel: [ val(name), path(fasta) ]
    ch_fai           // channel: [ path(fai) ]
    ch_ploidy_map    // channel: [ path(sample_ploidy.tsv) ]
    ch_ploidies      // channel: val(ploidy)                distinct ingroup ploidies

    main:

    // 1) Per-ploidy-group read-count matrices (EBG + updog formats)
    ch_extract_in = ch_candidate.combine(ch_ploidies)          // [ meta, vcf, ploidy ]
    GENOTYPE_EXTRACT ( ch_extract_in, ch_ploidy_map )

    // 2) Shared per-locus sequencing error (one mpileup over all ingroup BAMs)
    SAMTOOLS_MPILEUP (
        ch_candidate,
        ch_ingroup_bams,
        ch_ingroup_bais,
        ch_reference,
        ch_fai
    )
    ch_error_in = SAMTOOLS_MPILEUP.out.pileup
        .combine(SAMTOOLS_MPILEUP.out.positions)
        .map { pileup, positions -> [ [id:'genotyping'], pileup, positions ] }
    GENOTYPE_ERROR ( ch_error_in )

    // 3) EBG (per group) — group matrices + the shared error vector
    ch_ebg_in = GENOTYPE_EXTRACT.out.ebg.combine(GENOTYPE_ERROR.out.error)  // [ ploidy, tot, alt, error ]
    EBG_DISEQ ( ch_ebg_in )

    // 4) updog (per group)
    UPDOG ( GENOTYPE_EXTRACT.out.updog )

    // 5) Concordance comparison (once): gather every group's EBG genos, the
    //    sample/loci id files, and the updog calls into one process.
    ch_compare_in = EBG_DISEQ.out.genos.map { ploidy, genos -> genos }
        .mix( GENOTYPE_EXTRACT.out.ids.flatMap { ploidy, samples, loci -> [ samples, loci ] } )
        .mix( UPDOG.out.genos.map { ploidy, genos -> genos } )
        .collect()
        .map { files -> [ [id:'genotyping'], files ] }
    GENOTYPE_COMPARE ( ch_compare_in )

    emit:
    concordance = GENOTYPE_COMPARE.out.concordance
    summary     = GENOTYPE_COMPARE.out.summary
}
