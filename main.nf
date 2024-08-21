nextflow.enable.dsl=2

include {vep} from "./workflows/vep"
include {mutect2} from "./modules/mutect2"


workflow {
    tumor_bam_files = channel.fromPath("${params.test_data}/mutect2/input/*${params.mutect2.tumor_meta}*.bam")
    tumor_bam_index_files = channel.fromPath("${params.test_data}/mutect2/input/*${params.mutect2.tumor_meta}*.bai")

    if (!params.mutect2.tumor_only_mode) {
        normal_bam_files = channel.fromPath("${params.test_data}/mutect2/input/*${params.mutect2.normal_meta}*.bam")
        normal_bam_index_files = channel.fromPath("${params.test_data}/mutect2/input/*${params.mutect2.normal_meta}*.bai")
    }

    if (params.mutect2.tumor_only_mode) {
        mutect2_bams = tumor_bam_files
        mutect2_indexes = tumor_bam_index_files
    } else {
        mutect2_bams = tumor_bam_files.mix(normal_bam_files)
        mutect2_indexes = tumor_bam_index_files.mix(normal_bam_index_files)
    }

    mutect2(
        channel.value(params.mutect2.tumor_meta),
        mutect2_bams,
        mutect2_indexes,
        channel.fromPath(params.mutect2.intervals),
        params.mutect2.panel_of_normals ? channel.fromPath(params.mutect2.panel_of_normals) : channel.value('NO_PON'),
        params.mutect2.panel_of_normals_tbi ? channel.fromPath(params.mutect2.panel_of_normals_tbi) :  channel.value('NO_PON_TBI'),
        channel.value(params.mutect2.reference),
        channel.value(params.mutect2.gatk)
    )

    tumor_name = params.vep.tumorName
    reference = params.vep.reference
    normal_name = params.vep.normalName
    vep_tumor_only = params.vep.onlyTumor
    target_bed = params.vep.targetBed

    vep_vcf = channel.fromPath("${params.test_data}/vep/input/*${tumor_name}*.vcf.gz")

    vep(
        channel.value(tumor_name),
        vep_vcf,
        channel.value(reference),
        channel.value(normal_name),
        channel.value(vep_tumor_only),
        channel.value(target_bed)
    )
}