nextflow.enable.dsl=2

include {BCL2FASTQ} from './workflows/bcl2fastq'
include {bwaMem} from './modules/bwamem'
include {vep} from "./workflows/vep"
include {mutect2} from "./modules/mutect2"
include {delly} from './workflows/delly'

workflow {
    
    channel
        .of([meta: [id: params.bcl2fastq.meta_id, lane: params.bcl2fastq.meta_lane], samplesheet: file(params.bcl2fastq.samplesheet), run_dir: file(params.bcl2fastq.run_dir)])
        .set { bcl2fastq_input }
    
    BCL2FASTQ(bcl2fastq_input)

    fastq_input = channel
        .fromPath("${params.test_data}/bwamem/input/*${params.bwamem.meta}*.fastq.gz")
        .map { file -> 
            return tuple(params.bwamem.meta, file)
        }
        .groupTuple()

    bwaMem_reads = fastq_input
        .map { key, files ->
            def r1 = files.find { it.name.endsWith('_R1.fastq.gz') }
            def r2 = files.find { it.name.endsWith('_R2.fastq.gz') }
            tuple(params.bwamem.meta, r1, r2)
        }

    bwaMem(
        bwaMem_reads,
        channel.of(params.bwamem.readGroups),
        channel.of(params.bwamem.reference),
        channel.of(params.bwamem.sort_bam),
        channel.of(params.bwamem.threads),
        channel.of(params.bwamem.addParem)
    )
   
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
        mutect2_bams = tumor_bam_files.combine(normal_bam_files)
        mutect2_indexes = tumor_bam_index_files.combine(normal_bam_index_files)
    }
    mutect2_bams.view(it -> "mutect2_bams: $it")
    mutect2_indexes.view(it -> "mutect2_indexes: $it")

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

    tumor_bam_files = channel.fromPath("${params.test_data}/delly/input/*${params.delly.tumorName}*.bam")
    tumor_bam_index_files = channel.fromPath("${params.test_data}/delly/input/*${params.delly.tumorName}*.bam.bai")

    if (!params.delly.tumor_only_mode) {
        normal_bam_files = channel.fromPath("${params.test_data}/delly/input/*${params.delly.normalName}*.bam")
        normal_bam_index_files = channel.fromPath("${params.test_data}/delly/input/*${params.delly.normalName}*.bam.bai")
    }

    if (params.delly.tumor_only_mode) {
        delly_bams = tumor_bam_files
        delly_indexes = tumor_bam_index_files
    } else {
        delly_bams = tumor_bam_files.combine(normal_bam_files)
        delly_indexes = tumor_bam_index_files.combine(normal_bam_index_files)
    }
    
    delly (
        delly_bams,
        delly_indexes,
        params.delly.tumorName,
        params.delly.markdup,
        params.delly.reference,
        params.delly.picard_module
    )
}