nextflow.enable.dsl=2

include {BCL2FASTQ} from './workflows/bcl2fastq'
include {bwaMem} from './workflows/bwamem'
include {vep} from "./workflows/vep"
include {mutect2} from "./workflows/mutect2"
include {delly} from './workflows/delly'

workflow {
/*
    channel
        .of([meta: [id: params.bcl2fastq.meta_id, lane: params.bcl2fastq.meta_lane], samplesheet: file(params.bcl2fastq.samplesheet), run_dir: file(params.bcl2fastq.run_dir)])
        .set { bcl2fastq_input }
    
    BCL2FASTQ(bcl2fastq_input)

    fastq_input = channel
    .fromPath("${params.test_data}/inputs/*${params.bwamem.meta}*.fastq.gz")
    .map { file -> tuple(params.bwamem.meta, file) }
    .groupTuple()

    bwaMem_reads = fastq_input
        .flatMap { meta, files ->
        // Assume fastq file has format R{1,2}.fastq.gz
            def pairs = files.groupBy { it.name.tokenize('_').dropRight(1).join('_') }
            pairs.collect { prefix, pairFiles ->
                def r1 = pairFiles.find { it.name.endsWith('_R1.fastq.gz ') }
                def r2 = pairFiles.find { it.name.endsWith('_R2.fastq.gz ') }
                if (r1 && r2) {
                    tuple(meta, r1, r2)
                } else {
                    null
                }
            }.findAll { it != null }
        }
*/
// Read meta.txt and create a channel
fastq_inputs = Channel
    .fromPath("${params.test_data}/meta.txt")
    .splitCsv(sep: ' ')
    .map { row -> tuple(row[0], file("${params.test_data}/inputs/${row[1]}")) }
    .groupTuple()
    .map { meta, files -> 
        def r1 = files.find { it.name.endsWith("_R1.fastq.gz") }
        def r2 = files.find { it.name.endsWith("_R2.fastq.gz") }
        if (r1 && r2) {
            return tuple(meta, r1, r2)
        } else {
            return null
        }
    }
    .filter { it != null }

    bwaMem(
        fastq_inputs,
        channel.of(params.bwamem.readGroups),
        channel.of(params.bwamem.reference),
        channel.of(params.bwamem.sort_bam),
        channel.of(params.bwamem.threads),
        channel.of(params.bwamem.addParem)
    )

 
    def tumor_bam_files = bwaMem.out.bam
    .filter { it.toString().contains(params.mutect2.tumor_meta) }

    def tumor_bam_index_files = bwaMem.out.bai
    .filter { it.toString().contains(params.mutect2.tumor_meta) }
    //tumor_bam_files.view(it -> "tumor_bam: ${it}")

    if (params.mutect2.tumor_only_mode) {
        mutect2_bams = tumor_bam_files
        mutect2_indexes = tumor_bam_index_files
    } else {
        def normal_bam_files = bwaMem.out.bam
        .filter { it.toString().contains(params.mutect2.normal_meta) }

        def normal_bam_index_files = bwaMem.out.bai
            .filter { it.toString().contains(params.mutect2.normal_meta) }
        mutect2_bams = tumor_bam_files.combine(normal_bam_files)
        mutect2_indexes = tumor_bam_index_files.combine(normal_bam_index_files)
    }

    mutect2_bams.view(it -> "mutect2_bams: $it")


/*
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
   

/*   
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
    */
}