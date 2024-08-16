nextflow.enable.dsl=2

include {BCL2FASTQ} from './bcl2fastq/bcl2fastq'
include {BWA_MEM} from './bwaMem/bwaMem'
include {GATK4_MUTECT2} from "./mutect2/mutect2.nf"
include {ENSEMBLVEP_VEP} from "./vep/vep.nf"

workflow {
    Channel
        .of([meta: [id: params.meta_id, lane: params.meta_lane], samplesheet: file(params.samplesheet), run_dir: file(params.run_dir)])
        .set { bcl2fastq_input }
    BCL2FASTQ(bcl2fastq_input)

    bcl2fastq_out = BCL2FASTQ.out.fastq
    .filter { meta, files ->
        files.any { it.name.contains('GNLR_0026_Bn_U_PE_688_WG_1_S17') }
    }

    bwaMem_reads = bcl2fastq_out
        .map { meta, files ->
            def r1 = files.find { it.name.endsWith('_R1_001.fastq.gz') }
            def r2 = files.find { it.name.endsWith('_R2_001.fastq.gz') }
            tuple(meta, r1, r2)
        }
    BWA_MEM(
        bwaMem_reads,
        channel.of(params.ref_fasta),
        channel.of(params.sort_bam),
        channel.of(params.threads)
    )
   /* 
    params.normal_bam = ''
    params.normal_bam_index = ''
    bam_bai_inputs = tuple(
        params.tumor_meta,
        params.normal_bam ? [file(params.tumor_bam), file(params.normal_bam)] : [file(params.tumor_bam)],
        params.normal_bam_index ? [file(params.tumor_bam_index), file(params.normal_bam_index)] : [file(params.tumor_bam_index)]
    )
    */
    Channel
        .fromPath("${params.mutect2_inputdir}/*${params.tumor_meta}*.bam")
        .set { tumor_bam_files }

    Channel
        .fromPath("${params.mutect2_inputdir}/*${params.tumor_meta}*.bai")
        .set { tumor_bam_index_files }

    if (params.tumor_only_mode) {
        normal_bam_files = Channel.empty()
        normal_bam_index_files = Channel.empty()
    } else {
        normal_bam_files = Channel.fromPath("${params.mutect2_inputdir}/*${params.normal_meta}*.bam")
        normal_bam_index_files = Channel.fromPath("${params.mutect2_inputdir}/*${params.normal_meta}*.bai")
    }

    if (params.tumor_only_mode) {
        paired_bams_and_indexes = tumor_bam_files
            .combine(tumor_bam_index_files)
            .map { tumor_bam, tumor_bam_index ->
                return [
                    params.tumor_meta, 
                    [tumor_bam], 
                    [tumor_bam_index]
                ]
            }
    } else {
        paired_bams_and_indexes = tumor_bam_files
            .combine(normal_bam_files)
            .combine(tumor_bam_index_files)
            .combine(normal_bam_index_files)
            .map { tumor_bam, normal_bam, tumor_bam_index, normal_bam_index ->
                return [
                    params.tumor_meta, 
                    [tumor_bam, normal_bam], 
                    [tumor_bam_index, normal_bam_index]
                ]
            }
    }

    GATK4_MUTECT2(
        paired_bams_and_indexes,
        channel.fromPath(params.intervals),
        channel.fromPath(params.ref_fasta),
        channel.fromPath(params.ref_fai),
        channel.fromPath(params.ref_dict),
        params.germline_resource ? channel.fromPath(params.germline_resource) : channel.of([]),
        params.germline_resource_tbi ? channel.fromPath(params.germline_resource_tbi) : channel.of([]),
        params.panel_of_normals ? channel.fromPath(params.panel_of_normals) : channel.of([]),
        params.panel_of_normals_tbi ? channel.fromPath(params.panel_of_normals_tbi) : channel.of([])
    )

    Channel
        .fromPath("${params.vep_inputdir}/*${params.tumorName}*.vcf.gz")
        .map { file ->
            tuple(params.tumorName, file)
        }
        .set { vep_inputs }

    ENSEMBLVEP_VEP(
        vep_inputs,
        "hg38"
    )
}


