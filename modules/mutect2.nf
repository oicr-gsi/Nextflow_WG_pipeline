include {GATK4_MUTECT2} from "./gatk4_mutect2"

workflow mutect2 {

    take:

    tumor_meta
    input_bam
    input_index
    intervalFile
    pon
    ponIdx
    reference
    gatk

    main:

    def GenomeResources = [
        hg19: [
                refDict : "/.mounts/labs/gsi/modulator/sw/data/hg19-p13/hg19_random.dict",
                refFai : "/.mounts/labs/gsi/modulator/sw/data/hg19-p13/hg19_random.fa.fai",
                refFasta : "/.mounts/labs/gsi/modulator/sw/data/hg19-p13/hg19_random.fa",
                modules : "hg19/p13 samtools/1.9",
                gnomad: "no_gnomad",
                gnomadIdx: "no_gnomadIdx"
        ],
        hg38: [
                refDict : "/.mounts/labs/gsi/modulator/sw/data/hg38-p12/hg38_random.dict",
                refFai : "/.mounts/labs/gsi/modulator/sw/data/hg38-p12/hg38_random.fa.fai",
                refFasta : "/.mounts/labs/gsi/modulator/sw/data/hg38-p12/hg38_random.fa",
                gnomad: "/.mounts/labs/gsi/modulator/sw/data/hg38-gatk-gnomad-2.0/af-only-gnomad.hg38.vcf.gz",
                gnomadIdx: "/.mounts/labs/gsi/modulator/sw/data/hg38-gatk-gnomad-2.0/af-only-gnomad.hg38.vcf.gz.tbi",
                modules : "hg38/p12 samtools/1.9 hg38-gatk-gnomad/2.0"
        ],
        mm10: [
                refDict : "/.mounts/labs/gsi/modulator/sw/data/mm39-p6/mm10.dict",
                refFai : "/.mounts/labs/gsi/modulator/sw/data/mm39-p6/mm10.fa.fai",
                refFasta : "/.mounts/labs/gsi/modulator/sw/data/mm39-p6/mm10.fa",
                modules : "mm10/p6 samtools/1.9",
                gnomad: "no_gnomad",
                gnomadIdx: "no_gnomadIdx"
            ]
    ]
    reference
    .map { ref ->
        return [
            refDict: GenomeResources[ref]['refDict'],
            refFasta: GenomeResources[ref]['refFasta'],
            refFai: GenomeResources[ref]['refFai'],
            gnomad: GenomeResources[ref]['gnomad'],
            gnomadIdx: GenomeResources[ref]['gnomadIdx'],
            modules: GenomeResources[ref]['modules']
        ]
    }
    .set { mutect2_params }

    GATK4_MUTECT2(
        tumor_meta,
        input_bam,
        input_index,
        intervalFile,
        mutect2_params.refFasta,
        mutect2_params.refFai,
        mutect2_params.refDict,
        mutect2_params.gnomad,
        mutect2_params.gnomadIdx,
        pon,
        ponIdx,
        gatk
    )

    emit: 
    vcf = GATK4_MUTECT2.out.vcf
    tbi = GATK4_MUTECT2.out.tbi
    stats = GATK4_MUTECT2.out.stats
}

