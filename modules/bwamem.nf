nextflow.enable.dsl=2
include {BWA_MEM} from './run_bwamem'

workflow bwaMem {
    take:
    bwamem_reads
    readGroups
    reference
    sort_bam
    threads
    addParem

    main:

    def bwaMem_modules_by_genome = [
        hg19: "samtools/1.9 bwa/0.7.17 hg19-bwa-index/0.7.17",
        hg38: "samtools/1.9 bwa/0.7.17 hg38-bwa-index-with-alt/0.7.17",
        mm10: "samtools/1.9 bwa/0.7.17 mm10-bwa-index/0.7.17"
    ]
    def bwaMemRef_by_genome = [
        hg19: "/.mounts/labs/gsi/modulator/sw/data/hg19-bwa-index-0.7.17/hg19_random.fa",
        hg38: "/.mounts/labs/gsi/modulator/sw/data/hg38-bwa-index-with-alt-0.7.17/hg38_random.fa",
        mm10: "/.mounts/labs/gsi/modulator/sw/data/mm10-bwa-index-0.7.17/mm10.fa"
    ]

    bwaMemRef = reference.map { ref -> bwaMemRef_by_genome[ref]}
    modules = reference.map { ref -> bwaMem_modules_by_genome[ref]}

    BWA_MEM (
        bwamem_reads,
        readGroups,
        bwaMemRef,
        sort_bam,
        threads,
        modules,
        addParem
    )

    emit:
    bam = BWA_MEM.out.bam
    bai = BWA_MEM.out.bai
}