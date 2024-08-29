nextflow.enable.dsl=2

process PICARD_MARKDUPLICATES {
    tag "$tumorName"
    publishDir  "${params.test_data}/marDup/output", mode: "copy"

    input:
    path reads
    val picard_module
    val tumorName

    output:
    tuple val(tumorName), path("*.bam") , emit: bam,  optional: true
    tuple val(tumorName), path("*.bai") , emit: bai,  optional: true
    tuple val(tumorName), path("*.cram"), emit: cram, optional: true
    tuple val(tumorName), path("*.metrics.txt"), emit: metrics
    path  "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${tumorName}"
    def suffix = task.ext.suffix    ?: "${reads.getExtension()}"
    def avail_mem = 30
    
    if ("$reads" == "${prefix}.${suffix}") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"

    """
    module load $picard_module

    java -Xmx${avail_mem}G -jar /.mounts/labs/gsi/modulator/sw/Ubuntu20.04/picard-2.19.2/picard.jar MarkDuplicates \
        TMP_DIR=picardTmp \\
        ASSUME_SORTED=true \\
        VALIDATION_STRINGENCY=LENIENT \\
        $args \\
        INPUT=$reads \\
        OUTPUT=${prefix}.${suffix} \\
        METRICS_FILE=${prefix}.MarkDuplicates.metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard MarkDuplicates --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${tumorName}"
    def suffix = task.ext.suffix    ?: "${reads.getExtension()}"
    if ("$reads" == "${prefix}.${suffix}") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${prefix}.${suffix}
    touch ${prefix}.MarkDuplicates.metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard MarkDuplicates --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}