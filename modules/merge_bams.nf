process PICARD_MERGESAMFILES {
    tag "$meta.library_name"

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("*.bam"), emit: bam

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.library_name}"
    def bam_files = bams.sort()
    def avail_mem = 8
    if (!task.memory) {
        log.info '[Picard MergeSamFiles] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.giga*0.8).intValue()
    }
    if (bam_files.size() > 1) {
        """
        module load picard

        java -Xmx${avail_mem}G -jar \$PICARD_ROOT/picard.jar MergeSamFiles \\
            $args \\
            ${'--INPUT '+bam_files.join(' --INPUT ')} \\
            --OUTPUT ${prefix}.bam
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            picard: \$( echo \$(picard MergeSamFiles --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
        END_VERSIONS
        """
    } else {
        """
        ln -s ${bam_files[0]} ${prefix}.bam
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            picard: \$( echo \$(picard MergeSamFiles --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
        END_VERSIONS
        """
    }
}