nextflow.enable.dsl=2
process runDelly {
    tag "$sampleName"

    publishDir  "${params.test_data}/delly/output", mode: 'copy'

    input:
    val dellymode, 
    path bams_and_indexes,
    val callType,
    val sampleName,
    val rundelly_module,
    path rundelly_fasta,
    path rundelly_exclude_list

    output:
    tuple val(meta), path("*.bam")  , emit: bam,    optional: true
    tuple val(meta), path("*.bai")  , emit: bai,    optional: true
    tuple val(meta), path("*.cram") , emit: cram,   optional: true
    tuple val(meta), path("*.csi")  , emit: csi,    optional: true
    tuple val(meta), path("*.crai") , emit: crai,   optional: true
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def module_list = rundelly.split(',')
    def module_load_cmds = module_list.collect { module -> "module load ${module}" }.join('\n')
    Int mappingQuality = 30
    Int translocationQuality = 20
    Int insertSizeCutoff = 9
    Int minClip = 25
    Int minCliqueSize = 2
    Int minRefSeparation = 25
    Int maxReadSeparation = 40
    Int jobMemory = 16
    Int timeout = 20

    """
    ${module_load_cmds}
    delly call -t $dellyMode \\
      -x $rundelly_exclude_list \\
      -o " $sampleName.$dellyMode.$callType.bcf" \\
      -q $mappingQuality \\
      -s $insertSizeCutoff \\
      -r $translocationQuality \\
      -c $minClip \\
      -z $minCliqueSize \\
      -m $minRefSeparation \\
      -n $maxReadSeparation \\
      -g $refFasta \\
         $inputBams

    if [ "$callType" == "somatic" ]; then
    echo "Somatic mode requested, will run delly filtering for somatic SVs"
    bcftools view "$sampleName.$dellyMode.$callType.bcf" | grep ^# | tail -n 1 | \
    sed 's/.*FORMAT\t//' | awk -F "\t" '{print $1"\ttumor";print $2"\tcontrol"' > samples.tsv
    delly filter -f somatic -o "$sampleName.$dellyMode.$callType_filtered.bcf" -s samples.tsv \
                                "$sampleName.$dellyMode.$callType.bcf"
    bcftools view "$sampleName.$dellyMode.$callType_filtered.bcf" | \
    bgzip -c > "$sampleName.$dellyMode.$callType_filtered.vcf.gz"
    fi

    bcftools view "$sampleName.$dellyMode.$callType.bcf" | bgzip -c > "$sampleName.$dellyMode.$callType.vcf.gz"

    if [ -e "$sampleName.$dellyMode.$callType.vcf.gz" ]; then
    tabix -p vcf "$sampleName.$dellyMode.$callType.vcf.gz"
    fi
    if [ -e "$sampleName.$dellyMode.$callType_filtered.vcf.gz" ]; then
    tabix -p vcf "$sampleName.$dellyMode.$callType_filtered.vcf.gz"
    fi
"""

}