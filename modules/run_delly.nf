nextflow.enable.dsl=2

process RUN_DELLY {
    tag "$sampleName"

    input:
    val dellyMode
    path bams
    path indexes
    val callType
    val sampleName
    val rundelly_module
    path rundelly_fasta
    path rundelly_exclude_list

    output:
    path("${sampleName}.${dellyMode}.${callType}.vcf.gz"), emit: outVcf, optional: true
    path("${sampleName}.${dellyMode}.${callType}.vcf.gz.tbi"), emit: outTbi, optional: true
    path("${sampleName}.${dellyMode}.${callType}_filtered.vcf.gz"), emit: outVcf_filtered, optional: true
    path("${sampleName}.${dellyMode}.${callType}_filtered.vcf.gz.tbi"), emit: outTbi_filtered, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def module_list = rundelly_module.split(',')
    def module_load_cmds = module_list.collect { module -> "module load ${module}" }.join('\n')
    def mappingQuality = 30
    def translocationQuality = 20
    def insertSizeCutoff = 9
    def minClip = 25
    def minCliqueSize = 2
    def minRefSeparation = 25
    def maxReadSeparation = 40
    def jobMemory = 16
    def timeout = 20

    """
    ${module_load_cmds}

    delly call -t $dellyMode \\
      -x $rundelly_exclude_list \\
      -o "${sampleName}.${dellyMode}.${callType}.bcf" \\
      -q $mappingQuality \\
      -s $insertSizeCutoff \\
      -r $translocationQuality \\
      -c $minClip \\
      -z $minCliqueSize \\
      -m $minRefSeparation \\
      -n $maxReadSeparation \\
      -g $rundelly_fasta \\
         $bams

    if [ "$callType" == "somatic" ]; then
        echo "Somatic mode requested, will run delly filtering for somatic SVs"
        bcftools view "${sampleName}.${dellyMode}.${callType}.bcf" | grep ^# | tail -n 1 | \\
        sed 's/.*FORMAT\\t//' | awk -F "\\t" '{print \$1"\\ttumor";print \$2"\\tcontrol"}' > samples.tsv
        delly filter -f somatic -o "${sampleName}.${dellyMode}.${callType}_filtered.bcf" -s samples.tsv \\
                                    "${sampleName}.${dellyMode}.${callType}.bcf"
        bcftools view "${sampleName}.${dellyMode}.${callType}_filtered.bcf" | \\
        bgzip -c > "${sampleName}.${dellyMode}.${callType}_filtered.vcf.gz"
    fi

    bcftools view "${sampleName}.${dellyMode}.${callType}.bcf" | bgzip -c > "${sampleName}.${dellyMode}.${callType}.vcf.gz"

    if [ -e "${sampleName}.${dellyMode}.${callType}.vcf.gz" ]; then
        tabix -p vcf "${sampleName}.${dellyMode}.${callType}.vcf.gz"
    fi
    if [ -e "${sampleName}.${dellyMode}.${callType}_filtered.vcf.gz" ]; then
        tabix -p vcf "${sampleName}.${dellyMode}.${callType}_filtered.vcf.gz"
    fi
    """
}