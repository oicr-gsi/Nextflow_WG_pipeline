nextflow.enable.dsl = 2
process MERGE_AND_ZIP {
    tag "$sampleName"
    publishDir "${params.test_data}/delly/output", mode: 'copy'

    input:
    path inputVcfs
    path inputTbis
    val sampleName
    val callType
    val modules 
    val prefix
    val variantSupport

    output:
    path("${sampleName}.${callType}${prefix}.vcf.gz"), emit: dellyMergedVcf
    path("${sampleName}.${callType}${prefix}.vcf.gz.tbi"), emit: dellyMergedTabixIndex
    path("${sampleName}.${callType}${prefix}.pass.vcf.gz"), emit: dellyMergedPassVcf, optional: true
    path("${sampleName}.${callType}${prefix}.pass.vcf.gz.tbi"), emit: dellyMergedPassTabixIndex, optional: true

    script:
    """
    set -eu -o pipefail

    module load ${modules}

    vcf-concat ${inputVcfs.join(' ')} | vcf-sort | bgzip -c > "${sampleName}.${callType}${prefix}.vcf.gz"

    tabix -p vcf "${sampleName}.${callType}${prefix}.vcf.gz"

    if [ -e "${sampleName}.${callType}_filtered.vcf.gz" ]; then
        bcftools view -i '%FILTER="PASS" & INFO/PE>${variantSupport}' "${sampleName}.${callType}${prefix}.vcf.gz" -Oz -o "${sampleName}.${callType}${prefix}.pass.vcf.gz"
        tabix -p vcf "${sampleName}.${callType}${prefix}.pass.vcf.gz"
    fi
    """
}
