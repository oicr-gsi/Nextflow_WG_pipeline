nextflow.enable.dsl = 2
process MERGE_AND_ZIP {

    input:
    path inputVcfs                    // Array of input VCF files
    path inputTbis                    // Array of input TBI files (not used in the command, but added for consistency)
    val sampleName = 'SAMPLE'         // Default sample name
    val callType = 'unmatched'        // Type of call, e.g., unmatched or somatic
    val modules = 'bcftools/1.9 vcftools/0.1.16 tabix/0.2.6' // Modules required
    val prefix = ''                   // Prefix for output files
    val variantSupport = 0            // Minimum variant support
    val jobMemory = 10                // Memory allocation in GB

    output:
    path "${sampleName}.${callType}${prefix}.vcf.gz" into dellyMergedVcf
    path "${sampleName}.${callType}${prefix}.vcf.gz.tbi" into dellyMergedTabixIndex
    path("${sampleName}.${callType}${prefix}.pass.vcf.gz").optional() into dellyMergedPassVcf
    path("${sampleName}.${callType}${prefix}.pass.vcf.gz.tbi").optional() into dellyMergedPassTabixIndex

    script:
    """
    set -eu -o pipefail

    # Load required modules
    module load ${modules}

    # Merge, sort, and compress VCF files
    vcf-concat ${inputVcfs.join(' ')} | vcf-sort | bgzip -c > "${sampleName}.${callType}${prefix}.vcf.gz"

    # Index the compressed VCF
    tabix -p vcf "${sampleName}.${callType}${prefix}.vcf.gz"

    # If a filtered VCF exists, create a filtered output
    if [ -e "${sampleName}.${callType}_filtered.vcf.gz" ]; then
        bcftools view -i '%FILTER="PASS" & INFO/PE>${variantSupport}' "${sampleName}.${callType}${prefix}.vcf.gz" -Oz -o "${sampleName}.${callType}${prefix}.pass.vcf.gz"
        tabix -p vcf "${sampleName}.${callType}${prefix}.pass.vcf.gz"
    fi
    """
}
