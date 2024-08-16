process ENSEMBLVEP_VEP {
    tag "$tumorName"
    publishDir  "${params.vep_outdir}", mode: 'copy'

    input:
    tuple val(tumorName), path(vcf)
    val  reference

    output:
    tuple val(tumorName), path("*.vcf.gz")  , optional:true, emit: vcf
    tuple val(tumorName), path("*.tab.gz")  , optional:true, emit: tab
    tuple val(tumorName), path("*.json.gz") , optional:true, emit: json
    path "*.html"                      , optional:true, emit: report
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def GenomeResources = [
            hg19: [
                vep_modules: "vep/105.0 tabix/0.2.6 vep-hg19-cache/105 hg19/p13",
                vcf2maf_modules: "vcf2maf/1.6.21b tabix/0.2.6 hg19/p13 vep-hg19-cache/105",
                vepCacheDir: "/.mounts/labs/gsi/modulator/sw/data/vep-hg19-cache-105/.vep",
                referenceFasta: "/.mounts/labs/gsi/modulator/sw/data/hg19-p13/hg19_random.fa",
                species: "homo_sapiens",
                ncbiBuild: "GRCh37",
                vepPath: "/.mounts/labs/gsi/modulator/sw/Ubuntu20.04/vep-105.0/bin/"
            ],
            hg38: [
                vep_modules: "vep/105.0 tabix/0.2.6 vep-hg38-cache/105 hg38/p12",
                vcf2maf_modules: "vcf2maf/1.6.21b tabix/0.2.6 hg38/p12 vep-hg38-cache/105",
                vepCacheDir: "/.mounts/labs/gsi/modulator/sw/data/vep-hg38-cache-105/.vep",
                referenceFasta: "/.mounts/labs/gsi/modulator/sw/data/hg38-p12/hg38_random.fa",
                species: "homo_sapiens",
                ncbiBuild: "GRCh38",
                vepPath: "/.mounts/labs/gsi/modulator/sw/Ubuntu20.04/vep-105.0/bin/"
            ],
            mm39: [
                vep_modules: "vep/105.0 tabix/0.2.6 vep-mm39-cache/105 mm39/p6",
                vcf2maf_modules: "vcf2maf/1.6.21b tabix/0.2.6 mm39/p6 vep-mm39-cache/105",
                vepCacheDir: "/.mounts/labs/gsi/modulator/sw/data/vep-mm39-cache-105/.vep",
                referenceFasta: "/.mounts/labs/gsi/modulator/sw/data/mm39-p6/mm39.fa",
                species: "mus_musculus",
                ncbiBuild: "GRCm39",
                vepPath: "/.mounts/labs/gsi/modulator/sw/Ubuntu20.04/vep-105.0/bin/"
            ]
        ]

    def file_extension = args.contains("--vcf") ? 'vcf' : args.contains("--json")? 'json' : args.contains("--tab")? 'tab' : 'vcf'
    def compress_cmd = args.contains("--compress_output") ? '' : '--compress_output bgzip'
    def prefix = task.ext.prefix ?: "${tumorName}"
    def dir_cache = "${GenomeResources[reference]['vepCacheDir']}"
    def reference_fasta = "${GenomeResources[reference]['referenceFasta']}"
    def genome = "${GenomeResources[reference]['ncbiBuild']}"
    def species = "${GenomeResources[reference]['species']}"
    def vep_modules = "${GenomeResources[reference]['vep_modules']}"
    def module_list = vep_modules.split(' ')
    def module_load_cmds = module_list.collect { module -> "module load ${module}" }.join('\n')
    def human_only_command_line = species == "homo_sapiens" ? "--polyphen b --af --af_1kg --af_esp --af_gnomad" : ""

    """
    ${module_load_cmds}

    vep \\
        --offline \\
        --dir $dir_cache \\
        -i $vcf \\
        -o ${prefix}.${file_extension}.gz \\
        $args \\
        $compress_cmd \\
        --fasta $reference_fasta \\
        --assembly $genome \\
        --species $species \\
        --no_progress --sift b --ccds --uniprot --hgvs --symbol --numbers --domains --gene_phenotype --mane  \\
        --allele_number --no_escape --xref_refseq --failed 1 --flag_pick_allele --pick_order canonical,tsl,biotype,rank,ccds,length \\
        $human_only_command_line \\
        --pubmed --regulatory \\
        --fork $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${tumorName}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    echo "" | gzip > ${prefix}.tab.gz
    echo "" | gzip > ${prefix}.json.gz
    touch ${prefix}_summary.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """
}