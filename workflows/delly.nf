nextflow.enable.dsl=2
include {PICARD_MARKDUPLICATES} from "./mark_duplicates"

workflow delly {

    take:
    inputTumor
    inputNormal
    markdup
    outputFileNamePrefix
    reference
    picard_module

    main:

    GenomeResources = [
        hg19: [
            rundelly_module: "delly/0.9.1 bcftools/1.9 tabix/0.2.6 hg19/p13 hg19-delly/1.0",
            rundelly_fasta: "/.mounts/labs/gsi/modulator/sw/data/hg19-p13/hg19_random.fa",
            rundelly_exclude_list: "/.mounts/labs/gsi/modulator/sw/data/hg19-delly-1.0/human.hg19.excl.tsv"
        ],
        hg38: [ 
            rundelly_module: "delly/0.9.1 bcftools/1.9 tabix/0.2.6 hg38/p12 hg38-delly/1.0",
            rundelly_fasta: "/.mounts/labs/gsi/modulator/sw/data/hg38-p12/hg38_random.fa",
            rundelly_exclude_list: "/.mounts/labs/gsi/modulator/sw/data/hg38-delly-1.0/human.hg38.excl.tsv"
        ]
    ]

    inputBams = inputTumor.mix(inputNormal)
    if markdup {
        PICARD_MARKDUPLICATES (
            inputBam
            picard_module         
        ) 
    }
}