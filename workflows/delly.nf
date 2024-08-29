nextflow.enable.dsl = 2

include { PICARD_MARKDUPLICATES } from "../modules/mark_duplicates"
include { RUN_DELLY } from "../modules/run_delly.nf"

workflow delly {
    take:
        inputBams
        tumorName
        markdup
        reference
        picard_module

    main:
        def GenomeResources = [
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

        def markedBams = markdup ? PICARD_MARKDUPLICATES(inputBams, picard_module, tumorName) : inputBams

        def callType = (inputBams.size() == 1) ? 'unmatched' : 'somatic'
        def delly_modes = channel.from("DEL", "DUP", "INV", "INS", "BND")
        
        def delly_inputs = delly_modes
            .combine(channel.from(markedBams))
            .combine(channel.value(callType))
            .combine(channel.value(tumorName))
            .combine(channel.value(GenomeResources[reference]['rundelly_module']))
            .combine(channel.value(GenomeResources[reference]['rundelly_fasta']))
            .combine(channel.value(GenomeResources[reference]['rundelly_exclude_list']))

        RUN_DELLY(delly_inputs)
}
