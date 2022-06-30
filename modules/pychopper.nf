process PYCHOPPER {

    publishDir "results/${params.out_dir}/pychopper/"

    label "large"

    input:
        tuple val(id), path(fastq)
        path(txt)
        val(cdna_kit)

    output:
        val "$id", emit: id
        path "${id}_pychop.fq", emit: fastq
        path "$txt", emit: txt
        path "*pychopper*", emit: multiqc

    script:
        """
        cdna_classifier.py -t 20 \
            -k $cdna_kit \
            -r "${id}_pychopper_report.pdf" \
            -u "${id}_pychopper.unclassified.fq" \
            -w "${id}_pychopper.rescued.fq" \
            -S "${id}_pychopper.stats" \
            -A "${id}_pychopper.scores" \
            "${fastq}" "${id}_pychop.fq"
        """
}
