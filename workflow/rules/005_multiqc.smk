rule multiqc:
    message: 
        "Running MultiQC on raw data"
    conda: 
        "envs/001_QC.yml"
    input:
        expand("analysis/003_posttrim_qc/{sample}/{sample}_{lane}_{R}_trimmed_fastqc.html",
               sample=sample_mrn, lane=lane, R=read),
        expand("analysis/001_QC/{sample}/{sample}_{lane}_{R}_fastqc.html",
               sample=sample_mrn, lane=lane, R=read),
        expand("analysis/004_alignment/hisat2/{sample}_{lane}/{sample}_{lane}.bam.summary",
               sample=sample_mrn, lane=lane)
    output:
        directory("analysis/multiqc_raw")
    log:
        "logs/003_posttrim_qc/multiqc_raw.log"
    benchmark:
        "benchmarks/003_posttrim_qc/multiqc_raw.txt"
    shell:
        """
        multiqc analysis/ \
        -o {output} \
        > {log} 2>&1
        """