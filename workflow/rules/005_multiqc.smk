# rule multiqc:
#     message: 
#         "Running MultiQC on raw data"
#     conda: 
#         "envs/001_QC.yml"
#     input:
#         expand("analysis/003_posttrim_qc/{sample}/{sample}_{lane}_{R}_trimmed_fastqc.html",
#                sample=sample_mrn, lane=lane, R=read),
#         expand("analysis/001_QC/{sample}/{sample}_{lane}_{R}_fastqc.html",
#                sample=sample_mrn, lane=lane, R=read),
#         expand("analysis/004_alignment/{tool}/{sample}_{lane}/{sample}_{lane}.bam.summary",
#                tool=tool, sample=sample_mrn, lane=lane)
#     output:
#         "analysis/multiqc_raw/{tool}"
#     log:
#         "logs/multiqc/{tool}/multiqc_raw.log"
#     benchmark:
#         "benchmarks/multiqc/{tool}/multiqc_raw.txt"
#     shell:
#         """
#         multiqc analysis/ \
#         -o {output} \
#         > {log} 2>&1
#         """


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
        expand("analysis/006_count/{tool}/{sample}_{lane}/{sample}_{lane}.counts",
               tool=count_tool, sample=sample_mrn, lane=lane)
    output:
        "analysis/multiqc_raw/{tool}/multiqc_report.html"
    log:
        "logs/multiqc/{tool}/multiqc_raw.log"
    benchmark:
        "benchmarks/multiqc/{tool}/multiqc_raw.txt"
    shell:
        """
        multiqc analysis/001_QC analysis/003_posttrim_qc logs \
        -o {output} \
        > {log} 2>&1
        """