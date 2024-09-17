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
        expand("analysis/006_count/{tool}/{sample}_{lane}",
               tool=count_tool, sample=sample_mrn, lane=lane)
    output:
        "analysis/007_multiqc/{tool}/multiqc_report.html"
    params:
        output=lambda wildcards: "analysis/007_multiqc/{}".format(wildcards.tool)
        
    log:
        "logs/multiqc/{tool}/multiqc_raw.log"
    benchmark:
        "benchmarks/multiqc/{tool}/multiqc_raw.txt"
    shell:
        """
        multiqc analysis/001_QC analysis/003_posttrim_qc logs \
        -o {params.output} \
        > {log} 2>&1
        """