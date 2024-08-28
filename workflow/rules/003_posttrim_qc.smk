# A rule to perform post-trimming QC

rule posttrim_fastqc:
    message: "Running FastQC on trimmed data"
    conda: "envs/001_QC.yml"
    input:
        "analysis/002_trimming/{sample}/{sample}_{lane}_{R}_001_trimmed.fastq.gz"
    output:
        html="analysis/003_posttrim_qc/{sample}/{sample}_{lane}_{R}_001_trimmed_fastqc.html",
        zip="analysis/003_posttrim_qc/{sample}/{sample}_{lane}_{R}_001_trimmed_fastqc.zip"
    threads: config["threads"]
    params: 
        path=lambda wildcards: "analysis/003_posttrim_qc/{}".format(wildcards.sample)
    log:
        "logs/003_posttrim_qc/{sample}/{sample}_{lane}_{R}.log"
    benchmark:
        "benchmarks/003_posttrim_qc/{sample}/{sample}_{lane}_{R}.txt"
    shell:
        """
        mkdir -p {params.path}
        fastqc {input} \
        -t {threads} \
        -o {params.path} \
        > {log} 2>&1
        """


rule multiqc:
    message: "Running MultiQC on raw data"
    conda: "envs/001_QC.yml"
    input:
        expand("analysis/003_posttrim_qc/{sample}/{sample}_{lane}_{R}_001_trimmed_fastqc.html",
               sample=sample_mrn, lane=lane, R=read),
        expand("analysis/001_QC/{sample}/{sample}_{lane}_{R}_001_fastqc.html",
               sample=sample_mrn, lane=lane, R=read)
    output:
        directory("analysis/003_posttrim_qc/multiqc_raw")
    log:
        "logs/003_posttrim_qc/multiqc_raw.log"
    benchmark:
        "benchmarks/003_posttrim_qc/multiqc_raw.txt"
    shell:
        """
        multiqc analysis/001_QC/ \
        -o {output} \
        > {log} 2>&1
        """