# A rule to run FastQC on the raw data

rule raw_fastqc:
    message: "Running FastQC on raw data"
    conda: "envs/001_QC.yml"
    input:
        "samples/{sample}_{lane}_{R}_001.fastq.gz"
    output:
        html="analysis/001_QC/{sample}/{sample}_{lane}_{R}_001_fastqc.html",
        zip="analysis/001_QC/{sample}/{sample}_{lane}_{R}_001_fastqc.zip"
    threads: config["threads"]
    params: 
        path=lambda wildcards: "analysis/001_QC/{}".format(wildcards.sample)
    log:
        "logs/001_QC/{sample}/{sample}_{lane}_{R}.log"
    benchmark:
        "benchmarks/001_QC/{sample}/{sample}_{lane}_{R}.txt"
    shell:
        """
        mkdir -p {params.path}
        fastqc {input} \
        -t {threads} \
        -o {params.path} \
        > {log} 2>&1
        """