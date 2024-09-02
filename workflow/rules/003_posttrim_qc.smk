# A rule to perform post-trimming QC

rule posttrim_fastqc:
    message: 
        "Running FastQC on trimmed sample {wildcards.sample}_{lane}_{R}"
    conda: 
        "envs/001_QC.yml"
    input:
        "analysis/002_trimming/{sample}/{sample}_{lane}_{R}_trimmed.fastq.gz"
    output:
        html="analysis/003_posttrim_qc/{sample}/{sample}_{lane}_{R}_trimmed_fastqc.html",
        zip="analysis/003_posttrim_qc/{sample}/{sample}_{lane}_{R}_trimmed_fastqc.zip"
    threads: 
        config["threads"]
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