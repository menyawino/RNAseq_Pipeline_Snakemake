# A rule to run FastQC on the raw data

rule raw_fastqc:
    message: 
        "Running FastQC on sample {wildcards.sample}_{lane}_{wildcards.R}"
    conda: 
        "envs/001_QC.yml"
    input:
        "samples/{sample}_{lane}_{R}_001.fastq.gz"
    output:
        html="analysis/001_QC/{sample}/{sample}_{lane}_{R}_fastqc.html",
        zip="analysis/001_QC/{sample}/{sample}_{lane}_{R}_fastqc.zip"
    threads: 
        config["np_threads"]
    params: 
        path=lambda wildcards: "analysis/001_QC/{}".format(wildcards.sample)
    log:
        "logs/001_QC/{sample}/{sample}_{lane}_{R}.log"
    benchmark:
        repeat("benchmarks/001_QC/{sample}/{sample}_{lane}_{R}.txt", config["benchmark"])
    shell:
        """
        mkdir -p {params.path}
        
        fastqc {input} \
        -t {threads} \
        -o {params.path} \
        > {log} 2>&1
        mv {params.path}/{wildcards.sample}_{wildcards.lane}_{wildcards.R}_001_fastqc.html {output.html}
        mv {params.path}/{wildcards.sample}_{wildcards.lane}_{wildcards.R}_001_fastqc.zip {output.zip}
        """