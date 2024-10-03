# A rule to merge all lanes of trimmed reads for R1 and R2 then delete the original files

rule merge_reads:
    message:
        "Merging trimmed reads for sample {wildcards.sample}"
    input:
        fq1=lambda wildcards: expand("analysis/002_trimming/{sample}/{sample}_L{lane}_R1_trimmed.fastq.gz", sample=wildcards.sample, lane=["001", "002", "003", "004"]),
        fq2=lambda wildcards: expand("analysis/002_trimming/{sample}/{sample}_L{lane}_R2_trimmed.fastq.gz", sample=wildcards.sample, lane=["001", "002", "003", "004"])
    output:
        fq1_merged="analysis/002_trimming/{sample}/{sample}_merged_R1_trimmed.fastq.gz",
        fq2_merged="analysis/002_trimming/{sample}/{sample}_merged_R2_trimmed.fastq.gz"
    threads:
        config["threads"]
    log:
        "logs/002_trimming/{sample}_merge.log"
    benchmark:
        repeat("benchmarks/002_trimming/merge/{sample}_merge.txt", config["benchmark"])
    shell:
        """
        cat {input.fq1} > {output.fq1_merged} && \
        cat {input.fq2} > {output.fq2_merged} && \
        rm {input.fq1} {input.fq2}
        """


# A rule to perform post-trimming QC

rule posttrim_fastqc:
    message: 
        "Running FastQC on trimmed sample {wildcards.sample}"
    conda: 
        "envs/001_QC.yml"
    input:
        "analysis/002_trimming/{sample}/{sample}_merged_R{R}_trimmed.fastq.gz"
    output:
        html="analysis/003_posttrim_QC/{sample}/{sample}_merged_R{R}_trimmed_fastqc.html",
        zip="analysis/003_posttrim_QC/{sample}/{sample}_merged_R{R}_trimmed_fastqc.zip"
    threads:
        config["np_threads"]
    params:
        path=lambda wildcards: "analysis/003_posttrim_QC/{}".format(wildcards.sample)
    log:
        "logs/003_posttrim_QC/{sample}/{sample}_merged_R{R}_trimmed.log"
    benchmark:
        repeat("benchmarks/003_posttrim_QC/{sample}/{sample}_merged_R{R}_trimmed.txt", config["benchmark"])
    shell:
        """
        mkdir -p {params.path}

        fastqc {input} \
        -t {threads} \
        -o {params.path} \
        > {log} 2>&1
        mv {params.path}/{wildcards.sample}_merged_R{R}_trimmed_fastqc.html {output.html}
        mv {params.path}/{wildcards.sample}_merged_R{R}_trimmed_fastqc.zip {output.zip}
        """