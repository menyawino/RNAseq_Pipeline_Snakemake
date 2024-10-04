# A rule to run kallisto for read counting

rule stringtie_count:
    message:
        "Counting reads in sample {wildcards.sample} with stringtie"
    input:
        gtf="analysis/005_assembly/merged.gtf",
        bam="analysis/004_alignment/hisat2/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        "analysis/006_count/stringtie/{sample}/"
    conda:
        "envs/005_stringtie.yml"
    threads:
        config["threads"]
    # stringtie requires a lot of memory which causes the cluster to kill the job, this would ensure there is enough memory before running the job
    params:
        output=lambda wildcards: "analysis/006_count/stringtie/{}/".format(wildcards.sample)
    resources:
        mem_mb=20000
    log:
        "logs/006_count/stringtie/{sample}.log"
    benchmark:
        repeat("benchmarks/006_count/stringtie/{sample}/{sample}.txt", config["benchmark"])
    shell:
        """
        stringtie \
        -e \
        -B \
        -p {threads} \
        -G {input.gtf} \
        -o {params.output} \
        {input.bam} \
        2> {log}
        """

# A rule to run kallisto for read counting

rule kallisto_count:
    message:
        "Counting reads in sample {wildcards.sample} with kallisto"
    input:
        fastq1="analysis/002_trimming/{sample}/{sample}_R1_trimmed.fastq.gz",
        fastq2="analysis/002_trimming/{sample}/{sample}_R2_trimmed.fastq.gz"
    output:
        directory("analysis/006_count/kallisto/{sample}")
    conda:
        "envs/006_kallisto.yml"
    threads:
        config["threads"]
    params:
        index=config["aligner"]["index_kallisto"],
        # output dir with sample and
        output=lambda wildcards: "analysis/006_count/kallisto/{}".format(wildcards.sample)
    log:
        "logs/006_count/kallisto/{sample}/{sample}.log"
    benchmark:
        repeat("benchmarks/006_count/kallisto/{sample}/{sample}.txt", config["benchmark"])
    shell:
        """
        mkdir -p {params.output}
        
        kallisto quant \
        -i {params.index} \
        -o {params.output} \
        -t {threads} \
        -b 100 \
        {input.fastq1} {input.fastq2} \
        2> {log}
        """