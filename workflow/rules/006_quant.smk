rule stringtie_count:
    message: 
        "Counting reads in sample {wildcards.sample}_{lane} with stringtie"
    input:
        gtf="analysis/005_assembly/merged.gtf",
        bam=expand("analysis/004_alignment/hisat2/{sample}_{lane}/{sample}_{lane}_Aligned.sortedByCoord.out.bam",
            sample=sample_mrn, lane=lane)
    output:
        "analysis/006_count/stringtie/{sample}_{lane}/{sample}_{lane}.counts"
    conda: 
        "envs/005_stringtie.yml"
    threads: 
        config["threads"]
    log:
        "logs/006_count/stringtie/{sample}_{lane}.log"
    benchmark:
        "benchmarks/006_count/stringtie/{sample}/{sample}_{lane}.txt"
    shell:
        """
        stringtie \
        -e \
        -B \
        -p {threads} \
        -G {input.gtf} \
        -o {output} \
        {input.bam} \
        2> {log}
        """


rule kallisto_count:
    message:
        "Counting reads in sample {wildcards.sample}_{lane} with kallisto"
    input:
        fastq1="analysis/002_trimming/{sample}/{sample}_{lane}_R1_001_trimmed.fastq.gz",
        fastq2="analysis/002_trimming/{sample}/{sample}_{lane}_R2_001_trimmed.fastq.gz"
    output:
        directory("analysis/006_count/kallisto/{sample}_{lane}")
    conda:
        "envs/006_kallisto.yml"
    threads:
        config["threads"]
    params:
        index=config["aligner"]["index_kallisto"]
    log:
        "logs/006_count/kallisto/{sample}_{lane}.log"
    benchmark:
        "benchmarks/006_count/kallisto/{sample}_{lane}.txt"
    shell:
        """
        kallisto quant \
        -i {params.index} \
        -o {output} \
        -t {threads} \
        -b 100 \
        {input.fastq1} {input.fastq2} \
        2> {log}
        """