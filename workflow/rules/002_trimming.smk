# A rule to trim and remove adapters from the raw data using fastp

rule trimming:
    message: "Trimming and removing adapters from raw data"
    conda: "envs/001_QC.yml"
    input:
        fw="samples/{sample}_{lane}_R1_001.fastq.gz",
        rv="samples/{sample}_{lane}_R2_001.fastq.gz"
    output:
        fw="analysis/002_trimming/{sample}/{sample}_{lane}_R1_001_trimmed.fastq.gz",
        rv="analysis/002_trimming/{sample}/{sample}_{lane}_R2_001_trimmed.fastq.gz",
        html="analysis/002_trimming/{sample}/{sample}_{lane}_fastp.html"
    threads: config["threads"]
    log:
        "logs/002_trimming/{sample}/{sample}_{lane}.log"
    benchmark:
        "benchmarks/002_trimming/{sample}/{sample}_{lane}.txt"
    shell:
        """
        fastp \
        -i {input.fw} \
        -I {input.rv} \
        -o {output.fw} \
        -O {output.rv} \
        -j {log} \
        -w {threads} \
        -h {output.html} \
        --adapter_sequence TACACTCTTTCCCTACACGACGCTCTTCCGATCT \
        --adapter_sequence_r2 GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT \
        > {log} 2>&1
        """