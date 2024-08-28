# A rule to run stringtie on the aligned reads

rule stringtie_assembly:
    message: 
        "Building individual assembly for {wildcards.sample}_{lane}"
    input:
        bam="analysis/004_alignment/hisat2/{sample}_{lane}/{sample}_{lane}_Aligned.sortedByCoord.out.bam"
    output:
        "analysis/005_assembly/{sample}_{lane}/{sample}_{lane}.gtf"
    conda: 
        "envs/005_stringtie.yml"
    threads: 
        config["threads"]
    log:
        "logs/005_stringtie/assembly/{sample}_{lane}.log"
    benchmark:
        "benchmarks/005_stringtie/{sample}/{sample}_{lane}.txt"
    shell:
        """
        stringtie \
        {input.bam} \
        -o {output} \
        -p {threads} \
        2> {log}
        """

rule stringtie_merge:
    message: 
        "Merging the assembled transcripts"
    input:
        gtf=expand("analysis/005_assembly/{sample}_{lane}/{sample}_{lane}.gtf", sample=sample_mrn, lane=lane)
    output:
        "analysis/005_assembly/merged.gtf"
    conda: 
        "envs/005_stringtie.yml"
    threads: 
        config["threads"]
    log:
        "logs/005_stringtie/merge/merged.log"
    benchmark:
        "benchmarks/005_stringtie/merged.txt"
    shell:
        """
        stringtie \
        --merge \
        -o {output} \
        -p {threads} \
        {input.gtf} \
        2> {log}
        """