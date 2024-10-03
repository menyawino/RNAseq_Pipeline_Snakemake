# A rule to run stringtie on the aligned reads

rule stringtie_assembly:
    message: 
        "Building individual assembly for {wildcards.sample}"
    input:
        bam="analysis/004_alignment/hisat2/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        "analysis/005_assembly/{sample}/{sample}.gtf"
    conda: 
        "envs/005_stringtie.yml"
    threads: 
        config["threads"]
    log:
        "logs/005_stringtie/assembly/{sample}.log"
    benchmark:
        repeat("benchmarks/005_stringtie/{sample}/{sample}.txt", config["benchmark"])
    shell:
        """
        stringtie \
        {input.bam} \
        -o {output} \
        -p {threads} \
        2> {log}
        """

# A rule to merge the assembled transcripts

rule stringtie_merge:
    message: 
        "Merging the assembled transcripts"
    input:
        gtf=expand("analysis/005_assembly/{sample}/{sample}.gtf", sample=sample_mrn)
    output:
        "analysis/005_assembly/merged.gtf"
    conda: 
        "envs/005_stringtie.yml"
    threads: 
        config["threads"]
    log:
        "logs/005_stringtie/merge/merged.log"
    benchmark:
        repeat("benchmarks/005_stringtie/merged.txt", config["benchmark"])
    shell:
        """
        stringtie \
        --merge \
        -o {output} \
        -p {threads} \
        {input.gtf} \
        2> {log}
        """