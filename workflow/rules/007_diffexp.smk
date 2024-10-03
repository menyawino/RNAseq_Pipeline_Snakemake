# A rule to run Ballgown for differential expression analysis

rule ballgown_diffexp:
    message:
        "Running Ballgown for differential expression analysis"
    input:
        "analysis/006_count/stringtie/{sample}/{sample}.counts"
    output:
        "results/diffexp/stringtie/{sample}.diffexp.txt"
    conda:
        "envs/007_ballgown.yml"
    log:
        "logs/diffexp/{sample}.log"
    benchmark:
        repeat("benchmarks/diffexp/{sample}.txt", config["benchmark"])
    script:
        "scripts/ballgown.R"


# A rule to run Sleuth for differential expression analysis

rule sleuth_analysis:
    message:
        "Running Sleuth for differential expression analysis"
    input:
        samples=expand("analysis/006_count/kallisto/{sample}", sample=sample_mrn),
        script="workflow/scripts/sleuth.R"
    output:
        "results/sleuth/differential_expression_results.tsv",
        "results/sleuth/sleuth_report.pdf"
    conda:
        "envs/sleuth.yml"
    log:
        "logs/007_sleuth_analysis.log"
    benchmark:
        repeat("benchmarks/007_sleuth_analysis.txt", config["benchmark"])
    shell:
        """
        Rscript {input.script} \
        --output {output[0]} \
        --report {output[1]} \
        --log {log}
        """