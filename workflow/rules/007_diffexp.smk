rule ballgown_diffexp:
    message:
        "Running Ballgown for differential expression analysis"
    input:
        "analysis/006_count/stringtie/{sample}/{sample}_{lane}.counts"
    output:
        "results/diffexp/stringtie/{sample}_{lane}.diffexp.txt"
    conda:
        "envs/007_ballgown.yml"
    log:
        "logs/diffexp/{sample}_{lane}.log"
    benchmark:
        "benchmarks/diffexp/{sample}_{lane}.txt"
    script:
        "scripts/ballgown.R"


rule sleuth_analysis:
    message:
        "Running Sleuth for differential expression analysis"
    input:
        samples=expand("analysis/006_count/kallisto/{sample}_{lane}", sample=sample_mrn, lane=lane),
        script="workflow/scripts/sleuth.R"
    output:
        "results/sleuth/differential_expression_results.tsv",
        "results/sleuth/sleuth_report.pdf"
    conda:
        "envs/sleuth.yml"
    log:
        "logs/007_sleuth_analysis.log"
    benchmark:
        "benchmarks/007_sleuth_analysis.txt"
    shell:
        """
        Rscript {input.script} \
        --output {output[0]} \
        --report {output[1]} \
        --log {log}
        """