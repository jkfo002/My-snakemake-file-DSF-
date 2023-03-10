rule all:
    input:
        expand("{SAMPLE}_{REP}.{PAIR}.out", SAMPLE=["A", "B"], REP=["1", "2"], PAIR=["1", "2"])

rule A:
    input:
        "data/{sample}.txt" 
    output:
        "out/{sample}.txt"
    shell:
        "cat {input} > {output}"

##########################################################################################################
# 正则表达式选择文件执行
rule all:
    input:
        expand("output/A_{REP}.{PAIR}.out", REP=["1", "2"], PAIR=["1", "2"])

rule A:
    input:
        "data/A_{sample}.txt"
    output:
        "output/A_{sample}.out"
    shell:
        "cat {input} > {output}"

##########################################################################################################
# document: https://snakemake.readthedocs.io/en/stable/index.html
# {wildcards.wildcards_name} 仅用于shell：https://endrebak.gitbooks.io/the-snakemake-book/content/chapters/wildcards/wildcards.html
# this is a test
# this is a test