import os

configfile: "RNA_seq_config.yaml" 

# software
hisat2=config["software"]["hisat2"]
stringtie=config["software"]["stringtie"]
samtools=config["software"]["samtools"]

# genome
genome=config["genome"]

# raw data
RNAseq_path=config["RNAseq_path"]
gtfpath=config["gtfpath"]

# output prefix
outputname=config["outputname"]

# According to the unsure of the hisat2-build,
# highly recommand to build the index for your genome
# or do not use this pipline
indexp=config["indexp"]
index_f_list=os.listdir(indexp)
num = 0
for f in index_f_list:
    if f.endswith(".1.ht2"):
        index = indexp + f.strip(".1.ht2")
    if f.endswith("ht2"):
        num += 1

# threads
THREAD=config["THREADS"]

# samplelist
samplelist=list(set([i.split("_")[0] for i in os.listdir(RNAseq_path)]))

rule all:
    input: 
        expand("C400.{SAMPLE}.transcripts.round2.gtf", SAMPLE=samplelist),

rule Hisat2:
    input: 
        f1=RNAseq_path+"{sample}_1.clean.fq.gz",
        f2=RNAseq_path+"{sample}_2.clean.fq.gz",
    output: 
        sam=outputname + ".{sample}.hisat2.sam"
    params:
        i=index
    threads:
        THREAD
    shell:
        hisat2+" --dta  -p {threads} --max-intronlen 5000000 -x {params.i} -1 {input.f1} -2 {input.f2} -S {output.sam}"

rule Samtools:
    input:
        sam=outputname + ".{sample}.hisat2.sam"
    output:
        bam=outputname + ".{sample}.sorted.bam",
        bai=outputname + ".{sample}.sorted.bam.bai",
    params:
        temp=outputname+".{sample}.hisat2.sorted.tmp"
    threads:
        THREAD
    shell:
        "{samtools} view -F 4 -Su {input.sam} -@ {threads} | "
        "{samtools} sort -T {params.temp} -o {output.bam} -@ {threads} |"
        "{samtools} index {output.bam}"

rule r1_Stringtie:
    input: 
        bam=outputname + ".{sample}.sorted.bam",
        bai=outputname + ".{sample}.sorted.bam.bai"
    output: 
        gtf=outputname + ".{sample}.transcripts.gtf",
    params:
        n=outputname,
        refgtf=gtfpath
    threads:
        THREAD
    shell: 
        "{stringtie} {input.bam} -G {params.refgtf} -l {params.n} -o {output.gtf} -p {threads}"

rule create_mergelist:
    input: gtf=expand(outputname + ".{SAMPLE}.transcripts.gtf", SAMPLE=samplelist)
    output: "merge.list"
    shell:
        "ls {input.gtf} >> merge.list"

rule Stringtie_merge:
    input:
        list="merge.list"
    output:
        outgft="StringTie_merged.gtf"
    params:
        refgtf=gtfpath
    shell:
        "{stringtie} --merge -G {params.refgtf} -F 0.1 -T 0.1 -i -o {output.outgtf} {input.list}"

rule r2_Stringtie:
    input: 
        bam=outputname + ".{sample}.sorted.bam",
        bai=outputname + ".{sample}.sorted.bam.bai",
        merged_gft="StringTie_merged.gtf"
    output: 
        gtf=outputname + ".{sample}.transcripts.round2.gtf",
    params:
        n=outputname
    threads:
        THREAD
    shell: 
        "{stringtie} {input.bam} -G {input.merged_gft} -l {params.n} -o {output.gtf} -p {threads}"
