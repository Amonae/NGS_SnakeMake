

rule fastqc:
    input:
        "data/{read}.fastq.gz"
    output:
        html="QC/{read}.html",
        zip="QC/{read}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params:
        "--quiet"
    log:
        "logs/QC/{read}.log"
    threads: 1
    wrapper:
        "v1.20.0/bio/fastqc"


rule trimmomatic:
    input:
        r1="data/{sample}_R1.fastq.gz",
        r2="data/{sample}_R2.fastq.gz"
    output:
        r1="trimmed/{sample}.1.fastq.gz",
        r2="trimmed/{sample}.2.fastq.gz",
        # reads where trimming entirely removed the mate
        r1_unpaired="trimmed/{sample}.1.unpaired.fastq.gz",
        r2_unpaired="trimmed/{sample}.2.unpaired.fastq.gz"
    log:
        "logs/trimmomatic/{sample}.log"
    params:
        # list of trimmers (see manual)
        trimmer=["TRAILING:26"],
        # optional parameters
        extra="",
        compression_level="-9"
    threads:
        32
    resources:
        mem_mb=1024
    wrapper:
        "v1.20.0/bio/trimmomatic/pe"
        
        
  rule bwa_index:
    input:
        "data/{chromosome}.fa.gz",
    output:
        idx=multiext("index/{chromosome}", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa/index/{chromosome}.log",
    params:
        algorithm="bwtsw",
    wrapper:
        "v1.21.0/bio/bwa/index"      
        
        
  rule bwa_mem_sortsam:
    input:
        reads=["trimmed/{sample}.1.fastq.gz", "trimmed/{sample}.2.fastq.gz"],
        idx=multiext("index/{chromosome}", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        "mapped/{sample}.bam",
    log:
        "logs/bwa/mem_align/{sample}.log",
    params:
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sorting="samtools",  # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
        sort_extra="",  # Extra args for samtools/picard.
    threads: 8
    wrapper:
        "v1.21.0/bio/bwa/mem"
   
   
 rule mark_duplicates:
    input:
        bams="mapped/{sample}.bam",
    output:
        bam="mapped/{sample}_dedup.bam",
        metrics="mapped/{sample}_dedup.metrics.txt",
    log:
        "logs/dedup/{sample}.log",
    params:
        extra="--REMOVE_DUPLICATES true",
    resources:
        mem_mb=1024,
    wrapper:
        "v1.21.0/bio/picard/markduplicates" 
    
    
rule samtools_index:
    input:
        "mapped/{sample}_dedup.bam",
    output:
        "mapped/{sample}_dedup.sorted.bam.bai",
    log:
        "logs/samtools/index/{sample}.log",
    params:
        extra="",  # optional params string
    threads: 4 
    wrapper:
        "v1.21.0/bio/samtools/index"
        
    
 rule samtools_idxstats:
    input:
        bam="mapped/{sample}_dedup.bam",
        idx="mapped/{sample}_dedup.bam.bai",
    output:
        "mapped/{sample}.bam.idxstats",
    log:
        "logs/samtools/idxstats/{sample}.log",
    params:
        extra="",
    wrapper:
        "v1.21.0/bio/samtools/idxstats"
     
     
rule normalize_fasta
    input: 
        "data/{chromosome}.fa.gz"
    output: 
        "data/{chromosome}_normalized.fa.gz"
    shell:
        "picard NormalizeFasta I={input} O={output}"
rule unzip:
    input:
        "data/{sample}.fa.gz"
    output:
        "data/{sample}.fa"
    shell:
        "gunzip -c {input} >{output}"
        
        
rule samtools_fai:
    input:
       "data/{sample}.fa"
    output:
        "data/{sample}.fa.fai",
    log:
        "logs/samtools/index/{sample}.log",
    params:
        extra="",  # optional params string
    wrapper:
        "v1.21.0/bio/samtools/faidx"

rule freebayes:
    input:
        ref="data/{chromosome}_normalized.fa.gz",
        samples="mapped/{sample}_dedup.bam",
        # the matching BAI indexes have to present for freebayes
        indexes="mapped/{sample}_dedup.bam.bai",
    output:
        "calls/{sample}.vcf",  # either .vcf or .bcf
    log:
        "logs/freebayes/{sample}_dedup.log",
    params:
        extra="",  # optional parameters
        chunksize=100000,  # reference genome chunk size for parallelization (default: 100000)
        normalize=False,  # optional flag to use bcftools norm to normalize indels (Valid params are -a, -f, -m, -D or -d)
    threads: 2
    wrapper:
        "v1.21.0/bio/freebayes"
