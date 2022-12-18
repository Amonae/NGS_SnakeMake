

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
        "logs/bwa_index/{chromosome}.log",
    params:
        algorithm="bwtsw",
    wrapper:
        "v1.21.0/bio/bwa/index"      
        
        
  rule bwa_mem_sortsam:
    input:
        reads=["trimmed/{sample}.1.fastq.gz", "trimmed/{sample}.1.fastq.gz"],
        idx=multiext("index/{chromosome}", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        "mapped/{sample}.bam",
    log:
        "logs/bwa_mem/{sample}.log",
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
        bam="dedup/{sample}.bam",
        metrics="dedup/{sample}.metrics.txt",
    log:
        "logs/picard/dedup/{sample}.log",
    params:
        extra="--REMOVE_DUPLICATES true",
    resources:
        mem_mb=1024,
    wrapper:
        "v1.21.0/bio/picard/markduplicates" 
    
    
rule samtools_index:
    input:
        "dedup/{sample}.bam",
    output:
        "mapped/{sample}.sorted.bam.bai",
    log:
        "logs/samtools_index/{sample}.log",
    params:
        extra="",  # optional params string
    threads: 4 
    wrapper:
        "v1.21.0/bio/samtools/index"
        
    
 rule samtools_idxstats:
    input:
        bam="dedup/{sample}.bam",
        idx="mapped/{sample}.bam.bai",
    output:
        "mapped/{sample}.bam.idxstats",
    log:
        "logs/samtools/idxstats/{sample}.log",
    params:
        extra="",
    wrapper:
        "v1.21.0/bio/samtools/idxstats"
        
        
rule freebayes:
    input:
        ref="data/{chromosome}.fa.gz",
        samples="dedup/{sample}.bam",
        # the matching BAI indexes have to present for freebayes
        indexes="mapped/{sample}.bam.bai",
        # optional BED file specifying chromosomal regions on which freebayes
        # should run, e.g. all regions that show coverage
        #regions="path/to/region-file.bed"
    output:
        "calls/{sample}.vcf",  # either .vcf or .bcf
    log:
        "logs/freebayes/{sample}.log",
    params:
        extra="",  # optional parameters
        chunksize=100000,  # reference genome chunk size for parallelization (default: 100000)
        normalize=False,  # optional flag to use bcftools norm to normalize indels (Valid params are -a, -f, -m, -D or -d)
    threads: 2
    wrapper:
        "v1.21.0/bio/freebayes"
