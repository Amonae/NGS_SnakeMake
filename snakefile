

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
        
  
