READS = ["example_R1", "example_R2"]

rule fastqc:
    input:
        expand("data/{reads}.fastq.gz", reads=READS)
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
