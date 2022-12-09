rule fastQC:
  input:
    "data/{read1}.fastq.qz",
    "data/{read2}.fastq.qz"
  output:
    "QC/
  
