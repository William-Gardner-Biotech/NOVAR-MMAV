# NOVAseq Reads-Metagenomic Mapping Against Viruses

This Nextflow pipeline processes Illumina sequencing data to analyze viral genomes. It includes:

- Merging sequencing lanes
- Removing non-template control (NTC) contamination (optional)
- Preprocessing reads (quality, filtering)
- Aligning reads to a viral reference
- Calculating overall average coverage (OVAL) per sample
- Generating a metrics report (CSV)

## Usage

1. Install Nextflow and dependencies
2. Update nextflow.config:
     - inputData: FASTQ file path
     - genomeOfInterest: Viral reference genome
     - virus_name: Output file prefix
     - NTC: Set to true if you have non-template controls
3. Run the pipeline:
   ```nextflow run main.nf```

## Output
- {virus_name}_BAMs: Mapped BAM files
- Combined_Reads.csv: Metrics (mapped reads, total reads, OVAL)
- oval_value.txt: OVAL per sample

## Dependencies
- Nextflow
- BBtools, samtools, python3
