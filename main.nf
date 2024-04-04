#!/opt/homebrew/Caskroom/miniconda/base/envs/aegis/bin/nextflow

nextflow.enable.dsl=2

workflow {

    println "Virus name: ${params.virus_name}"
    println "Reference Seq: ${params.genomeOfInterest}"



    // This is designed for multi-lane Illumina reads
    ch_fastqs_RPIPs = Channel.fromFilePairs ("${params.inputData}/**/*{L001,L002,L003,L004}_R1_001.fastq.gz", size:4)

    MERGE_LANES (
        ch_fastqs_RPIPs
    )

    // Check for "NTC" and then branch the run if condition met
    if (params.NTC != false) { 

        CONVERT_NTC (
            MERGE_LANES.out.fastq
        )

        REMOVE_NTC_CONTAMINANTS (
            CONVERT_NTC.out.fasta,
            MERGE_LANES.out.fastq
        )

        READ_PREPROCESS (
            REMOVE_NTC_CONTAMINANTS.out
        )
    }
    else {

        READ_PREPROCESS (
            MERGE_LANES.out.fastq
        )

    }

    MAPPING ( READ_PREPROCESS.out )

    CalculateOval ( MAPPING.out.stats )
        .set { ch_oval_plus_map_stats }

    NUM_M_READS ( ch_oval_plus_map_stats )

    COMBINE_READS ( NUM_M_READS.out.collect() )
}

process MERGE_LANES {

    tag "${sample_id}"

    cpus 4

    input:
    tuple val(sampleID_raw), path(reads)

    output:
    tuple path("${sample_id}.fastq.gz"), val(sample_id), emit: fastq
    val(sample_id), emit: sample_id

    script:
    String[] parts = sampleID_raw.split('-')
    // handles multiple runs having same name collision from CPC or NTC
    if (sampleID_raw.contains("NTC") || sampleID_raw.contains("CPC")) {
            sample_id =parts[0]+"_"+parts[1]
        }
    else {sample_id = parts[0]}

    """
    cat ${reads} > ${sample_id}.fastq.gz
    rm ${reads}
    """
}

process CONVERT_NTC {

    tag "${sample_id}"

    input:
    tuple path(input_seq), val(sample_id)

    output:
    path("${sample_id}.fasta"), emit: fasta

    when:
    sample_id.contains("NTC")

    script:
    """
    bbduk.sh in=${input_seq} out=duk_${input_seq} minlen=100 hdist=2 ftm=5 maq=10
    seqkit fq2fa duk_${input_seq} -o ${sample_id}.fasta
    dedupe.sh -Xmx8g in=${sample_id}.fasta out=${sample_id}_dedupe.fasta
    """
}

// Contamination filtering upstream of the trimming to save time
process REMOVE_NTC_CONTAMINANTS {

    tag "${sample_id}"

    cpus 4
    memory '8 GB'

    input:
    each path(ntc)
    tuple path(input_seq), val(sample_id)

    // unmapped reads are whatever doesn't map to the negative control AKA contamination
    output:
    tuple path("${sample_id}.fastq.gz"), val(sample_id)

    script:
    """
    bbmap.sh -Xmx8g perfectmode=f in=${input_seq} out=${sample_id}.sam ref=${ntc}
    reformat.sh unmappedonly=t in=${sample_id}.sam \
    ref=${ntc} \
    out=${sample_id}.fastq.gz
    """
}

// Picks up from the MERGE_FASTQs process and begins to process the files
process READ_PREPROCESS {

    tag "${sample_id}"

    input:
    tuple path(input_seq), val(sample_id)

    output:
    tuple path(pre_processed), val(sample_id)

    cpus 4

    when:
    !sample_id.contains("NTC")

    // publishDir "preprocessed", mode: 'copy'

    //BBDUK ref https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/

    // Take it further by k-trimming the ends
    //bbduk.sh in=processed.fastq out=processed_further.fastq ref=~/Bioinformatics/bbmap/resources/adapters.fa hdist=2 k=21 mink=11 ktrim=r
    script:
    pre_processed = ("pre_processed_" + sample_id)
    // Keep the type of pre_processed consitent so downstream processes don't get angry with improper methods.
    // new is a java method that will instantiate a new File/path variable

    """
    bbduk.sh in=${input_seq} out=${pre_processed} minlen=100 hdist=2 ftm=5 maq=10
    """
}

// bbduk is working as expected 27Oct23

/*
minlen = read length cutoff
hdist = how much memory to give k-mers to allow for more combinations
ftm = modulus divide the read to attempt to clean dirty ends
maq = mean quality score removal
*/

// output mapped is a bam file
process MAPPING {

    tag "${sample_id}"

    publishDir "${params.virus_name}_BAMs", mode: 'copy'

    input:
    tuple path(processed_seqs), val(sample_id)

    output:
    tuple path(mapped), val(sample_id), emit: mapped_reads // Bam file location
    tuple path(map_stats), val(sample_id), path(cov_stats), emit: stats // bbmap output log and covstat.txt file

    cpus 4

    memory '8 GB'

    //bbmap will generate me a pileup/coverage file that should be retrieved
    script: 
    mapped = sample_id + ".bam"
    map_stats = "STATS_" + sample_id + ".txt"
    cov_stats = sample_id + ".txt"

    """
    bbmap.sh -Xmx8g in=${processed_seqs} out=${sample_id}.sam ref=${params.genomeOfInterest} maxindel=100 minid=0.9 > ${map_stats} 2>&1 covstats=${cov_stats}
    samtools view -bS ${sample_id}.sam > ${mapped}
    rm ${sample_id}.sam
    """

    
    // 2>&1 moves the stderror to the stout position so I can capture the stats
    // Setting minid to 0.9 reduces specificity to hopefully capture more variants of Wuhan-1
    // maxindel set low because of viral genome and getting a huge insert wouldn't be captured by our read and the deletion really would inactivate virus of this size anyways
}
// bbmap is working as expected

process CalculateOval {
    
    tag "${sample_id}"

    input:
    tuple path(stats_file), val(sample_id), path(cov_stats)

    output:
    tuple path(stats_file), val(sample_id), path("oval_value.txt")

    // conditional that will perform the oval calc on segmented genomes
    script:
    """
    linecount=\$(awk 'END {print NR}' ${cov_stats})
    if [ "\$linecount" -eq 2 ]; then
        awk -F'\t' 'NR==2 {print \$2}' ${cov_stats} > oval_value.txt
    else
        python3 ${projectDir}/bin/oval.py -i ${cov_stats} > oval_value.txt
    fi
    """
}

// process extracts the number of reads mapped
// Split the process up to 
process NUM_M_READS {

    tag "${sample_id}"
    
    input:
    tuple path(stats_file), val(sample_id), path(oval_file)

    // not a tuple as it will be collected
    output:
    path output_file

    //publishDir "Indv_stats", mode: 'symlink'

    script:
    output_file=("${sample_id}_Reads.csv")

    """
    number=\$(grep "Mapped reads:" ${stats_file} | awk '{print \$3}')
    total_reads=\$(grep "Reads:" ${stats_file} | awk '{print \$2}')
    oval=\$(cat ${oval_file})
    echo "${sample_id},\$number,\$total_reads,\$oval" > ${output_file}
    """
}
// mixing bash variables into the script which must be escaped using "/"

process COMBINE_READS {

    input:
    path (reads_stats)

    output:
    path "${params.virus_name}_Reads.csv"

    publishDir "Combined_Reads", mode: 'copy'

    script:
    """
    echo "Accession,${params.virus_name}_mapped_reads,${params.virus_name}_total_reads,${params.virus_name}_overall_coverage" > ${params.virus_name}_Reads.csv
    cat ${reads_stats} >> ${params.virus_name}_Reads.csv
    """
}