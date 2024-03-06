#!/usr/bin/env nextflow
/*
* USAGE: nextflow run myscript.nf -qs 8
* This script creates hard links to data that exists in nextflows work directory.
*/

 
// Path to raw fastq files
Channel
    .fromFilePairs("${params.rawDataPath}", flat: true)
    .ifEmpty {error "Cannot find any reads matching: ${params.reads}"}
    .set {reads}



/* Code that should not change - STARTS HERE */

params.trimPath = "${params.outputDir}/trimmomatic"
params.kraken2Path = "${params.outputDir}/kraken2"
params.fastqcPath = "${params.outputDir}/fastqc"
params.alignPath = "${params.outputDir}/bowtie2"
params.viremaPath = "${params.outputDir}/Virema"
params.depthPath = "${params.outputDir}/depth"
params.varPath = "${params.outputDir}/var"
params.conPath = "${params.outputDir}/consensus"

/*
* Step 1. Trimming
* WARNING: considers '1' a valid exit status to get around wrapper error
*/

process trimmomatic {
    cpus params.trimCPU
    memory "$params.trimMemory GB"
    publishDir params.trimPath
    validExitStatus 0,1

    input:
    set val(id), file(read1), file(read2) from reads

    output:
    set val(id), "${read1.baseName}_qualtrim_paired.fastq.gz", "${read2.baseName}_qualtrim_paired.fastq.gz" into kraken2Channel
    set val(id), "${read1.baseName}_qualtrim_paired.fastq.gz", "${read2.baseName}_qualtrim_paired.fastq.gz" into extractChannel
    file "*_qualtrim_unpaired.fastq.gz"
    stdout trim_out

    """
    trimmomatic PE \
    -threads $params.trimCPU -phred33 $read1 $read2 \
    ${read1.baseName}_qualtrim_paired.fastq.gz ${read1.baseName}_qualtrim_unpaired.fastq.gz \
    ${read2.baseName}_qualtrim_paired.fastq.gz ${read2.baseName}_qualtrim_unpaired.fastq.gz \
    $params.trimOptions
    """
}

/*
* Step 2. Kraken2
*/

process kraken2 {
    cpus params.kraken2CPU
    memory "${params.kraken2Memory} GB"
    publishDir params.kraken2Path
    validExitStatus 0

    input:
    set val(id), file(read1), file(read2) from kraken2Channel

    output:
    file "*_report.txt"
    file "*_kraken.txt"
    set val(id), "${read1.simpleName}_iav.fastq.gz", "${read2.simpleName}_iav.fastq.gz" into fastqcChannel
    set val(id), "${read1.simpleName}_iav.fastq.gz", "${read2.simpleName}_iav.fastq.gz" into catChannel

    """
    kraken2 \
	--db $params.kraken2DB \
	--report ${read1.baseName}_report.txt \
	--gzip-compressed \
	--paired \
	$read1 \
	$read2 \
	> ${read1.simpleName}_kraken.txt
    
    python3 ${params.kraken2ScriptPath}/extract_kraken_reads.py \
	-k ${read1.simpleName}_kraken.txt \
	-s1 $read1 \
	-s2 $read2 \
	-o ${read1.simpleName}_iav.fastq \
	-o2 ${read2.simpleName}_iav.fastq \
	--fastq-output \
	--include-children \
	-t 11320 \
	-r ${read1.baseName}_report.txt \
	--max 1000000000 

    gzip ${read1.simpleName}_iav.fastq
    gzip ${read2.simpleName}_iav.fastq
    """

}

/*
* Step 2. FASTQC of trimmed reads
*/
process runFASTQC {
    cpus 1
    memory "$params.fastqcMemory GB"
    publishDir params.fastqcPath

    input:
    set pair_id, file(read1), file(read2) from fastqcChannel


    output:
    file "*.html"
    file "*.zip"

    """
    fastqc -t 2 -o ./ --noextract $read1 $read2
    """
}

/*
* Step 3. Combine FASTQ pairs
*/
process combineFASTQ {
    publishDir params.trimPath

    input:
    set pair_id, file(read1), file(read2) from catChannel

    output:
    file  "*_qualtrim_paired_iav.fastq.gz" into bowtie2Channel

    """
    cat $read1 $read2 > ${pair_id}_qualtrim_paired_iav.fastq.gz
    """
}

/*
* Step 4. Bowtie2 global alignment
*/
process runbowtie2 {
    cpus params.bowtie2CPU
    memory "$params.bowtie2Mem GB"
    publishDir params.alignPath

    input:
    file in_cat from bowtie2Channel

    output: 
    file "*_unaligned.fastq" into viremaChannel
    file "*.bam" into depthChannel, vcfChannel

    """
    bowtie2 \
            -p $params.bowtie2CPU \
            -x $params.bowtie2_index \
            -U $in_cat \
            --score-min $params.scoreMin \
            --un ${in_cat.simpleName}_unaligned.fastq | \
        samtools sort -n -@ $params.bowtie2CPU -O sam - |  \
        samtools fixmate -@ $params.bowtie2CPU -O sam -m - - | \
        samtools sort -@ $params.bowtie2CPU -O sam | \
        samtools markdup -@ $params.bowtie2CPU -r -s -O bam - ${in_cat.simpleName}.bam
    """
}

/*
Step 6. Count depth
*/
process depth {
    cpus params.bowtie2CPU
    memory "$params.bowtie2Mem GB"
    publishDir params.depthPath, mode: 'copy'

    input:
    file in_bam from depthChannel    

    output:
    file "*_dp.tsv" 


    """
    samtools depth -aa  ${in_bam} > ${in_bam.simpleName}_dp.tsv
    """
}


/*
Step 7. Call variants
*/
process vcf {
    cpus params.bowtie2CPU
    memory "$params.bowtie2Mem GB"
    publishDir params.varPath, mode: 'copy'

    input:
    file in_bam from vcfChannel    

    output:
    file "${in_bam.simpleName}.vcf.gz" into vcfgzChannel
    file "${in_bam.simpleName}.vcf.gz.tbi" into tbiChannel
    file "*.vcf"

    """
    bcftools mpileup \
            --threads $params.bowtie2CPU \
            -A \
            -d 1000 \
            -q 20 \
            -Q 20 \
            -I \
            -f ${params.bowtie2_index}.fasta ${in_bam} |\
        bcftools call \
            -mv \
            --threads $params.bowtie2CPU \
            -Ov | \
        bcftools norm \
            -m- \
            -Ov \
            -o ${in_bam.simpleName}.vcf
    
    bgzip -c ${in_bam.simpleName}.vcf > ${in_bam.simpleName}.vcf.gz
    tabix ${in_bam.simpleName}.vcf.gz
    """
}

/*
Step 8. Individual level consensus
*/
process consensus {
    cpus params.bowtie2CPU
    memory "$params.bowtie2Mem GB"
    publishDir params.conPath, mode: 'copy'

    input:
    file vcf_gz from vcfgzChannel
    file vcf_gz_tbi from tbiChannel

    output:
    file "*.ebwt" into viremaIndexChannel
    file "*.fasta"

    """
    bash $params.conApp ${params.bowtie2_index}.fasta  $vcf_gz > ${vcf_gz.simpleName}.fasta

    # need to index sample-wise reference
    bowtie-build ${vcf_gz.simpleName}.fasta ${vcf_gz.simpleName}
    """
}


/*
* Step 5. ViReMa
*/
process runVirema {
    cpus params.viremaCPU
    memory "$params.viremaMem GB"
    publishDir params.viremaPath

    input:
    file unalign from viremaChannel
    file index from viremaIndexChannel

    output:
    file "*.results"
    file "*Virus_Recombination_Results.txt" into viremaSum
    file "*tions.txt"
    file "*unaligned*.txt"
    file "*_rename.fastq"
    file "*_Virema_log.out" 
    file "DeDuped*.results" into resultsBam

    """
    awk '{print (NR%4 == 1) ? "@1_" ++i : \$0}' $unalign >  ${unalign.simpleName}_rename.fastq
    
    python ${params.viremaApp}/ViReMa.py \
        ${index[1].simpleName} \
        ${unalign.simpleName}_rename.fastq     \
        ${unalign.simpleName}.results \
        --MicroInDel_Length $params.micro \
        -DeDup \
        --Defuzz 3 \
        --Seed ${params.seed} \
        --N ${params.mismatch} \
        --X ${params.X} \
        --Output_Tag ${unalign.simpleName} \
        -ReadNamesEntry \
        --p $params.viremaCPU \
        > ${unalign.simpleName}_Virema_log.out
    """
}


/*
* Step 6. Compress ViReMa reads
*/
process runSummary {
    memory "2 GB"
    publishDir params.viremaPath

    input:
    file in_file from resultsBam

    output:
    file "*.bam"

    """
    samtools view -bS ${in_file} > ${in_file.baseName}.bam
    """
}

/*
* Step 6. ViReMa Summary of results (w/ perl scripts)
*/
process runSummary {
    memory "2 GB"
    publishDir params.viremaPath, mode: 'copy'

    input:
    file in_file from viremaSum

    output:
    file "*.par*"

    """
    perl $params.summaryApp -i $in_file -o ${in_file.baseName}.par -d 1
    perl $params.summaryApp -i $in_file -o ${in_file.baseName}.par5 -d 5
    """
}
