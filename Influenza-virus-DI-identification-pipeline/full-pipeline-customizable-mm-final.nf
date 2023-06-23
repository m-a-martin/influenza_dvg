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
    set val(id), "${read1.baseName}.qualtrim.paired.fastq.gz", "${read2.baseName}.qualtrim.paired.fastq.gz" into kraken2Channel
    set val(id), "${read1.baseName}.qualtrim.paired.fastq.gz", "${read2.baseName}.qualtrim.paired.fastq.gz" into extractChannel
    file "*.qualtrim.unpaired.fastq.gz"
    stdout trim_out

    """
    trimmomatic PE \
    -threads $params.trimCPU -phred33 $read1 $read2 \
    ${read1.baseName}.qualtrim.paired.fastq.gz ${read1.baseName}.qualtrim.unpaired.fastq.gz \
    ${read2.baseName}.qualtrim.paired.fastq.gz ${read2.baseName}.qualtrim.unpaired.fastq.gz \
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
    set val(id), "${read1.baseName}.iav.fastq.gz", "${read2.baseName}.iav.fastq.gz" into fastqcChannel
    set val(id), "${read1.baseName}.iav.fastq.gz", "${read2.baseName}.iav.fastq.gz" into catChannel

    """
    kraken2 \
	--db $params.kraken2DB \
	--report ${read1.baseName}_report.txt \
	--gzip-compressed \
	--paired \
	$read1 \
	$read2 \
	> ${read1.baseName}_kraken.txt
    
    python3 ${params.kraken2ScriptPath}/extract_kraken_reads.py \
	-k ${read1.baseName}_kraken.txt \
	-s1 $read1 \
	-s2 $read2 \
	-o ${read1.baseName}.iav.fastq \
	-o2 ${read2.baseName}.iav.fastq \
	--fastq-output \
	--include-children \
	-t 11320 \
	-r ${read1.baseName}_report.txt \
	--max 1000000000 

    gzip ${read1.baseName}.iav.fastq
    gzip ${read2.baseName}.iav.fastq
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
    file  "*both.fq.gz" into bowtie2Channel

    """
    cat $read1 $read2 > ${pair_id}both.fq.gz
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
    file "*_unaligned.fq" into viremaChannel
    file "*.bam" into depthChannel

    """
    bowtie2 -p $params.bowtie2CPU -x $params.bowtie2_pad_index -U $in_cat --score-min $params.scoreMin \
    --un ${in_cat.baseName}_unaligned.fq | samtools view -bS - > ${in_cat.baseName}.bam
    """
}

/*
Step 5. Remove duplicates, count depth
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
    samtools sort -n -@ $params.bowtie2CPU -O sam ${in_bam} |  \
        samtools fixmate -@ $params.bowtie2CPU -O sam -m - - | \
        samtools sort -@ $params.bowtie2CPU -O sam | \
        samtools markdup -@ $params.bowtie2CPU -r -s -O sam - - | \
        samtools idxstats  - > ${in_bam.baseName}_dp.tsv
    """
}

/*
* Step 5. ViReMa
*/
process runVirema {
    cpus params.viremaCPU
    memory "$params.viremaMem GB"
    publishDir params.viremaPath, mode: 'copy'

    input:
    file unalign from viremaChannel

    output:
    file "*.results"
    file "*Virus_Recombination_Results.txt" into viremaSum
    file "*tions.txt"
    file "*unaligned*.txt"
    file "*_rename.fq"
    file "*_Virema_log.out" 

    """
    awk '{print (NR%4 == 1) ? "@1_" ++i : \$0}' $unalign >  ${unalign.baseName}_rename.fq
    
    python ${params.viremaApp}/ViReMa.py $params.virema_index ${unalign.baseName}_rename.fq ${unalign.baseName}.results \
    --MicroInDel_Length $params.micro -DeDup --Defuzz 3 --Seed ${params.seed} \
    --N ${params.mismatch} --X ${params.X} --Output_Tag $unalign.baseName -ReadNamesEntry --p $params.viremaCPU > ${unalign.baseName}_Virema_log.out
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
