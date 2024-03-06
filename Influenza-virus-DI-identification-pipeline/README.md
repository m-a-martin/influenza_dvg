# Influenza-virus-DI-identification-pipeline

# DESCRIPTION

Adapted from: [https://github.com/BROOKELAB/Influenza-virus-DI-identification-pipeline](https://github.com/BROOKELAB/Influenza-virus-DI-identification-pipeline).

This is an Illumina-based sequencing framework and bioinformatics pipeline capable of generating highly accurate and reproducible profiles of DIP-associated junction sequences originally published [here](https://doi.org/10.1128/jvi.00354-19) and adopted by Michael A. Martin, Nick Berg, and Katia Koelle. Included configuration files are for analysis of the influenza A H3N2 data published by McCrone et al. [here](https://doi.org/10.7554/eLife.35962). The pipeline is linux-based.


# DEPENDENCIES

This pipeline expects the following tools/languages to be installed as *linux modules* (http://modules.sourceforge.net/) and be available in your path:

- <b>Nextflow</b>    tested with v19.01.0.5050 ( download page https://github.com/nextflow-io/nextflow/releases )
- <b>Trimmomatic</b> tested with v0.38 ( download page https://github.com/timflutre/trimmomatic )
- <b>FastQC</b>      tested with v0.11.8 ( download page https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ )
- <b>Bowtie2</b>     tested with v2.3.4.3 ( download page  http://bowtie-bio.sourceforge.net/bowtie2/index.shtml )
- <b>Samtools</b>    tested with v1.9 ( download page http://www.htslib.org )
- <b>Bcftools</b>    tested with v1.8 ( download page http://www.htslib.org )
- <b>ViReMa          tested with modified v0.25 ( included in this repo)
- <b>Bowtie</b>      tested with version v1.2.2 ( download page http://bowtie-bio.sourceforge.net/index.shtml )
- <b>Perl</b>        tested with version v5.26.2  ( download page https://www.perl.org/ )

# WORKFLOW OF THE PIPELINE

In short, starting from a set of illumina paired sequencing reads the pipeline performs the following steps: 
1. Trim reads using Trimmomatic
2. Identify only influenza A reads using Kraken2 and the k2_pluspf_16gb database
3. Perform fastqc on the reads
4. Combine read pairs
5. Align reads to the reference genome in global mode
6. Count read depth, call consensus (>50%) variants, and generate run-specific consensus genome based on the reads that align in global mode
7. Run ViReMa on the unaligned reads using the run-specific consensus
8. Summarize results


# RUNNING THE PIPELINE

- This pipeline expects each sample of viral RNA to be made up of short Illumina reads (paired-end). We tested the pipeline with the data from [McCrone et al. eLife 2018](https://doi.org/10.7554/eLife.35962). 
- The sample(s) to be analyzed by this pipeline must be placed together in the same folder and must have a similar file naming patter, included config files expect .fastq.gz
- Configuration files provided in <b>new-config-files/</b>
- To run the pipeline type this command at the prompt: 

<pre>
nextflow -c config.file full-pipeline-customizable-mm-final.nf
</pre>

# OUTPUTS

Nextflow generates two folders to keep track of execution progress. You can delete them once the execution ends successfully. They are called <i>.nextflow/ </i> and <i>work/ </i>

The actual results of the pipeline are placed in these folders:

- <b>trimmomatic</b> contains the trimmed reads
- <b>fastqc</b>      contains the results of FastQC on  raw and trimmed reads
- <b>kraken2</b>	 contains the output of the kraken analysis to subset reads to just influenza A reads
- <b>bowtie2</b>	 contains the Bowtie2 output from the global alignment
- <b>depth</b>		 contains the read depth across the genome
- <b>var</b>		 contains the variants identified in each sequencing run
- <b>consensus</b>	 contains the run specific consensus genome including consensus variants and 210 nucleotide poly-A pads on either end
- <b>Virema</b>		 contains the results of running ViReMa on the unaligned reads to detect DIP-associated deletion junctions. Each sample will have several files with intermediary and final results. The final results are the files ending in <i> *.par </i>*.

# LICENSE

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>



