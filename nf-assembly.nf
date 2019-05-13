#!/home/johannes/bin nextflow
VERSION = "1.2.0"

//+++++++++++++++++ SETUP++++++++++++++++++++++++
// setup paths to programs 
params.trimmomatic = "/opt/trimmomatic/current/trimmomatic*.jar"
params.fastp = "fastp"
params.spades = "/opt/spades/current/bin/spades.py"
params.genemark = "/opt/genemark-ES/gmes_petap.pl"
params.infernal_cmpress = "cmpress"
params.infernal_cmscan = "cmscan"
params.trnascan = "tRNAscan-SE"
params.quast = "/opt/quast/current/quast.py"

// point this to the scripts directory of this repository
params.scripts = "/home/johannes/rdrive/Johannes-DEBLEJ-SE00276/bioinformatics/nf-assembly/scripts"

// Path to raw fastq.gz files
// This pipeline can be run either targeted at individual isolates
// or untargeted if you want to assemble and annotate everything
// within the input folder

// --> uncomment this for individual isolates
// ID of isolate you want to assemble
 params.isolate = "P94-24" 
 params.input = "/home/johannes/rdrive/PPG_SEQ_DATA-LICHTJ-SE00182/Basespace/**/LinaSample*/${params.isolate}*_R{1,2}_*fastq.gz"

// --> uncomment this for all isolates in input folder
// params.input = "/home/johannes/rdrive/PPG_SEQ_DATA-LICHTJ-SE00182/Basespace/**/LinaSample*/*_R{1,2}_*fastq.gz"

// Path where output data shall go
params.outputdir = "/home/johannes/rdrive/PPG_SEQ_DATA-LICHTJ-SE00182/johannes/notebook/test"
//+++++++++++++++++++++++++++++++++++++++++++++++




// Create channel that provides the sampleID and the raw read files 
reads = Channel
.fromFilePairs(params.input)
.map {sampleID, fwdrevreads -> [sampleID.tokenize('_')[0], fwdrevreads]}

log.info "=============================================================================="
log.info "Illumina assembly, RNA annotation and Genemark prediction version " + VERSION
log.info "Isolate : ${params.isolate}"
log.info "Output  : ${params.outputdir}/${params.isolate}"
log.info "=============================================================================="

reads
.groupTuple()
.map {sampleID, ary -> [sampleID, ary.transpose()]}
.map {sampleID, ary -> [sampleID, ary[0], ary[1]]}
.set {rawReads}

// Create a file that contains the version numbers of the tools run at the time
// the pipeline was run
process versions {
    tag {sampleID}
    publishDir "${params.outputdir}/${params.isolate}", mode: 'copy'

    output:
    file('versions.txt')

    """
    echo "Programs and versions used in this pipeline:" >> versions.txt
    date >> versions.txt
    echo "============================================" >> versions.txt
    echo "trimmomatic:" >> versions.txt
    java -jar ${params.trimmomatic} -version >> versions.txt
    echo "--------------------------------------------" >> versions.txt
    echo "spades:" >> versions.txt
    ${params.spades} -v >> versions.txt
    echo "--------------------------------------------" >> versions.txt
    echo "GeneMarkES:" >> versions.txt
    ${params.genemark} | grep 'GeneMark-ES Suite version' >> versions.txt
    echo "--------------------------------------------" >> versions.txt
    echo "infernal:" >> versions.txt
    ${params.infernal_cmpress} -h | grep INFERNAL >> versions.txt
    echo "--------------------------------------------" >> versions.txt
    echo "trnaScanSE:" >> versions.txt
    ${params.trnascan} -h 2>&1 | head -2 | grep tRNAscan >> versions.txt
    echo "--------------------------------------------" >> versions.txt
    echo "Quast:" >> versions.txt
    ${params.quast} -v >> versions.txt
    """
}

// Combine individual read files into one for forward and one for reverse reads

process combineReads {
    tag {sampleID}

    publishDir "${params.outputdir}/${sampleID}/00-rawReads/", mode: 'copy'

    input:
    set sampleID, "fwd.*.fastq.gz", "rev.*.fastq.gz" from rawReads

    output:
    set sampleID, "${sampleID}.fwd.fastq.gz", "${sampleID}.rev.fastq.gz" into rawCombinedReads

    """
    cat fwd.*.fastq.gz > ${sampleID}.fwd.fastq.gz
    cat rev.*.fastq.gz > ${sampleID}.rev.fastq.gz \
    """
}

// Quality trimming with Trimmomatic
process trimReadsWithTrimmomatic {
    tag {sampleID}
    cpus 5

    publishDir "${params.outputdir}/${sampleID}/01-trimmedReads/trimmomatic/", mode: 'copy'
        
    input:
    set sampleID, "fwd.fastq.gz", "rev.fastq.gz" from rawCombinedReads

    output:
    set sampleID, "${sampleID}.trimmed.fwd.fastq.gz", "${sampleID}.trimmed.rev.fastq.gz", "${sampleID}.singletons.fastq.gz" into trimmedReads

    """
    cp /opt/trimmomatic/current/adapters/NexteraPE-PE.fa .

    java -jar ${params.trimmomatic} PE \
    -threads ${task.cpus} \
    fwd.fastq.gz rev.fastq.gz \
    ${sampleID}.trimmed.fwd.fastq.gz singles.fwd.fastq.gz \
    ${sampleID}.trimmed.rev.fastq.gz singles.rev.fastq.gz \
    ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50

    cat singles.*.fastq.gz > ${sampleID}.singletons.fastq.gz
    """
}

// SPAdes short read assembly from Trimmomatic trimmed reads
// Since Trimmomatic only trims the reads, we want SPAdes to perfom
// some quality filtering. It does this by default
process assemblySpades {
    tag {sampleID}
    
    memory '25G'
    cpus 10

    input:
    set sampleID, "paired.fwd.fastq.gz", "paired.rev.fastq.gz", "singletons.fastq.gz" from trimmedReads

    output:
    set sampleID, "${sampleID}.contigs.fasta" into contigsForCleanup
    set sampleID, "${sampleID}.scaffolds.fasta" into scaffoldsForCleanup

    """
    ${params.spades} \
    -k 21,33,55,77 \
    --memory ${task.memory.toGiga()} \
    --threads ${task.cpus} \
    --careful \
    -1 paired.fwd.fastq.gz \
    -2 paired.rev.fastq.gz \
    -s singletons.fastq.gz \
    -o .
    mv scaffolds.fasta ${sampleID}.scaffolds.fasta
    mv contigs.fasta ${sampleID}.contigs.fasta

    """
}

process cleanupSpadesOutputScaffolds {
    tag {sampleID}
    
    input:
	set sampleID, "${sampleID}.scaffolds.fasta" from scaffoldsForCleanup

    output:
    set sampleID, "${sampleID}.scaffolds.clean.fasta" into scaffoldsForHeaderadjustment

    """
    #!/usr/bin/env ruby

    require 'optparse'
    require 'bio'
    require 'pp'

    options = {:prefix => nil}
    OptionParser.new do |opts|
      opts.banner = "Usage: fix_names.rb [options] input.fasta"

      opts.on("-p", "--prefix [NAME]", "Strain prefix to prepend to contig name") do |p|
        options[:prefix] = p + "_"
      end
    end.parse!
    out = File.open("${sampleID}.scaffolds.clean.fasta", 'w')
    Bio::FlatFile
    .open("${sampleID}.scaffolds.fasta")
    .sort_by{|entry| entry.length}
    .reverse
    .each
    .with_index(1)
    .map{|entry, i| [entry, "%sscf%d" % [options[:prefix], i]]}
    .each{|entry, name| out.puts entry.seq.to_fasta(name, 80)}
    out.close
    """
}

process cleanupSpadesOutputContigs {
    tag {sampleID}
    
    input:
	set sampleID, "${sampleID}.contigs.fasta" from contigsForCleanup

    output:
    set sampleID, "${sampleID}.contigs.clean.fasta" into contigsForHeaderadjustment

    """
    #!/usr/bin/env ruby

    require 'optparse'
    require 'bio'
    require 'pp'

    options = {:prefix => nil}
    OptionParser.new do |opts|
      opts.banner = "Usage: fix_names.rb [options] input.fasta"

      opts.on("-p", "--prefix [NAME]", "Strain prefix to prepend to contig name") do |p|
        options[:prefix] = p + "_"
      end
    end.parse!
    out = File.open("${sampleID}.contigs.clean.fasta", 'w')
    Bio::FlatFile
    .open("${sampleID}.contigs.fasta")
    .sort_by{|entry| entry.length}
    .reverse
    .each
    .with_index(1)
    .map{|entry, i| [entry, "%sctg%d" % [options[:prefix], i]]}
    .each{|entry, name| out.puts entry.seq.to_fasta(name, 80)}
    out.close
    """
}

process addSpeciesNameToFastaHeadersScaffolds {
    tag {sampleID}
    publishDir "${params.outputdir}/${sampleID}/02-assembly/spades/", mode: 'copy'

    input:
    set sampleID, "${sampleID}.scaffolds.clean.fasta" from scaffoldsForHeaderadjustment

    output:
    set sampleID, "${sampleID}.scaffolds.clean.fasta" into scaffoldsRawForGenemark
    set sampleID, "${sampleID}.scaffolds.clean.fasta" into scaffoldsRawForTRNAscan
    set sampleID, "${sampleID}.scaffolds.clean.fasta" into scaffoldsRawForInfernal
    set sampleID, "${sampleID}.scaffolds.clean.fasta" into scaffoldsRawForQuast
    set sampleID, "${sampleID}.scaffolds.clean.fasta" into scaffoldsRawForFastQC

    """
    sed 's,>,&${sampleID}_,g' -i ${sampleID}.scaffolds.clean.fasta
    """
}

process addSpeciesNameToFastaHeadersContigs {
    tag {sampleID}
    publishDir "${params.outputdir}/${sampleID}/02-assembly/spades/", mode: 'copy'

    input:
    set sampleID, "${sampleID}.contigs.clean.fasta" from contigsForHeaderadjustment

    output:
    set sampleID, "${sampleID}.contigs.clean.fasta" into contigsRawForGenemark
    set sampleID, "${sampleID}.contigs.clean.fasta" into contigsRawForTRNAscan
    set sampleID, "${sampleID}.contigs.clean.fasta" into contigsRawForInfernal
    set sampleID, "${sampleID}.contigs.clean.fasta" into contigsRawForFastQC
    set sampleID, "${sampleID}.contigs.clean.fasta" into contigsRawForQuast

    """
    sed 's,>,&${sampleID}_,g' -i ${sampleID}.contigs.clean.fasta
    """
}

// de novo gene annotation with GenemarkES
process annotation_genemark_scaffolds {
    tag {sampleID}
    cpus 10
    publishDir "${params.outputdir}/${sampleID}/03-annotation/", mode: 'copy'

    input:
	set sampleID, "${sampleID}.scaffolds.clean.fasta" from scaffoldsRawForGenemark

    output:
	set sampleID, "${sampleID}.scaffolds.genemark.gtf"

    """
    ${params.genemark} --ES --fungus --cores ${task.cpus} --sequence ${sampleID}.scaffolds.clean.fasta
    mv genemark.gtf ${sampleID}.scaffolds.genemark.gtf
    """
}

process annotation_genemark_contigs {
    tag {sampleID}
    cpus 10
    publishDir "${params.outputdir}/${sampleID}/03-annotation/", mode: 'copy'

    input:
	set sampleID, "${sampleID}.contigs.clean.fasta" from contigsRawForGenemark

    output:
	set sampleID, "${sampleID}.contigs.genemark.gtf"

    """
    /opt/genemark-ES/gmes_petap.pl --ES --fungus --cores ${task.cpus} --sequence ${sampleID}.contigs.clean.fasta
    mv genemark.gtf ${sampleID}.contigs.genemark.gtf
    """
}

// tRNA annotation with tRNAscanSE
process annotation_trnascan_scaffolds {
    tag {sampleID}
    publishDir "${params.outputdir}/${sampleID}/03-annotation/", mode: 'copy'

    input:
	set sampleID, "${sampleID}.scaffolds.clean.fasta" from scaffoldsRawForTRNAscan

    output:
	set sampleID, "${sampleID}.scaffolds.trnascanSE.gff3"
    """
    ${params.trnascan} -o trnascanoutput.out ${sampleID}.scaffolds.clean.fasta 
    ${params.scripts}/convert_tRNAScanSE_to_gff3.pl --input=trnascanoutput.out > ${sampleID}.scaffolds.trnascanSE.gff3
    """
}

process annotation_trnascan_contigs {
    tag {sampleID}
    publishDir "${params.outputdir}/${sampleID}/03-annotation/", mode: 'copy'

    input:
	set sampleID, "${sampleID}.contigs.clean.fasta" from contigsRawForTRNAscan

    output:
	set sampleID, "${sampleID}.contigs.trnascanSE.gff3"
    """
    ${params.trnascan} -o trnascanoutput.out ${sampleID}.contigs.clean.fasta 
    ${params.scripts}/convert_tRNAScanSE_to_gff3.pl --input=trnascanoutput.out > ${sampleID}.contigs.trnascanSE.gff3
    """
}

// RNA annotation with infernal
process annotaton_infernal_scaffolds {
    tag {sampleID}
    cpus 10
    
    input:
	set sampleID, "${sampleID}.scaffolds.clean.fasta" from scaffoldsRawForInfernal

    output:
	set sampleID, "scaffolds.cmscan.tbl" into infernalToGff3scaffolds
    """
    wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/14.1/Rfam.clanin
    wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/14.1/Rfam.cm.gz
    gzip -d Rfam.cm.gz
    ${params.infernal_cmpress} Rfam.cm
    ${params.infernal_cmscan} --rfam --cpu ${task.cpus} --cut_ga --nohmmonly --tblout scaffolds.cmscan.tbl --fmt 2 --clanin Rfam.clanin Rfam.cm ${sampleID}.scaffolds.clean.fasta > infernal.cmscan
    """
}

process annotaton_infernal_contigs {
    tag {sampleID}
    cpus 10
    
    input:
	set sampleID, "${sampleID}.contigs.clean.fasta" from contigsRawForInfernal

    output:
	set sampleID, "contigs.cmscan.tbl" into infernalToGff3contigs
    """
    wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/13.0/Rfam.clanin
    wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/13.0/Rfam.cm.gz
    gzip -d Rfam.cm.gz
    ${params.infernal_cmpress} Rfam.cm
    ${params.infernal_cmscan} --rfam --cpu ${task.cpus} --cut_ga --nohmmonly --tblout contigs.cmscan.tbl --fmt 2 --clanin Rfam.clanin Rfam.cm ${sampleID}.contigs.clean.fasta > infernal.cmscan
    """
}

// infernal output conversion to GFF3
process infernalToGff3scaffolds {
    tag {sampleID}
    publishDir "${params.outputdir}/${sampleID}/03-annotation/", mode: 'copy'

    input:
	set sampleID, "scaffolds.cmscan.tbl" from infernalToGff3scaffolds

    output:
	set sampleID, "${sampleID}.scaffolds.infernal.gff3"
    """
    grep -v ^# scaffolds.cmscan.tbl > scaffolds.cmscan.clean.tbl && awk '{printf "%s\tinfernal\t%s\t%d\t%d\t%s\t%s\t.\tNote=RfamID-%s\\n" ,\$4,\$2,\$10,\$11,\$17,\$12,\$3}'  scaffolds.cmscan.clean.tbl > ${sampleID}.scaffolds.infernal.gff3
    """
}

process infernalToGff3contigs {
    tag {sampleID}
    publishDir "${params.outputdir}/${sampleID}/03-annotation/", mode: 'copy'

    input:
	set sampleID, "contigs.cmscan.tbl" from infernalToGff3contigs

    output:
	set sampleID, "${sampleID}.contigs.infernal.gff3"
    """
    grep -v ^# contigs.cmscan.tbl > contigs.cmscan.clean.tbl && awk '{printf "%s\tinfernal\t%s\t%d\t%d\t%s\t%s\t.\tNote=RfamID-%s\\n" ,\$4,\$2,\$10,\$11,\$17,\$12,\$3}'  contigs.cmscan.clean.tbl > ${sampleID}.contigs.infernal.gff3
    """
}

process quastScaffolds {
    tag {sampleID}
    publishDir "${params.outputdir}/${sampleID}/04-QC/scaffolds/", mode: 'copy'

    input: 
    set sampleID, "${sampleID}.fasta" from scaffoldsRawForQuast

    output:
    file("quast_results/*") into quastResultsScaffolds

    """
    ${params.quast} --fungus --conserved-genes-finding --min-contig 0 --threads 10 --split-scaffolds ${sampleID}.fasta
    """
}

process quastContigs {
    tag {sampleID}
    publishDir "${params.outputdir}/${sampleID}/04-QC/contigs/", mode: 'copy'

    input: 
    set sampleID, "${sampleID}.fasta" from contigsRawForQuast

    output:
    file("quast_results/*") into quastResultsContigs

    """
    ${params.quast} --fungus --conserved-genes-finding --min-contig 0 --threads 10 ${sampleID}.fasta
    """
}


workflow.onComplete {
    log.info "========================================================"
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'Failed' }"
    log.info "========================================================"
}
