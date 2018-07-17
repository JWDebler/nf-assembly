#!/home/johannes/bin nextflow
VERSION = "1.0.0"

//+++++++++++++++++ SETUP++++++++++++++++++++++++
//ID of isolate you want to assemble 
params.isolate = "P94-24" 

//Path to raw fastq.gz files
//This pipeline can be run either targeted at individual isolates
//or untargeted if you want to assemble and annotate everything
//within the input folder

//--> uncomment this for individual isolates
params.input = "/home/johannes/rdrive/PPG_SEQ_DATA-LICHTJ-SE00182/Basespace/**/LinaSample*/${params.isolate}*_R{1,2}_*fastq.gz"

// --> uncomment this for all isolates in input folder
//params.input = "/home/johannes/rdrive/PPG_SEQ_DATA-LICHTJ-SE00182/Basespace/**/LinaSample*/*_R{1,2}_*fastq.gz"

//Path where output data shall go
params.outputdir = "/home/johannes/rdrive/PPG_SEQ_DATA-LICHTJ-SE00182/johannes/notebook/test"
//+++++++++++++++++++++++++++++++++++++++++++++++

//Create channel that provides the sampleID and the raw read files 
reads = Channel.fromFilePairs(params.input).map {sampleID, fwdrevreads -> [sampleID.tokenize('_')[0], fwdrevreads]}


log.info "====================================================================="
log.info "Assembly, RNA annotation and Genemark prediction version " + VERSION
log.info "Isolate : ${params.isolate}"
log.info "Output  : ${params.outputdir}/${params.isolate}"
log.info "====================================================================="


reads
.groupTuple()
.map {sampleID, ary -> [sampleID, ary.transpose()]}
.map {sampleID, ary -> [sampleID, ary[0], ary[1]]}
.set {rawReads}

//Combine individual read files into one for forward and one for reverse reads
process combineReads {
    tag {sampleID}

    publishDir "${params.outputdir}/${sampleID}/01-combinedReads", mode: 'copy'

    input:
    set sampleID, "fwd.*.fastq.gz", "rev.*.fastq.gz" from rawReads

    output:
    set sampleID, "fwd.fastq.gz", "rev.fastq.gz" into rawCombinedReads

    """
    cat fwd.*.fastq.gz > fwd.fastq.gz
    cat rev.*.fastq.gz > rev.fastq.gz \
    """
}

//Trimmomatic
process trimReads {
    tag {sampleID}
    cpus 5
    publishDir "${params.outputdir}/${sampleID}/02-trimmedReads/", mode: 'copy'
        
    input:
    set sampleID, "fwd.fastq.gz", "rev.fastq.gz" from rawCombinedReads

    output:
    set sampleID, "paired.fwd.fastq.gz", "paired.rev.fastq.gz", "singletons.fastq.gz" into trimmedReads

    """
    cp /opt/trimmomatic/current/adapters/NexteraPE-PE.fa .

    java -jar /opt/trimmomatic/current/trimmomatic.jar PE \
    -threads ${task.cpus} \
    fwd.fastq.gz rev.fastq.gz \
    paired.fwd.fastq.gz singles.fwd.fastq.gz \
    paired.rev.fastq.gz singles.rev.fastq.gz \
    ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50

    cat singles.*.fastq.gz > singletons.fastq.gz
    """
}

//SPAdes assembly
process assemble {
    tag {sampleID}
    memory '25G'
    cpus 10

    input:
    set sampleID, "paired.fwd.fastq.gz", "paired.rev.fastq.gz", "singletons.fastq.gz" from trimmedReads

    output:
    set sampleID, "contigs.fasta" into contigsRaw
    set sampleID, "scaffolds.fasta" into cleanup

    """
    /opt/spades/current/bin/spades.py \
    -k 21,33,55,77,99,127 \
    --memory ${task.memory.toGiga()} \
    --threads ${task.cpus} \
    --careful \
    --pe1-1 paired.fwd.fastq.gz \
    --pe1-2 paired.rev.fastq.gz \
    --pe1-s singletons.fastq.gz \
    -o .
    """
}

process cleanupSpadesOutput {
    tag {sampleID}
    publishDir "${params.outputdir}/${sampleID}/03-assembly/spades/", mode: 'copy'

    input:
	set sampleID, "scaffolds.fasta" from cleanup

    output:
	set sampleID, "scaffolds.clean.fasta" into scaffoldsRawForGenemark
    set sampleID, "scaffolds.clean.fasta" into scaffoldsRawForTRNAscan
    set sampleID, "scaffolds.clean.fasta" into scaffoldsRawForInfernal

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
    out = File.open("scaffolds.clean.fasta", 'w')
    Bio::FlatFile
    .open("scaffolds.fasta")
    .sort_by{|entry| entry.length}
    .reverse
    .each
    .with_index(1)
    .map{|entry, i| [entry, "%sscf%d" % [options[:prefix], i]]}
    .each{|entry, name| out.puts entry.seq.to_fasta(name, 80)}
    out.close
    """
}

//GenemarkES annotation
process annotation_genemark {
    tag {sampleID}
    cpus 10
    publishDir "${params.outputdir}/${sampleID}/04-annotation/", mode: 'copy'

    input:
	set sampleID, "scaffolds.clean.fasta" from scaffoldsRawForGenemark

    output:
	set sampleID, "${sampleID}.genemark.gtf"

    """
    /opt/genemark-ES/gmes_petap.pl --ES --fungus --cores ${task.cpus} --sequence scaffolds.clean.fasta
    """
}

//tRNA annotation with tRNAscanSE
process annotation_trnascan {
    tag {sampleID}
    publishDir "${params.outputdir}/${sampleID}/04-annotation/", mode: 'copy'

    input:
	set sampleID, "scaffolds.clean.fasta" from scaffoldsRawForTRNAscan

    output:
	set sampleID, "${sampleID}.trnascanSE.gff3"
    """
    /home/johannes/local/bin/tRNAscan-SE -o trnascanoutput.out scaffolds.clean.fasta 
    /home/johannes/scripts/convert_tRNAScanSE_to_gff3.pl --input=trnascanoutput.out > trnascanSE.gff3
    """
}

//RNA annotation with infernal
process annotaton_infernal {
    tag {sampleID}
    cpus 10
    
    input:
	set sampleID, "scaffolds.clean.fasta" from scaffoldsRawForInfernal

    output:
	set sampleID, "scaffolds.cmscan.tbl" into infernalToGff3
    """
    wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/13.0/Rfam.clanin
    wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/13.0/Rfam.cm.gz
    gzip -d Rfam.cm.gz
    cmpress Rfam.cm
    cmscan --rfam --cpu ${task.cpus} --cut_ga --nohmmonly --tblout scaffolds.cmscan.tbl --fmt 2 --clanin Rfam.clanin Rfam.cm scaffolds.clean.fasta > infernal.cmscan
    """
}

//infernal output conversion to GFF3
process infernalToGff3 {
    tag {sampleID}
    publishDir "${params.outputdir}/${sampleID}/04-annotation/", mode: 'copy'

    input:
	set sampleID, "scaffolds.cmscan.tbl" from infernalToGff3

    output:
	set sampleID, "${sampleID}.infernal.gff3"
    """
    grep -v ^# scaffolds.cmscan.tbl > scaffolds.cmscan.clean.tbl && awk '{printf "%s\tinfernal\t%s\t%d\t%d\t%s\t%s\t.\tNote=RfamID-%s\\n" ,\$4,\$2,\$10,\$11,\$17,\$12,\$3}'  scaffolds.cmscan.clean.tbl > infernal.gff3
    """
}

workflow.onComplete {
    log.info "========================================================"
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'Failed' }"
    log.info "========================================================"
}
