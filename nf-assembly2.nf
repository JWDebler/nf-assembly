#!/home/johannes/bin nextflow
VERSION = "1.1.0"

//+++++++++++++++++ SETUP++++++++++++++++++++++++
// setup paths to programs 
params.trimmomatic = "/opt/trimmomatic/current/trimmomatic*.jar"
params.fastp = "fastp"
params.spades = "/opt/spades/current/bin/spades.py"
params.genemark = "/opt/genemark-ES/gmes_petap.pl"
params.infernal_cmpress = "cmpress"
params.infernal_cmscan = "cmscan"
params.trnascan = "tRNAscan-SE"
params.interproscan = "/opt/interproscan/current/interproscan.sh"

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
.into {rawReads; versions}

// Create a file that contains the version numbers of the tools run at the time
// the pipeline was run
process versions {
    tag {sampleID}
    publishDir "${params.outputdir}/${params.isolate}", mode: 'copy'

    input:
    set sampleID, "fwd.*.fastq.gz", "rev.*.fastq.gz" from versions

    output:
    file('versions.txt')

    """
    echo "Programs and versions used in this pipeline:" >> versions.txt
    date >> versions.txt
    echo "============================================" >> versions.txt
    echo "trimmomatic:" >> versions.txt
    java -jar ${params.trimmomatic} -version >> versions.txt
    echo "--------------------------------------------" >> versions.txt
    echo "FastP:" >> versions.txt
    ${params.fastp} -v 2>&1 | head >> versions.txt
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
    """
}

// Combine individual read files into one for forward and one for reverse reads

process combineReads {
    tag {sampleID}

    publishDir "${params.outputdir}/${sampleID}/00-rawReads/", mode: 'copy'

    input:
    set sampleID, "fwd.*.fastq.gz", "rev.*.fastq.gz" from rawReads

    output:
    set sampleID, "fwd.fastq.gz", "rev.fastq.gz" into rawCombinedReads

    """
    cat fwd.*.fastq.gz > fwd.fastq.gz
    cat rev.*.fastq.gz > rev.fastq.gz \
    """
}
/*
// Quality trimming with Trimmomatic
process trimReadsWithTrimmomatic {
    tag {sampleID}
    cpus 5

    publishDir "${params.outputdir}/${sampleID}/00-trimmedReads/trimmomatic/", mode: 'copy'
        
    input:
    set sampleID, "fwd.fastq.gz", "rev.fastq.gz" from rawCombinedReads

    output:
    set sampleID, "paired.fwd.fastq.gz", "paired.rev.fastq.gz", "singletons.fastq.gz" into trimmedReads

    """
    cp /opt/trimmomatic/current/adapters/NexteraPE-PE.fa .

    java -jar ${params.trimmomatic} PE \
    -threads ${task.cpus} \
    fwd.fastq.gz rev.fastq.gz \
    paired.fwd.fastq.gz singles.fwd.fastq.gz \
    paired.rev.fastq.gz singles.rev.fastq.gz \
    ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50

    cat singles.*.fastq.gz > singletons.fastq.gz
    """
}

// SPAdes short read assembly from Trimmomatic trimmed reads
// Since Trimmomatic only trims the reads, we want SPAdes to perfom
// some quality filtering. It does this by default
process assembleTrimmomatic {
    tag {sampleID}
    publishDir "${params.outputdir}/${sampleID}/01-assembly/spades/", mode: 'copy'
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

*/

quality = Channel.from(28)
rawCombinedReads
.combine(quality)
.set{trimInput}

// Quality trimming and error correction with FastP
// -l, --length_required: reads shorter than length_required will be discarded, default is 15. (int [=15])
// -e, --average_qual: if one read's average quality score <avg_qual, then this read/pair is discarded. Default 0 means no requirement (int [=0])
// -p, --overrepresentation_analysis: enable overrepresented sequence analysis.
// -c, --correction: enable base correction in overlapped regions (only for PE data), default is disabled
// -y, --low_complexity_filter: enable low complexity filter. The complexity is defined as the percentage of base that is different from its next base (base[i] != base[i+1]).
// -x, --trim_poly_x: enable polyX trimming in 3' ends.
process trimReadsWithFastP {
    tag {"${sampleID} avg quality ${quality}"}
    cpus 5

    publishDir "${params.outputdir}/${sampleID}/fastp${quality}/01-trimmedReads/", mode: 'copy'
        
    input:
    set sampleID, "fwd.fastq.gz", "rev.fastq.gz", quality from trimInput

    output:
    set sampleID, quality, "paired.fwd.fastq.gz", "paired.rev.fastq.gz", "singletons.fastq.gz" into trimmedReads
    set sampleID, "report.${sampleID}.html"

    """
    fastp \
    -i fwd.fastq.gz \
    -I rev.fastq.gz \
    -o paired.fwd.fastq.gz \
    -O paired.rev.fastq.gz \
    --unpaired1 singletons.fastq.gz \
    --unpaired2 singletons.fastq.gz \
    --detect_adapter_for_pe \
    -l 50 \
    -e ${quality} \
    -x \
    -y \
    -c \
    -p \
    -h "report.${sampleID}.html"
    """
}

// SPAdes short read assembly from fastP corrected files
// Since fastP perfomres quality trimming and correction, we run
// Spades with the '--only-assembler' flag
process assembleFastP {
    tag {"${sampleID} avg quality ${quality}"}
    publishDir "${params.outputdir}/${sampleID}/fastp${quality}/02-assembly/", mode: 'copy'
    memory '30G'
    cpus 10

    input:
    set sampleID, quality, "paired.fwd.fastq.gz", "paired.rev.fastq.gz", "singletons.fastq.gz" from trimmedReads

    output:
    set sampleID, quality, "${sampleID}.${quality}.contigs.fasta" into contigsForCleanup
    set sampleID, quality, "${sampleID}.${quality}.scaffolds.fasta" into scaffoldsForCleanup

    """
    ${params.spades} \
    -k 21,33,55,77 \
    --memory ${task.memory.toGiga()} \
    --threads ${task.cpus} \
    --careful \
    --only-assembler \
    -1 paired.fwd.fastq.gz \
    -2 paired.rev.fastq.gz \
    -s singletons.fastq.gz \
    -o .
    mv scaffolds.fasta ${sampleID}.${quality}.scaffolds.fasta
    mv contigs.fasta ${sampleID}.${quality}.contigs.fasta

    """
}

process cleanupSpadesOutputScaffolds {
    tag {"${sampleID} avg quality ${quality}"}
    
    input:
	set sampleID, quality, "${sampleID}.${quality}.scaffolds.fasta" from scaffoldsForCleanup

    output:
    set sampleID, quality, "${sampleID}.${quality}.scaffolds.clean.fasta" into scaffoldsForHeaderadjustment

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
    out = File.open("${sampleID}.${quality}.scaffolds.clean.fasta", 'w')
    Bio::FlatFile
    .open("${sampleID}.${quality}.scaffolds.fasta")
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
    tag {"${sampleID} avg quality ${quality}"}
    
    input:
	set sampleID, quality, "${sampleID}.${quality}.contigs.fasta" from contigsForCleanup

    output:
    set sampleID, quality, "${sampleID}.${quality}.contigs.clean.fasta" into contigsForHeaderadjustment

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
    out = File.open("${sampleID}.${quality}.contigs.clean.fasta", 'w')
    Bio::FlatFile
    .open("${sampleID}.${quality}.contigs.fasta")
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
    tag {"${sampleID} avg quality ${quality}"}

    publishDir "${params.outputdir}/${sampleID}/fastp${quality}/02-assembly/", mode: 'copy'

    input:
    set sampleID, quality, "${sampleID}.${quality}.scaffolds.clean.fasta" from scaffoldsForHeaderadjustment

    output:
    set sampleID, quality, "${sampleID}.${quality}.scaffolds.clean.fasta" into scaffoldsRawForGenemark
    set sampleID, quality, "${sampleID}.${quality}.scaffolds.clean.fasta" into scaffoldsRawForTRNAscan
    set sampleID, quality, "${sampleID}.${quality}.scaffolds.clean.fasta" into scaffoldsRawForInfernal

    """
    sed 's,>,&${sampleID}_,g' -i ${sampleID}.${quality}.scaffolds.clean.fasta
    """
}

process addSpeciesNameToFastaHeadersContigs {
    tag {"${sampleID} avg quality ${quality}"}

    publishDir "${params.outputdir}/${sampleID}/fastp${quality}/02-assembly/", mode: 'copy'

    input:
    set sampleID, quality, "${sampleID}.${quality}.contigs.clean.fasta" from contigsForHeaderadjustment

    output:
    set sampleID, quality, "${sampleID}.${quality}.contigs.clean.fasta" into contigsRawForGenemark
    set sampleID, quality, "${sampleID}.${quality}.contigs.clean.fasta" into contigsRawForTRNAscan
    set sampleID, quality, "${sampleID}.${quality}.contigs.clean.fasta" into contigsRawForInfernal

    """
    sed 's,>,&${sampleID}_,g' -i ${sampleID}.${quality}.contigs.clean.fasta
    """
}

// de novo gene annotation with GenemarkES
process annotation_genemark_scaffolds {
    tag {"${sampleID} avg quality ${quality}"}
    cpus 10
    publishDir "${params.outputdir}/${sampleID}/fastp${quality}/03-annotation/", mode: 'copy'

    input:
	set sampleID, quality, "${sampleID}.${quality}.scaffolds.genemark.fasta" from scaffoldsRawForGenemark

    output:
	set sampleID, quality, "${sampleID}.${quality}.scaffolds.genemark.gtf", "${sampleID}.${quality}.scaffolds.genemark.fasta" into proteinsFromGenemarkScaffolds

    """
    ${params.genemark} --ES --fungus --cores ${task.cpus} --sequence ${sampleID}.${quality}.scaffolds.genemark.fasta
    mv genemark.gtf ${sampleID}.${quality}.scaffolds.genemark.gtf
    """
}

process annotation_genemark_contigs {
    tag {"${sampleID} avg quality ${quality}"}
    cpus 10
    publishDir "${params.outputdir}/${sampleID}/fastp${quality}/03-annotation/", mode: 'copy'

    input:
	set sampleID, quality, "${sampleID}.${quality}.contigs.genemark.fasta" from contigsRawForGenemark

    output:
	set sampleID, quality, "${sampleID}.${quality}.contigs.genemark.gtf", "${sampleID}.${quality}.contigs.genemark.fasta"  into proteinsFromGenemarkContis

    """
    /opt/genemark-ES/gmes_petap.pl --ES --fungus --cores ${task.cpus} --sequence ${sampleID}.${quality}.contigs.genemark.fasta
    mv genemark.gtf ${sampleID}.${quality}.contigs.genemark.gtf
    """
}

process extractProteinsFromGenemarkContigs {
  tag {"${sampleID} avg quality ${quality}"}

  input:
  set sampleID, quality, "${sampleID}.genemark.gtf", "input.fasta" from proteinsFromGenemarkContis

  output:
  set sampleID, quality, "${sampleID}.genemark.proteins.fasta" into proteinsFromGenemarkContigs

  """
  /opt/genemark-ES/get_sequence_from_GTF.pl ${sampleID}.genemark.gtf input.fasta
  mv prot_seq.faa ${sampleID}.genemark.proteins.fasta
  """
}

process interproscan {
  tag {"${sampleID} avg quality ${quality}"}
  publishDir "${params.outputdir}/${sampleID}/fastp${quality}/03-annotation/", mode: 'copy'
  cpus 12

  input:
  set sampleID, quality, "proteins.fasta" from proteinsFromGenemarkContigs

  output:
  file "${sampleID}.interproscan.tsv"

  """
  ${params.interproscan} \
  --applications SignalP_EUK,Pfam,TMHMM,PANTHER,PRINTS,ProDom,ProSitePatterns,ProSiteProfiles,MobiDBLite\
  --cpu ${task.cpus} \
  --seqtype p \
  --disable-precalc \
  --goterms \
  --pathways \
  --iprlookup\
  --input proteins.fasta \
  --output-file-base ${sampleID}.interproscan \
  --format tsv
  """

}

// tRNA annotation with tRNAscanSE
process annotation_trnascan_scaffolds {
    tag {"${sampleID} avg quality ${quality}"}
    publishDir "${params.outputdir}/${sampleID}/fastp${quality}/03-annotation/", mode: 'copy'

    input:
	set sampleID, quality, "${sampleID}.${quality}.scaffolds.clean.fasta" from scaffoldsRawForTRNAscan

    output:
	set sampleID, quality, "${sampleID}.${quality}.scaffolds.trnascanSE.gff3"
    """
    ${params.trnascan} -o trnascanoutput.out ${sampleID}.${quality}.scaffolds.clean.fasta 
    ${params.scripts}/convert_tRNAScanSE_to_gff3.pl --input=trnascanoutput.out > ${sampleID}.${quality}.scaffolds.trnascanSE.gff3
    """
}

process annotation_trnascan_contigs {
    tag {"${sampleID} avg quality ${quality}"}
    publishDir "${params.outputdir}/${sampleID}/fastp${quality}/03-annotation/", mode: 'copy'

    input:
	set sampleID, quality, "${sampleID}.${quality}.contigs.clean.fasta" from contigsRawForTRNAscan

    output:
	set sampleID, quality, "${sampleID}.${quality}.contigs.trnascanSE.gff3"
    """
    ${params.trnascan} -o trnascanoutput.out ${sampleID}.${quality}.contigs.clean.fasta 
    ${params.scripts}/convert_tRNAScanSE_to_gff3.pl --input=trnascanoutput.out > ${sampleID}.${quality}.contigs.trnascanSE.gff3
    """
}

// RNA annotation with infernal
process annotaton_infernal_scaffolds {
    tag {"${sampleID} avg quality ${quality}"}
    cpus 10
    
    input:
	set sampleID, quality, "${sampleID}.${quality}.scaffolds.clean.fasta" from scaffoldsRawForInfernal

    output:
	set sampleID, quality, "scaffolds.cmscan.tbl" into infernalToGff3scaffolds
    """
    wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/14.1/Rfam.clanin
    wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/14.1/Rfam.cm.gz
    gzip -d Rfam.cm.gz
    ${params.infernal_cmpress} Rfam.cm
    ${params.infernal_cmscan} --rfam --cpu ${task.cpus} --cut_ga --nohmmonly --tblout scaffolds.cmscan.tbl --fmt 2 --clanin Rfam.clanin Rfam.cm ${sampleID}.${quality}.scaffolds.clean.fasta > infernal.cmscan
    """
}

process annotaton_infernal_contigs {
    tag {"${sampleID} avg quality ${quality}"}
    cpus 10
    
    input:
	set sampleID, quality, "${sampleID}.${quality}.contigs.clean.fasta" from contigsRawForInfernal

    output:
	set sampleID, quality, "contigs.cmscan.tbl" into infernalToGff3contigs
    """
    wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/13.0/Rfam.clanin
    wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/13.0/Rfam.cm.gz
    gzip -d Rfam.cm.gz
    ${params.infernal_cmpress} Rfam.cm
    ${params.infernal_cmscan} --rfam --cpu ${task.cpus} --cut_ga --nohmmonly --tblout contigs.cmscan.tbl --fmt 2 --clanin Rfam.clanin Rfam.cm ${sampleID}.${quality}.contigs.clean.fasta > infernal.cmscan
    """
}

// infernal output conversion to GFF3
process infernalToGff3scaffolds {
    tag {"${sampleID} avg quality ${quality}"}
    publishDir "${params.outputdir}/${sampleID}/fastp${quality}/03-annotation/", mode: 'copy'

    input:
	set sampleID, quality, "scaffolds.cmscan.tbl" from infernalToGff3scaffolds

    output:
	set sampleID, quality, "${sampleID}.${quality}.scaffolds.infernal.gff3"
    """
    grep -v ^# scaffolds.cmscan.tbl > scaffolds.cmscan.clean.tbl && awk '{printf "%s\tinfernal\t%s\t%d\t%d\t%s\t%s\t.\tNote=RfamID-%s\\n" ,\$4,\$2,\$10,\$11,\$17,\$12,\$3}'  scaffolds.cmscan.clean.tbl > ${sampleID}.${quality}.scaffolds.infernal.gff3
    """
}

process infernalToGff3contigs {
    tag {"${sampleID} avg quality ${quality}"}
    publishDir "${params.outputdir}/${sampleID}/fastp${quality}/03-annotation/", mode: 'copy'

    input:
	set sampleID, quality, "contigs.cmscan.tbl" from infernalToGff3contigs

    output:
	set sampleID, quality, "${sampleID}.${quality}.contigs.infernal.gff3"
    """
    grep -v ^# contigs.cmscan.tbl > contigs.cmscan.clean.tbl && awk '{printf "%s\tinfernal\t%s\t%d\t%d\t%s\t%s\t.\tNote=RfamID-%s\\n" ,\$4,\$2,\$10,\$11,\$17,\$12,\$3}'  contigs.cmscan.clean.tbl > ${sampleID}.${quality}.contigs.infernal.gff3
    """
}

workflow.onComplete {
    log.info "========================================================"
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'Failed' }"
    log.info "========================================================"
}
