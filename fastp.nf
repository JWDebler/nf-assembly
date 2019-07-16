#!/home/johannes/bin nextflow
VERSION = "1.1.0"

//+++++++++++++++++ SETUP++++++++++++++++++++++++
// setup paths to programs 
params.trimmomatic = "/opt/trimmomatic/current/trimmomatic*.jar"
params.fastp = "fastp"
params.fastqc = "/opt/fastQC/current/fastqc"

// --> uncomment this for individual isolates
// ID of isolate you want to assemble
 params.isolate = "P94-24" 
 params.input = "/home/johannes/rdrive/PPG_SEQ_DATA-LICHTJ-SE00182/Basespace/**/LinaSample*/${params.isolate}*_R{1,2}_*fastq.gz"

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
.into {rawReads}

quality = Channel.from(21,22,23,24,25,26,27,28,29,30)

// Combine individual read files into one for forward and one for reverse reads

process combineReads {
    tag {sampleID}

    input:
    set sampleID, "fwd.*.fastq.gz", "rev.*.fastq.gz" from rawReads

    output:
    set sampleID, "fwd.fastq.gz", "rev.fastq.gz" into rawCombinedReads

    """
    cat fwd.*.fastq.gz > fwd.fastq.gz
    cat rev.*.fastq.gz > rev.fastq.gz \
    """
}

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

    publishDir "${params.outputdir}/${sampleID}/comparisons", mode: 'copy'
        
    input:
    set sampleID, "fwd.fastq.gz", "rev.fastq.gz", quality from trimInput

    output:
    set sampleID, quality, "${quality}.paired.fwd.fastq.gz", "${quality}.paired.rev.fastq.gz", "${quality}.singletons.fastq.gz" into trimmedReads
    set sampleID, "report.${sampleID}.${quality}.html"

    """
    fastp \
    -i fwd.fastq.gz \
    -I rev.fastq.gz \
    -o ${quality}.paired.fwd.fastq.gz \
    -O ${quality}.paired.rev.fastq.gz \
    --unpaired1 ${quality}.singletons.fastq.gz \
    --unpaired2 ${quality}.singletons.fastq.gz \
    --detect_adapter_for_pe \
    -l 50 \
    -e ${quality} \
    -x \
    -y \
    -c \
    -p \
    -h "report.${sampleID}.${quality}.html"
    """
}

process fastQC {
    tag {"${sampleID} avg quality ${quality}"}

    publishDir "${params.outputdir}/${sampleID}/comparisons/fastQC", mode: 'copy'

    input:
    set sampleID, quality, "${quality}.paired.fwd.fastq.gz", "${quality}.paired.rev.fastq.gz", "${quality}.singletons.fastq.gz" from trimmedReads

    output:
    file('*fastqc*')

    """
    ${params.fastqc} *.gz
    """

}

workflow.onComplete {
    log.info "========================================================"
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'Failed' }"
    log.info "========================================================"
}
