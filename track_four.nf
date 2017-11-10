#!/usr/bin/env nextflow

def rev = workflow.revision ?: workflow.commitId ?: workflow.scriptId.substring(0,10)

log.info ''
log.info "================================================="
log.info " Bifrost assembly module version ${rev}"
log.info "================================================="
log.info "Reads                   : ${params.reads}"
log.info "#files in read set      : ${params.setsize}"
log.info "Results can be found in : ${params.out_dir}"
log.info "================================================="
log.info ""

preCmd = """
if [ -f /cluster/bin/jobsetup ];
then set +u; source /cluster/bin/jobsetup; set -u; fi
"""

// First, define the input data that go into input channels
Channel
    .fromFilePairs( params.reads, size:params.setsize )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into{ fastqc_reads, trimmomatic_reads }

// Second is to send all of through fastqc

process run_raw_fastqc {
    publishDir "${params.out_dir}/${params.fastqc}", mode: 'copy'

    tag { pair_id }

    input:
    set pair_id, file(reads) from fastqc_reads

    output:
    file "$pair_id" into raw_fastqc_eval

    """
    mkdir ${pair_id}
    $task.fastqc -q ${reads} -o ${pair_id} -t $task.threads
    """
}

process run_raw_fastqc_eval {
    publishDir "${params.out_dir}/${params.trimmed_fastqc_eval}", mode: 'copy'

    tag { pair_id }

    input:
    file "fastqc_output/*" from raw_fastqc_eval.toSortedList()

    output:
    file "raw_fastqc_eval.results"

    """
    fastqc_eval.py -d fastqc_output -o raw_fastqc_eval.results
    """
}

/*
 * Remove adapter sequences and low quality base pairs with Trimmomatic
 */
process run_trim {
	publishDir "${params.out_dir}/${params.trimmed}", mode: "copy"

	tag { pair_id }

    input:
    set pair_id, file(reads) from trimmomatic_reads

    set pair_id, file("${pair_id}_raw") from trimmomatic_reads

    output:
    set pair_id, file("${pair_id}_trimmed") into trimmed_fastqc_pairs

    """
    ${preCmd}
    mkdir ${pair_id}_trimmed
    $task.trimmomatic PE -threads $task.threads -trimlog ${pair_id}_trim.log ${pair_id}_raw/*${params.file_ending} \
        -baseout ${pair_id}_trimmed ILLUMINACLIP:$task.adapter_dir/${params.adapters}:{params.clipvalues} \
        LEADING:${params.leading} TRAILING:${params.trailing} \
        SLIDINGWINDOW:${params.slidingwindow} MINLEN:${params.minlen} &> ${pair_id}_run.log
    mv ${pair_id}_trimmed_1P ${pair_id}_trimmed/R1_trimmed${params.file_ending}
    mv ${pair_id}_trimmed_2P ${pair_id}_trimmed/R2_trimmed${params.file_ending}
    cat ${pair_id}_trimmed_1U ${pair_id}_trimmed_2U > ${pair_id}_trimmed/single${params.file_ending}
    """
}

process run_trim_fastqc {
    publishDir "${params.out_dir}/${params.trimmed_fastqc}", mode: 'copy'

    tag { pair_id }

    input:
    set pair_id, file("${pair_id}_trimmed") from trimmed_fastqc_pairs

    output:
    file "$pair_id" into trimmed_fastqc_results

    """
    mkdir ${pair_id}
    $task.fastqc -q ${reads} -o ${pair_id} -t $task.threads
    """
}


process run_trim_fastqc_eval {
    publishDir "${params.out_dir}/${params.trimmed_fastqc_eval}", mode: 'copy'

    tag { pair_id }

    input:
    file "fastqc_output/*" from fastqc_results.toSortedList()

    output:
    file "fastqc_eval.results"

    """
    fastqc_eval.py -d fastqc_output -o fastqc_eval.results
    """
}