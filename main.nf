#!/usr/bin/env nextflow

params.out = "results"
params.accession = "M21012"
params.inputDir = "hepatitis/"
params.out_combined = "combined.fasta"
params.mafft_out = "aligned_out.fasta"
params.trimal_out = "trimmed_alignment.fasta"
params.trimal_report_html = "trimal_report.html"

process fetch_reference {
    conda 'bioconda::entrez-direct=24.0'

    input:
        val accession
    
    output:
        path "${accession}.fasta", emit: ref_fasta
    
    script:
        """
        esearch -db nucleotide -query "$accession" \
        | efetch -format fasta \
        > "${accession}.fasta"
        """
}

process combined_fasta {
    input:
        path inputDir_path
        path reference_fasta
    
    output:
        path "${params.out_combined}", emit: combined_out
    
    script:
    """
    cat "${inputDir_path}"/*.fasta "${reference_fasta}" > "${params.out_combined}"
    """
}

process mafft_alignment {
    conda 'bioconda::mafft=7.525'

    input:
        path input_fasta

    output:
        path "${params.mafft_out}", emit: aligned_fasta

    script:
    """
    mafft --auto --thread -1 "${input_fasta}" > "${params.mafft_out}"
    """
}

process trim_alignment {
    conda 'bioconda::trimal=1.5.0'

    publishDir "${params.out}/trimal_results", mode: 'copy', pattern: '*'

    input:
        path input_aligned_fasta

    output:
        path "${params.trimal_out}", emit: trimmed_fasta
        path "${params.trimal_report_html}", emit: trimal_html_report

    script:
    """
    trimal -in "${input_aligned_fasta}" -out "${params.trimal_out}" -automated1 -htmlout "${params.trimal_report_html}"
    """
}

workflow {
    fetch_reference(params.accession)

    def inputDir_channel = Channel
        .fromPath(params.inputDir, type: 'dir')

    combined_fasta(inputDir_channel, fetch_reference.out.ref_fasta)

    mafft_alignment(combined_fasta.out.combined_out)

    trim_alignment(mafft_alignment.out.aligned_fasta)
}
