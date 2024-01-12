#!/usr/bin/env nextflow

// Define a channel for the input FASTA file
inputFileChannel = file(params.inputFile)

process checkGC {
    input:
    file fasta from inputFileChannel
    
    output:
    file 'output.txt' into outputChannel
    
    script:
    """
    #!/usr/bin/env python3
    # Calculate GC content for each sequence

    import re
    from Bio import SeqIO
    
    cutoff = ${params.cutoff}
    
    with open('${fasta}') as fasta_file:
        list = []
        for record in SeqIO.parse(fasta_file, 'fasta'):
            sequence = str(record.seq)
            gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence)
            
            if gc_content > cutoff:
                list.append(f'>{record.id}\n{sequence}\n')
    """
}

// Define the workflow
workflow {
    fastaRecord = Channel.fromPath('data/bacterial_dna.fasta')
    // Call a process that iterates over the FASTA sequences
    // and writes sequences with GC content greater than the cutoff to output.txt

    
    // Define the output channel for the output.txt file
    outputChannel = result
}
