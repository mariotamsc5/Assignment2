#!/usr/bin/env nextflow

params.inputFile = "your_input_file.fasta" 
params.cutoff =  null 

// Define a channel for the input FASTA file
inputFileChannel = file(params.inputFile)

process checkGC {
    input:
    file fasta from inputFileChannel
    val cutoff from params.cutoff
    
    output:
    file 'output.txt'
    
    script:
    """
    #!/usr/bin/env python3
    # Calculate GC content for each sequence

    from Bio import SeqIO
    
    with open('${fasta}') as fasta_file, open('output.txt', 'w') as output_file:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            sequence = str(record.seq)
            gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence)
            
            if gc_content > cutoff:
                output_file.write(f'>{record.id}\\n{sequence}\\n')
    """
}

workflow {
    // Use Channel.fromPath to create a channel from a file path
    fastaRecord = Channel.fromPath(params.inputFile)
    // Call a process that iterates over the FASTA sequences
    // and writes sequences with GC content greater than the cutoff to output.txt
    checkGC(fastaRecord)
}
