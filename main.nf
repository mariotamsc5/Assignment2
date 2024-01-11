#!/usr/bin/env nextflow

params.inputFile = file
params.cutoff = float

// Define a channel for the input FASTA file
inputFileChannel = file(params.inputFile).val

// Define the workflow
workflow {
    // Call a process that iterates over the FASTA sequences
    // and writes sequences with GC content greater than the cutoff to output.txt
    process checkGC {
        input:
        file fasta from inputFileChannel
        
        output:
        file 'output.txt' into outputChannel
        
        script:
        """
        # Calculate GC content for each sequence
        python3 << 'EOF'
import re
from Bio import SeqIO

cutoff = ${params.cutoff}

with open('${fasta}') as fasta_file, open('output.txt', 'w') as output_file:
    for record in SeqIO.parse(fasta_file, 'fasta'):
        sequence = str(record.seq)
        gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence)

        if gc_content > cutoff:
            output_file.write(f'>{record.id}\n{sequence}\n')
EOF
        """
    }
    
    // Define the output channel for the output.txt file
    outputChannel = result
}
    
    // Define the output channel for the output.txt file
    outputChannel = result
}
