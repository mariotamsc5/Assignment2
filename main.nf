#!/usr/bin/env nextflow

// Define the workflow
workflow {
    // Input parameters
    input:
    path inputFile
    float cutoff

    // Create a channel for the input FASTA file
    input_fasta = channel.fromPath(inputFile)

    // Process to calculate GC content and filter sequences
    process checkGC {
        input:
        path fastaFile from input_fasta
        
        output:
        path 'output.txt'

        // Script to calculate GC content and filter sequences
        script:
        """
        # Calculate GC content for each sequence
        python3 << 'EOF'
import sys

def calculate_gc_content(seq):
    gc_count = sum(seq.count(base) for base in ['G', 'C', 'g', 'c'])
    return gc_count / len(seq)

cutoff = float(${cutoff})
output_file = open("output.txt", "w")

with open("${fastaFile}") as fasta:
    header = ""
    sequence = ""
    for line in fasta:
        if line.startswith('>'):
            if header:
                gc_content = calculate_gc_content(sequence)
                if gc_content > cutoff:
                    output_file.write(header.strip() + "\\n" + sequence.strip() + "\\n")
            header = line
            sequence = ""
        else:
            sequence += line.strip()

    # Process the last sequence in the file
    gc_content = calculate_gc_content(sequence)
    if gc_content > cutoff:
        output_file.write(header.strip() + "\\n" + sequence.strip() + "\\n")

output_file.close()
EOF
        """
    }
}
