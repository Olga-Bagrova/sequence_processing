import os


def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = ''):
    """
    Convert multiline fasta-file to oneline file.
    arguments:
        - input_fasta: input fasta-file for convertation
        - output_fasta: converted fasta-file
    """
    if output_fasta == '':
        output_fasta = 'oneline_' + input_fasta
    if output_fasta[-6:] != '.fasta':
        output_fasta = output_fasta + '.fasta'
    if not os.path.isfile(input_fasta):
        raise ValueError('Input file is not exist. Please, check up and try again!')
    if os.path.isfile(output_fasta):
        raise ValueError('File with such name is exist. Please, use another name for your output file.')
    with open(input_fasta) as input_file:
        with open(output_fasta, mode = 'w') as output_file:
            line = input_file.readline().rstrip()
            while True:
                seq = ''
                if len(line) == 0:
                    break
                seq = seq + line + '\n'
                line = input_file.readline().rstrip()
                while (not (line.startswith('>'))) and (len(line) > 0):
                    seq = seq + line
                    line = input_file.readline().rstrip()
                seq += '\n'
                output_file.write(seq)

convert_multiline_fasta_to_oneline('example_data/example_multiline_fasta.fasta.txt', 'uxotdput')
