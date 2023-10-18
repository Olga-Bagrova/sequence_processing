import modules.filter_fastq
import modules.run_dna_rna_tools
import modules.run_protein_tools
import os


def read_fastq(input_path)->dict:
    """
    Read data from fastq-file by the path and return a dictionary with sequences names as keys and lists of sequences and quality strings.
    arguments:
        - input_path (str): path to the fastq-file for reading
    return:
        - seqs (dict): fastq-sequences the key is the name of the sequence, the value is sequence and its quality
    """
    if not os.path.isfile(input_path):
        raise ValueError('File is not exist. Please, check up and try again!')
    seqs = dict()
    with open(input_path) as file:
        while True:
            name = file.readline().rstrip()
            if len(name) == 0:
                break
            sequence = file.readline().rstrip()
            file.readline().rstrip()
            quality_string = file.readline().rstrip()
            seqs[name] = (sequence, quality_string)                
    return seqs


def filter_fastq(input_path: str, output_filename: str, gc_bounds: tuple = (0, 100), length_bounds: tuple = (0, 2**32), quality_threshold: int = 0)->dict:
    """
    Filter fastq-sequences based on its gc-content, length and quality.
    arguments:
        - input_path (str): path to the fastq-file for reading
        - output_filename (str): path for saving file with result (filtered fastq)
        - gc_bounds (tuple): lower and upper borders for gc-content of sequence
        - length_bounds (tuple): lower and upper borders for length of sequence
        - quality_threshold (int): lower border for quality of sequence
    return:
        - dict: fastq-sequences with True or False depending on the sufficiency of the quality of sequence
    """
    seqs = read_fastq(input_path)
    suitable_seqs = dict()
    if not isinstance(gc_bounds, tuple):
        gc_bounds = tuple((0, int(gc_bounds)))
    if not isinstance(length_bounds, tuple):
        length_bounds = tuple((0, int(length_bounds)))
    for key, value in seqs.items():
        if (modules.filter_fastq.is_gc_enough(value[0], gc_bounds)) and (modules.filter_fastq.is_length_enough(value[0], length_bounds)) and (modules.filter_fastq.is_quality_enough(value[1], quality_threshold)):
            suitable_seqs[key] = value    
    return suitable_seqs


def run_dna_rna_tools(*args)->list:
    """
    Process DNA- and RNA-sequence(s)
    arguments:
        - args: sequence or sequences and operation for it or them
    operations:
        -'transcribe': transcribe the DNA-sequence(s)
        -'reverse': reverse the sequence(s)
        -'complement': return the complemetary sequence(s)
        -'reverse_complement': return the reversed complemetary sequence(s)
    return:
        - list (or str): processed sequence(s)
    """
    *seqs, procedure = args
    res = []
    functions = {
        'transcribe': modules.run_dna_rna_tools.transcribe,
        'reverse': modules.run_dna_rna_tools.reverse,
        'complement': modules.run_dna_rna_tools.complement,
        'reverse_complement': modules.run_dna_rna_tools.reverse_complement
        }
    for seq in seqs:
        if modules.run_dna_rna_tools.is_dna(seq) or modules.run_dna_rna_tools.is_rna(seq):
            after = functions[procedure](seq)
            res.append(after)
        else:
            raise ValueError('Incorrect input sequence(s), please try again')
    if len(res) == 1:
        res = res[0]
    return res


def run_protein_tools(*args)->list:
    *seqs, procedure = args
    res = []
    functions = {
        'length': modules.run_protein_tools.count_length,
        'percentage': modules.run_protein_tools.count_percentage,
        'rename': modules.run_protein_tools.rename_another_letter_entry,
        'to_dna': modules.run_protein_tools.transform_to_DNA_code,
        'property': modules.run_protein_tools.aa_property,
        'coiled_coil': modules.run_protein_tools.coiled_coil_find
        }
    for seq in seqs:
        if modules.run_protein_tools.is_seq_one_letter_protein(seq) or modules.run_protein_tools.is_seq_three_letter_protein(seq):
            after = functions[procedure](seq)
            res.append(after)
        else:
            raise ValueError('Incorrect input sequence(s), please try again')
    if len(res) == 1:
        res = res[0]
    return res
