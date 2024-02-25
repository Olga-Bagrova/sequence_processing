import os

from Bio import SeqIO
from Bio import SeqUtils
from Bio import SeqRecord
from abc import ABC, abstractmethod


#-----------------4 обновленное---------------------------
def filter_fastq(input_path: str, output_filename: str = '', gc_bounds: tuple = (0, 100), length_bounds: tuple = (0, 2**32), quality_threshold: int = 0)->dict:
    """
    Filter fastq-sequences based on its gc-content, length and quality.
    arguments:
        - input_path (str): path to the fastq-file for reading
        - output_filename (str): path for saving file with result (filtered fastq)
        - gc_bounds (tuple): lower and upper borders for gc-content of sequence
        - length_bounds (tuple): lower and upper borders for length of sequence
        - quality_threshold (int): lower border for quality of sequence
    return:
        - dict of SeqRecords: filtered fastq-sequences
    """

    
    #seqs = SeqIO.parse(input_path, 'fastq')#seqs = read_fastq(input_path)
    if not isinstance(gc_bounds, tuple):
        gc_bounds = tuple((0, int(gc_bounds)))
    if not isinstance(length_bounds, tuple):
        length_bounds = tuple((0, int(length_bounds)))

    
    suitable_seqs = list(
        record
        for record in SeqIO.parse(input_path, "fastq")
        if ((gc_bounds[0] <= SeqUtils.gc_fraction(record.seq)*100 <= gc_bounds[1]) and 
            (length_bounds[0] <= len(record) <= length_bounds[1]) and 
            (sum(record.letter_annotations['phred_quality']) / len(record) >= quality_threshold))
    )  

    
    if output_filename == '':
        output_filename = 'good_quality.fastq'
    dir_name = 'fastq_filtrator_results'
    if not output_filename[-6:] == '.fastq':
        output_filename = output_filename + '.fastq'
    if not os.path.isdir(dir_name):
        os.mkdir(dir_name)
    output_path = os.path.join(dir_name, output_filename)
    if os.path.isfile(output_path):
        raise ValueError('File with such name is exist. Please, use another name for your output file.')
    SeqIO.write(suitable_seqs, output_path, 'fastq')
    #count = SeqIO.write(suitable_seqs, output_path, 'fastq')
    #print(count)
    print('The result of fastq-filtering is saved in', output_path)
    return suitable_seqs




#-----------------5-----------------------

class BiologicalSequence(ABC, str):
    def __init__(self, seq: str):
        self.seq = seq

    @abstractmethod
    def check_alphabet(self):
        return set(self.seq).issubset(self.alphabet)
        
class NucleicAcidSequence(BiologicalSequence):
    def __init__(self, seq: str):
        super().__init__(seq = seq)

    def complement(self):
        """
        Return the complemetary DNA or RNA sequence
        arguments:
            - seq (str): DNA or RNA sequence for transcription
        return:
            - str: complemetary DNA or RNA sequence
        """

        complementary_seq = ''.join([self.comlementation[nucleotide] for nucleotide in self])
        return type(self)(complementary_seq)


    def gc_content(self):
        return 100*(self.count('G') + self.count('C')) / len(self)#len(seq)

class DNASequence(NucleicAcidSequence):
    def __init__(self, seq: str):
        super().__init__(seq = seq)
        self.alphabet = {'A', 'T', 'G', 'C', 'a', 't', 'g', 'c'}
        self.comlementation = {
            'A': 'T',
            'T': 'A',
            'G': 'C',
            'C': 'G',
            'a': 't',
            't': 'a',
            'g': 'c',
            'c': 'g'
        }
        self.transcription = {
            'A': 'A',
            'T': 'U',
            'G': 'G',
            'C': 'C',
            'a': 'a',
            't': 'u',
            'g': 'g',
            'c': 'c'
        }
        
    def transcribe(self):
        """
        Transcribe the DNA-sequence
        arguments:
            - seq (str): DNA-sequence for transcription
        return:
            - str: transcribed DNA-sequence
        """
        transcribed_seq = ''.join([self.transcription[nucleotide] for nucleotide in self.seq])
        return type(self)(transcribed_seq)

class RNASequence(NucleicAcidSequence):
    def __init__(self, seq: str):
        super().__init__(seq = seq)
        self.alphabet = {'A', 'U', 'G', 'C', 'a', 'u', 'g', 'c'}
        self.comlementation = {
            'A': 'U',
            'U': 'A',
            'G': 'C',
            'C': 'G',
            'a': 'u',
            'u': 'a',
            'g': 'c',
            'c': 'g'
        }


class AminoAcidSequence(BiologicalSequence):
    def __init__(self, seq: str):
        super().__init__(seq = seq)
        self.alphabet = {
            'A', 'R', 'N', 'D', 'V', 
            'H', 'G', 'Q', 'E', 'I', 
            'L', 'K', 'M', 'P', 'S', 
            'Y', 'T', 'W', 'F', 'C', 
            'a', 'r', 'n', 'd', 'v', 
            'h', 'g', 'q', 'e', 'i', 
            'l', 'k', 'm', 'p', 's', 
            'y', 't', 'w', 'f', 'c'
        }
        
    def count_percentage(self):
        """
        Count percentage of each amino acid in sequence
        arguments:
            - seq (str): sequence for counting
        return:
            - dict: dictionary with counted percentage    
        """
        l = len(self)
        percentages = {}
        for i in range(0, l):
            aa = self[i:i+1]
            if aa not in percentages:
                percentages[aa] = 1
            else:
                percentages[aa] += 1
        percentages.update((key, round(value / l * 100, 2)) for key, value in percentages.items())
        percentages = {key: value for key, value in sorted(percentages.items(), key=lambda item: item[1], reverse=True)}
        return percentages

