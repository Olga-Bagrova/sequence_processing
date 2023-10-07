DNA = {'A', 'T', 'G', 'C'}
RNA = {'A', 'U', 'G', 'C'}
TRANSCRIPTION = {
    'A': 'A',
    'T': 'U',
    'G': 'G',
    'C': 'C',
    'a': 'a',
    't': 'u',
    'g': 'g',
    'c': 'c'
}
COMPLEMETATION_RNA = {
    'A': 'U',
    'U': 'A',
    'G': 'C',
    'C': 'G',
    'a': 'u',
    'u': 'a',
    'g': 'c',
    'c': 'g'
}
COMPLEMETATION_DNA = {
    'A': 'T',
    'T': 'A',
    'G': 'C',
    'C': 'G',
    'a': 't',
    't': 'a',
    'g': 'c',
    'c': 'g'
}


def is_dna(seq: str)->bool:
    """
    Check whether the sequence is DNA 
    arguments:
        - seq (str): sequence for checking
    return:
        - bool: True or False depending on whether the sequence is DNA
    """
    seqset = set(seq.upper())
    return (seqset <= DNA)


def is_rna(seq: str)->bool:
    """
    Check whether the sequence is RNA 
    arguments:
        - seq (str): sequence for checking
    return:
        - bool: True or False depending on whether the sequence is RNA
    """
    seqset = set(seq.upper())
    return (seqset <= RNA)


def transcribe(seq: str)->str:
    """
    Transcribe the DNA-sequence
    arguments:
        - seq (str): DNA-sequence for transcription
    return:
        - str: transcribed DNA-sequence
    """
    res = ''
    if is_dna(seq):
        res = ''.join([TRANSCRIPTION[nucleotide] for nucleotide in seq])
    else:
        raise ValueError('Only DNA-sequence are able to be transcribed, please entry DNA-sequence')
    return res

