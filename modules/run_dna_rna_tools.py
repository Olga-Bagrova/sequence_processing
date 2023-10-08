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


def reverse(seq: str)->str:
    """
    Reverse the DNA or RNA sequence
    arguments:
        - seq (str): DNA or RNA sequence for transcription
    return:
        - str: reversed sequence
    """
    return seq[::-1]


def complement(seq: str)->str:
    """
    Return the complemetary DNA or RNA sequence
    arguments:
        - seq (str): DNA or RNA sequence for transcription
    return:
        - str: complemetary DNA or RNA sequence
    """
    res = ''
    if is_rna(seq):
        res = ''.join([COMPLEMETATION_RNA[nucleotide] for nucleotide in seq])
    else:
        res = ''.join([COMPLEMETATION_DNA[nucleotide] for nucleotide in seq])
    return res


def reverse_complement(seq: str)->str:
    """
    Return the reversed complemetary DNA or RNA sequence
    arguments:
        - seq (str): DNA or RNA sequence for transcription and reversion
    return:
        - str: reversed complemetary DNA or RNA sequence
        """
    return reverse(complement(seq))


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
        'transcribe': transcribe,
        'reverse': reverse,
        'complement': complement,
        'reverse_complement': reverse_complement
        }
    for seq in seqs:
        if is_dna(seq) or is_rna(seq):
            after = functions[procedure](seq)
            res.append(after)
        else:
            raise ValueError('Incorrect input sequence(s), please try again')
    if len(res) == 1:
        res = res[0]
    return res
