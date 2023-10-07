ONE_LETTER_ENTRY = {
    'A', 'R', 'N', 'D', 'V', 
    'H', 'G', 'Q', 'E', 'I', 
    'L', 'K', 'M', 'P', 'S', 
    'Y', 'T', 'W', 'F', 'C', 
    'a', 'r', 'n', 'd', 'v', 
    'h', 'g', 'q', 'e', 'i', 
    'l', 'k', 'm', 'p', 's', 
    'y', 't', 'w', 'f', 'c'
}
THREE_LETTER_ENTRY = {
    'ALA', 'ARG', 'ASN', 'ASP', 'VAL', 
    'HIS', 'GLY', 'GLN', 'GLU', 'ILE', 
    'LEU', 'LYS', 'MET', 'PRO', 'SER', 
    'TYR', 'THR', 'TRP', 'PHE', 'CYS', 
    'ala', 'arg', 'asn', 'asp', 'val', 
    'his', 'gly', 'gln', 'glu', 'ile', 
    'leu', 'lys', 'met', 'pro', 'ser', 
    'tyr', 'thr', 'trp', 'phe', 'cys'
}
TO_THREE_LETTER_DICT = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'V': 'VAL', 
    'H': 'HIS', 'G': 'GLY', 'Q': 'GLN', 'E': 'GLU', 'I': 'ILE', 
    'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'P': 'PRO', 'S': 'SER', 
    'Y': 'TYR', 'T': 'THR', 'W': 'TRP', 'F': 'PHE', 'C': 'CYS', 
    'a': 'ala', 'r': 'arg', 'n': 'asn', 'd': 'asp', 'v': 'val', 
    'h': 'his', 'g': 'gly', 'q': 'gln', 'e': 'glu', 'i': 'ile', 
    'l': 'leu', 'k': 'lys', 'm': 'met', 'p': 'pro', 's': 'ser', 
    'y': 'tyr', 't': 'thr', 'w': 'trp', 'f': 'phe', 'c': 'cys'
}
TO_ONE_LETTER_DICT = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'VAL': 'V', 
    'HIS': 'H', 'GLY': 'G', 'GLN': 'Q', 'GLU': 'E', 'ILE': 'I', 
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PRO': 'P', 'SER': 'S', 
    'TYR': 'Y', 'THR': 'T', 'TRP': 'W', 'PHE': 'F', 'CYS': 'C', 
    'ala': 'a', 'arg': 'r', 'asn': 'n', 'asp': 'd', 'val': 'v', 
    'his': 'h', 'gly': 'g', 'gln': 'q', 'glu': 'e', 'ile': 'i', 
    'leu': 'l', 'lys': 'k', 'met': 'm', 'pro': 'p', 'ser': 's', 
    'tyr': 'y', 'thr': 't', 'trp': 'w', 'phe': 'f', 'cys': 'c',
}
RETRANSLATE = {
        'F': 'TTC', 'f': 'ttc',
        'L': 'TTA', 'l': 'tta',
        'S': 'TCG', 's': 'tcg',
        'Y': 'TAC', 'y': 'tac',
        'C': 'TGC', 'c': 'tgc',
        'W': 'TGG', 'w': 'tgg',
        'P': 'CCC', 'p': 'ccc',
        'H': 'CAT', 'h': 'cat',
        'Q': 'GAA', 'q': 'gaa',
        'R': 'CGA', 'r': 'cga',
        'I': 'ATT', 'i': 'att',
        'M': 'ATG', 'm': 'atg',
        'T': 'ACC', 't': 'acc',
        'N': 'AAT', 'n': 'aat',
        'K': 'AAA', 'k': 'aaa',
        'V': 'GTT', 'v': 'gtt',
        'A': 'GCA', 'a': 'gca',
        'D': 'GAT', 'd': 'gca',
        'E': 'GAG', 'e': 'gag',
        'G': 'GGG', 'g': 'ggg'
}
AA_PROPERTY = {
    'R':'P','H':'P','K':'P','D':'N','E':'N',
    'S':'L','T':'L','N':'L','Q':'L','C':'S',
    'G':'S','P':'S','A':'H','V':'H','I':'H',
    'L':'H','M':'H','F':'H','Y':'H','W':'H',
    'r':'p','h':'p','k':'p','d':'n','e':'n',
    's':'l','t':'l','n':'l','q':'l','c':'s',
    'g':'s','p':'s','a':'h','v':'h','i':'h',
    'l':'h','m':'h','f':'h','y':'h','w':'h'
}


def is_seq_three_letter_protein(seq: str)->bool:
    """
    Check whether the protein sequence is three-letter entry
    arguments:
        - seq (str): sequence for checking
    return:
        - bool: True or False depending on whether protein sequence is three-letter entry
    """
    seqset = {seq[i:i+3]  for i in range(0, len(seq), 3)}
    return (seqset <= THREE_LETTER_ENTRY)
    

def is_seq_one_letter_protein(seq: str)->bool:
    """
    Check whether the protein sequence is one-letter entry
    arguments:
        - seq (str): sequence for checking
    return:
        - bool: True or False depending on whether protein sequence is one-letter entry
    """
    seqset = set(seq.upper())
    return (seqset <= ONE_LETTER_ENTRY)


def count_length(seq: str)->int:
    """
    Count the length of the protein sequence
    arguments:
        - seq (str): sequence for length counting
    return:
        - int: length of protein sequence  
    """
    if is_seq_one_letter_protein(seq):
        return len(seq)
    else:
        return int(len(seq)/3)

