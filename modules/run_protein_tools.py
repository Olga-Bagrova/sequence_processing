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


def count_percentage(seq: str)->dict:
    """
    Count percentage of each amino acid in sequence
    arguments:
        - seq (str): sequence for counting
    return:
        - dict: dictionary with counted percentage    
    """
    l = count_length(seq)
    percentages = {}
    if is_seq_one_letter_protein(seq):
        step = 1
        #print('step=', step)
    else:
        step = 3
        #print('step=', step)
    for i in range(0, len(seq), step):
        aa = seq[i:i+step]
        if aa not in percentages:
            percentages[aa] = 1
        else:
            percentages[aa] += 1
    percentages.update((key, round(value/l*100, 2)) for key, value in percentages.items())
    percentages={key: value for key, value in sorted(percentages.items(), key=lambda item: item[1], reverse=True)}
    return percentages


def rename_another_letter_entry (seq: str)->str:
    """
    Transform to three-letter amino acids entry if initial is one-letter, and to one-letter if initial is three-letter
    arguments:
        - seqs (str): sequence for transforming to one(three)-letter entire
    return:
        - str: transformed sequence
    """
    renamed_seq = ''
    if is_seq_one_letter_protein(seq):
        step = 1
        to_entry = TO_THREE_LETTER_DICT
    else:
        step = 3
        to_entry = TO_ONE_LETTER_DICT
    for i in range(0, len(seq), step):
        aa = seq[i:i+step]
        renamed_seq = renamed_seq + to_entry[aa]
    return renamed_seq


def transform_to_DNA_code(seq):
    """
    Transforming of an amino acid sequence/protein to DNA sequence
    arguments:
        - seqs (str): amino acid sequence of protein
    return:
        - str: sequence of protein in the DNA sequence form  
    """
    if is_seq_three_letter_protein(seq):
        seq = rename_another_letter_entry(seq)
    return ''.join([RETRANSLATE[aa] for aa in seq])


def aa_property (seq: str)->str:
    """
    Print properties for each amino acids in string
    arguments:
        - seq (str): sequence for property-transforming
    return:
        - str: sequence of properties: P or p - positive side chain (R, H, K), N or n - negative (D, E), L or l - polar (S, T, N, Q), S or s - special (C, G, P), H or h - hydrophobic (A, V, I, L, M, F, Y, W). 
    """
    if is_seq_three_letter_protein(seq):
        seq = rename_another_letter_entry(seq)
    property_string = ''
    for aa in seq:
        property_string = property_string + AA_PROPERTY[aa]
    return property_string  


def coiled_coil_find(seq: str)->list:
    """
    Find for the positions of the colied coil motif
    Coiled coils usually contain a repeating hxxhcxc pattern of hydrophobic (h) and charged (c) amino acid residues -  'hxxhcxc'
    arguments:
        - seq (str): sequence for coiled coil searching
    return:
        - list: list of positions with found coiled coil
    """
    finds = []
    pattern = 'hxxhcxc'
    seq = seq.lower()
    if is_seq_three_letter_protein(seq):
        seq = rename_another_letter_entry(seq)
    property_string = aa_property(seq)
    for pos in range (len(property_string)-len(pattern)+1):
        if (property_string[pos] == 'h') and (property_string[pos+3] == 'h') and ((property_string[pos+4] == 'p') or (property_string[pos+4] == 'n')) and ((property_string[pos+6] == 'p') or (property_string[pos+6] == 'n')):
            finds.append(pos)
    if len(finds) == 0:
        finds = None
    return finds


def run_protein_tools(*args)->list:
    *seqs, procedure = args
    res = []
    functions = {
        'length': count_length,
        'percentage': count_percentage,
        'rename': rename_another_letter_entry,
        'to_dna': transform_to_DNA_code,
        'property': aa_property,
        'coiled_coil': coiled_coil_find
        }
    for seq in seqs:
        if is_seq_one_letter_protein(seq) or is_seq_three_letter_protein(seq):
            after = functions[procedure](seq)
            res.append(after)
        else:
            raise ValueError('Incorrect input sequence(s), please try again')
    if len(res) == 1:
        res = res[0]
    return res
