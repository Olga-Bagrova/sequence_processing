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


def read_gbk_to_list(input_gbk: str)->list:
    """
    Return list of information about genes from gbk-file
    arguments:
        - input_gbk (str): path to the gbk-file for reading
    return:
        - list: information from gbk-format in list type
    """
    with open(input_gbk) as file:
        line = file.readline()
        while not line.startswith('FEATURES'):
            line = file.readline()
        file.readline()
        line = file.readline()
        info_start = ''
        for char in line:
            if not char.isalpha():
                info_start += char
            else:
                break
        genes = []
        while True:
            if not line:
                break
            gene_info = []
            while line.startswith(info_start):
                info_line = line.rsplit(info_start)[1]
                if info_line.endswith('"\n') or info_line[-2].isnumeric():
                    gene_info.append(info_line.rstrip())
                else:
                    line = file.readline()
                    while not (line.endswith('"\n') or line[-1].isnumeric()):
                        info_line = info_line.rstrip()
                        info_line = info_line + line.rsplit(info_start[:-1])[1]
                        line = file.readline()
                        if not line:
                            break
                    gene_info.append(info_line.rstrip())
                line = file.readline()
            line = file.readline()
            genes.append(gene_info)
    return genes


def extract_genes(gbk_genes_info: list, inp_genes: list, n_before: int = 1, n_after: int = 1)->list:
    """
    Extract genes from gbk-list 
    arguments:
        - gbk_genes_info (list): information from gbk-format in type of list
        - inp_genes (list): gene names, which are interested
        - n_before (int): additional genes before the interesting one
        - n_after (int): additional genes after the interesting one
    return:
        - list: extracted information from gbk-format about interesting and additional genes
    """
    genes_to_file = []
    if not isinstance(inp_genes, list):
        genes = [inp_genes]
    for gene in genes:
        gene = 'gene="' + gene + '"'
        matching = [gbk_genes_info.index(gene_info) for gene_info in gbk_genes_info if gene in gene_info]
    for matched_index in matching:
        if matched_index-n_before <0:
            start = 0
        else: start = matched_index-n_before
        if (matched_index+n_after+1) > len(gbk_genes_info):
            end = len(gbk_genes_info)
        else:
            end = matched_index+n_after+1
        for i in range (start, end):
            genes_to_file.append(gbk_genes_info[i])
    return genes_to_file


def write_to_fasta(extracted_genes: list, output_fasta: str):
    """
    Write genes sequences in gbk-format to fasta-file
    arguments:
        - extracted_genes (list): information about genes in gbk-format in type of list
        - output_fasta (str): name of fasta-file
    """
    with open(output_fasta, mode = 'w') as file:
        for gene_info in extracted_genes:
            name_string = '>'
            for info in gene_info:
                if not info.startswith('translation='):
                    name_string = name_string + info + '|'
                else:
                    file.write(name_string[:-1]+'\n')
                    file.write(info.split('translation="')[1]+'\n')


def select_genes_from_gbk_to_fasta(input_gbk, genes, n_before = 1, n_after = 1, output_fasta = ''):
    """
    Extract the chosen genes from gbk-file and write them to fasta-file
    arguments:
        - input_gbk (str): path to the gbk-file for processing
        - genes (list): genes names, which are interested
        - n_before (int): additional genes before the interesting one
        - n_after (int): additional genes after the interesting one
        - output_fasta (str): name of fasta-file for saving the result
    """
    gbk_genes_info = read_gbk_to_list(input_gbk)
    extracted_genes = extract_genes(gbk_genes_info, genes, n_before=1, n_after=1)
    if output_fasta == '':
        output_fasta = 'out_' + input_gbk
    output_fasta += '.fasta'
    write_to_fasta(extracted_genes, output_fasta)
