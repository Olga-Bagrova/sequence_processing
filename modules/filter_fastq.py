def is_gc_enough(seq: str, gc_bounds: tuple)->bool:
    """
    Check if there is enough gc-content in the sequence
    arguments:
        - seq (str): sequence for checking
        - gc_bounds (tuple): lower and upper borders for gc-content of sequence
    return:
        - bool: True or False depending on the sufficiency of the gc content in the sequence
    """
    gc_content = 100*(seq.count('G') + seq.count('C'))/len(seq)
    if (gc_content >= gc_bounds[0]) and (gc_content <= gc_bounds[1]):
        return True
    else:
        return False


def is_length_enough(seq: str, length_bounds: tuple)->bool:
    """
    Check whether the length of the sequence satisfies the boundaries
    arguments:
        - seq (str): sequence for checking
        - length_bounds (tuple): lower and upper borders for length of sequence
    return:
        - bool: True or False depending on the sufficiency of the length of sequence
    """ 
    if (len(seq) >= length_bounds[0]) and (len(seq) <= length_bounds[1]):
        return True
    else:
        return False


def is_quality_enough(quality_string: str, quality_threshold: int)->bool:
    """
    Check whether the quality of the sequence satisfies the boundaries
    arguments:
        - seq (str): sequence for checking
        - quality_threshold (int): lower border for quality of sequence
    return:
        - bool: True or False depending on the sufficiency of the quality of sequence
    """
    quality = list(quality_string)
    q_score = []
    for symbol in quality:
        q_score.append(ord(symbol)-33)
    average_quality = sum(q_score) / len(q_score)
    if average_quality >= quality_threshold:
        return True
    else:
        return False


def filter_fastq(seqs: dict, gc_bounds: tuple = (0, 100), length_bounds: tuple = (0, 2**32), quality_threshold: int = 0)->dict:
    """
    Filter fastq-sequences based on its gc-content, length and quality.
    arguments:
        - seqs (dict): fastq-sequences for  (key (str) - name of the sequence, value (tuple) - sequence and its quality)
        - gc_bounds (tuple): lower and upper borders for gc-content of sequence
        - length_bounds (tuple): lower and upper borders for length of sequence
        - quality_threshold (int): lower border for quality of sequence
    return:
        - dict: fastq-sequences with True or False depending on the sufficiency of the quality of sequence
    """
    suitable_seqs = dict()
    if not isinstance(gc_bounds, tuple):
        gc_bounds = tuple((0, int(gc_bounds)))
    if not isinstance(length_bounds, tuple):
        length_bounds = tuple((0, int(length_bounds)))
    for key, value in seqs.items():
        if (is_gc_enough(value[0], gc_bounds)) and (is_length_enough(value[0], length_bounds)) and (is_quality_enough(value[1], quality_threshold)):
            suitable_seqs[key] = value    
    return suitable_seqs
