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


