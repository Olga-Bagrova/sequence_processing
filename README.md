# sequence_processing
> *This is the repo for the fifth homework of the BI Python 2023 course*

**sequence_processing** is a multifunctional tool for fastq-sequence, DNA-sequence, RNA-sequence and protein-sequence processing. It can be used to filter fastq-sequences, build complementary DNA and RNA chains, translate from a one-letter protein record to a three-letter one and others.

# Usage

To run the sequence_processing, first import it as module

`import sequence_processing`

Then use one or more its function: `filter_fastq`, `run_dna_rna_tools`, `run_protein_tools`.

## filter_fastq

**filter_fastq** is a function for filtering fastq-sequences. It is possible to filter sequences by GC-content, length and quality.

**Inputs**

* `seqs ` - dictionary of fastq sequences. The key is a string, the name of the sequence. The value is a tuple of two strings: sequence and quality. 
* `gc_bounds` - GC-content interval (in percent) for filtering (*default* `gc_bounds = (0, 100)`). If only one value is entered, the interval from 0 to the entered value is considered. Both borders are included.
* `length_bounds` - length interval for filtering (*default* `length_bounds = (0, 2**32)`). If only one value is entered, the interval from 0 to the entered value is considered. Both borders are included.
* `quality_threshold` - threshold value (phred33 scale) of the average quality of the read for filtering (*default* `quality_threshold = 0`. Reads with average quality for all nucleotides below the threshold are discarded.

**Outputs**

A dictionary with the original structure, but with sequences that satisfy the filtering conditions.

**Usage example**

```python
filter_fastq(fastq_dict, gc_bounds = (35, 80), length_bounds = (70, 88), quality_threshold = 32)
```


## run_dna_rna_tools

**run_dna_rna_tools** is a multifunctional algorithm for processing RNA and DNA sequences. It can be used to obtain transcribed, reversed and complementary (in two directions) sequences.

**Inputs**

The input of the program is string with DNA- or RNA-sequences. You can enter from one to several sequences of any register. The last argument must specify the selected sequence processing procedure (one of `'transcribe'`, `'reverse'`, `'complement'`, `'reverse_complement'`).

**Procedures** 

- `transcribe` — return the transcribed sequence
- `reverse` — return an inverted sequence
- `complement` — return a complementary sequence
- `reverse_complement` — return the reverse complementary sequence

**Outputs**

The output of the program is a string or a list of strings with processed sequences (case-preserving).

**Usage example**

```python
run_dna_rna_tools('ATG', 'transcribe') # 'AUG'
run_dna_rna_tools('ATG', 'aT', 'reverse') # ['GTA', 'Ta']
run_dna_rna_tools('AtG', 'complement') # 'TaC'
run_dna_rna_tools('ATg', 'reverse_complement') # 'cAT'
```


## run_protein_tools

Function `run_protein_tools.py` is for protein sequences. It can be used to obtain length of sequence, percentage of each amino acid in sequences, DNA-sequence by protein-sequence, sequence properties, coiled coil pattern positions and rename from a one-letter entry to a three-letter entry and vice versa.

**Inputs**

The input of the function is strings with protein-sequences. You can enter from one to several sequences of any register. The last argument must specify the selected sequence processing procedure (one of `'length'`, `'percentage'`, `'rename'`, `'to_dna'`, `'property'`, `'coiled_coil'`). You may enter the sequence either in one-letter form or in three-letter form.

**Procedures** 

- `length` — count the length of the sequence
- `percentage` —  list of dictionaries with the percentages of the corresponding amino acids in each sequence. The dictionary contains only amino acid residues whose percentage in the sequence is not equal to 0 (which are contained in the sequence at all). Also, the dictionary is ordered from the largest percentage of content to the smallest. Cases of amino acid residues are taken into account.
> :warning: Attention: We use rounding to 2 decimal places. In some cases, **the sum of percentages** of all amino acid residues for sequence **may not be exactly 100%** due to rounding.
- `rename` — return renamed sequence (from three-letter to one-letter if the initial is three-letter entry and vice versa if it is one-letter entry)
- `to_dna` — transform protein sequence(s) to DNA sequence(s). This uses only one variant of the triplet
- `property` - return string or list of string with properties: P or p - positive side chain (R, H, K), N or n - negative (D, E), L or l - polar (S, T, N, Q), S or s - special (C, G, P), H or h - hydrophobic (A, V, I, L, M, F, Y, W).
- `coiled_coil` - return list of positions in initial sequence, where coiled coil pattern (*hxxhcxc*: h - hydrophobic, c - charged) starts

**Outputs**

The output of the program is a string or a list of strings with processed sequences (case-preserving) or another data depends on chosen procedure.

**Usage example**

```python
run_protein_tools('AaafgTThgh', 'ALAARGASNASPVALALAalaALAHISGLYGLNGLULYSLEUlysMETPROSER', 'length') #[10, 18]
run_protein_tools('AaafgTThgh', 'ALAARGASNASPVALALAalaALAHISGLYGLNGLULYSLEUlysMETPROSER', 'percentage') #[{'a': 20.0, 'g': 20.0, 'T': 20.0, 'h': 20.0, 'A': 10.0, 'f': 10.0}, {'ALA': 16.67, 'ARG': 5.56, 'ASN': 5.56, 'ASP': 5.56, 'VAL': 5.56, 'ala': 5.56, 'HIS': 5.56, 'GLY': 5.56, 'GLN': 5.56, 'GLU': 5.56, 'LYS': 5.56, 'LEU': 5.56, 'lys': 5.56, 'MET': 5.56, 'PRO': 5.56, 'SER': 5.56}]
run_protein_tools('AaafgTThgh', 'ALAARGASNASPVALALAalaALAHISGLYGLNGLULYSLEUlysMETPROSER', 'rename') #['ALAalaalapheglyTHRTHRhisglyhis', 'ARNDVAaAHGQEKLkMPS']
run_protein_tools('AaafgTThgh', 'ALAARGASNASPVALALAalaALAHISGLYGLNGLULYSLEUlysMETPROSER', 'to_dna') #['GCAgcagcattcgggACCACCcatgggcat', 'GCACGAAATGATGTTGCAgcaGCACATGGGGAAGAGAAATTAaaaATGCCCTCG']
run_protein_tools('AaafgTThgh', 'ALAARGASNASPVALALAalaALAHISGLYGLNGLULYSLEUlysMETPROSER', 'property') #['HhhhsLLpsp', 'HPLNHHhHPSLNPHpHSL']
run_protein_tools('AaafgTThgh', 'ALAARGASNASPVALALAalaALAHISGLYGLNGLULYSLEUlysMETPROSER', 'AvcMDdErwgdgyqAvcMDdE', 'coiled_coil') #[None, None, [0, 14]]

```
This part of the tool was developed in the team of [Gleb Bobkov](https://github.com/GlebBobkov) and [Dmitry Matach](https://github.com/zmitserbio). See our work [here](https://github.com/GlebBobkov/HW4_Bobkov)

# Troubleshooting

It might be arised errors in the next cases:
* If you are not entering DNA or RNA sequences in `run_dna_rna_tools`
* If you are trying to transcribe a non-DNA sequence in  `run_dna_rna_tools`
* If you enter neither a one-letter nor a three-letter protein sequence in `run_protein_tools`



:star:*Thanks for reading!*:star:
