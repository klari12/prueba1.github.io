## Python: Data structures part1

Lists are formed of elements that we can modify. Strings, on the contrary, do not allow modification of any of their parts.

Lists are also objects. Let's imagine we make a list with the expression of a gene but we want to change one of its elements because you have made a mistake:



```python
gene_expression = ["SHOX1", 2.234, 500]
print(gene_expression)
del gene_expression[0]
print(gene_expression)
gene_expression.extend(["TRA11"])
gene_expression.reverse()
print(gene_expression)
```

    ['SHOX1', 2.234, 500]
    [2.234, 500]
    ['TRA11', 500, 2.234]
    

## Python: Data structures part 2

Tuple: a number of values separated by commas, another standard sequence data type, like strings and lists. The values of the tuple cannot be changed. Usually, they are used to contain heterogeneus elements.


```python
t=1,2,3
t
```

A set is an unordered collection with no duplicate elements. Set objects support mathematical operations like union, intersection and difference.


```python
brca1={"DNA repair", "zinc ion binding", "cancer"}
brca2={"nucleoplasma", "double-strand break repair", "cancer"}
# Union operation
brca1 | brca2
# intersection (to see what they have in common)
brca1 & brca2

```




    {'cancer'}



A dictionary is an unordered set of key and value pairs, with the requirement that the keys are unique (within one dictionary)


```python
TF_motif={"SP1":"gggccgg", "c-Myc":"atgatg", "Oct-1":"atgcaaat"}
TF_motif["c-Myc"]
# If we want to add a key:valure to the dictionary:
TF_motif["AP1"] = "gctagcg"
TF_motif
```




    {'SP1': 'gggccgg', 'c-Myc': 'atgatg', 'Oct-1': 'atgcaaat', 'AP1': 'gctagcg'}



#### IFS

We can create an "if" statement like this:


```python
DNA = input("Insert your DNA sequence ")

```

    Insert your DNA sequence acgttcgacgtNNN
    


```python
if "n" in DNA :
    nbases = DNA.count("n")
    print("This DNA sequence has \n%d undefined bases" %nbases)
elif "N" in DNA: 
    Nbases = DNA.count("N")
    print("This DNA sequence has \n%d undefined bases" %Nbases)
else:
    print("This sequences has no undefined bases")
```

    This DNA sequence has 
    5 undefined bases
    

This statements are based on boolean expressions, which are formed with the help of comparison (==, !=, >, etc), identity (is and is not) and membership (in and not in) operators.


```python
if "n" in DNA or "N" in DNA :
    nbases = DNA.count("n") + DNA.count("N")
    print("This DNA sequence has \n%d undefined bases" %nbases)
else:
    print("This sequences has no undefined bases")
```

    This DNA sequence has 
    4 undefined bases
    

#### LOOPS

Very important for analysing several sequences in an automatic way.


```python
dna = input("Insert your DNA sequence ")
```

    Insert your DNA sequence acgatygctagcatcgatctcacgatgactagcctca
    


```python
pos = dna.find("gt", 0)
pos2 = dna.find("ctca", 0)
```


```python
while pos>-1: 
    print("Donor splice site candidate at position %d"%pos)
    pos = dna.find("gt", pos+1)
```

    Donor splice site candidate at position 2
    


```python
len(dna)
```




    31




```python
while pos2>-1: 
    pos2list = pos2
    print("Braching points detected at position: %d"%pos2list)
    pos2 = dna.find("ctca", pos2+1)

```

    Braching points are at positions: 18
    Braching points are at positions: 33
    

The "continue" statement causes the program to continue with the next iteration of the nearest enclosing loop, skipping the rest of the code in the loop.
For example, it can be used to delete all invalid aminoacid characters from a protein sequence.
The continue statement can improve the readibility of the code.


```python
protein = "SCGSFDJJCBSTA"
corrected_protein = ""
for i in range(len(protein)):
    if protein[i] not in "ABCDEFGHIKLMNPQRSTVWXYZ":
        continue
    corrected_protein=corrected_protein+protein[i]
print("Corrected protein sequence: %s"%corrected_protein)
```

    Corrected protein sequence: SCGSFDCBSTA
    

Python allows the use of else in loops.  

The "pass" statement is a placeholder: it does nothing.

#### FUNCTIONS

Write a program that checks if a given DNA sequence contains an in frame stop codon


```python
 dna = input("Enter your sequence:")
```

    Enter your sequence:acgagctagcatgcacgatataaacgatcga
    


```python

def has_stop_codon(dna,frame=0):
    "If we write frame=0, then 0 will be its default value, but we can change it when calling the function"
    "The function detects stop condons in the specified reading frame"
    stop_codon_found= False
    stop_codons=["tga", "tag", "taa"]
    for i in range (frame, len(dna), 3):
        "We use the .lower method to transform all letter into lowercase to avoid mistakes"
        codon=dna[i:i+3].lower()
        if codon in stop_codons :
            stop_codon_found=True
            break
    return(stop_codon_found)

```


```python
has_stop_codon("agctaataatatataatacgtgagatcagctacgtacg", 3)
```




    True



#### Reading and writing files

Exercise: create a dictionary containing the sequences from a FASTA file


```python
FASTA_file = open('A.fumigatus.genome.fna')
# This function is for counting how many names of sequences are in the FASTA file
sequencenames = 0
for line in FASTA_file:
    if line[0] == ">":
        sequencenames += 1
print(sequencenames)
```

    228
    


```python
FASTA_file = open('A.fumigatus.genome.fna')
sequence_dict = {}
for line in FASTA_file:
    line= line.rstrip()
    if line[0]==">":
        words=line.split()
        name=words[0][1:]
        sequence_dict[name]=''
    else:
        sequence_dict[name] = sequence_dict[name] + line
FASTA_file.close()
```

There are multiple possible ways to do this. I did the way below but it gives que whole name for the sequence, instead of only the fist part as when using the code chunk from above.


```python
FASTA_file = open('A.fumigatus.genome.fna')
sequence_dictionary = {}
count=0
Lines = FASTA_file.readlines()
for line in Lines:
    if line[0] == ">":
        key = line.strip()
        sequence_dictionary[key] = ''
    else:
        sequence_dictionary[key] = sequence_dictionary[key]+line
        

```


```python
list(sequence_dictionary)
sequence_dictionary['>MCQI02000050.1 Aspergillus fumigatus strain LMB-35Aa scaffold103_size1834, whole genome shotgun sequence']
```




    'AAACCTCTAAACCCTATAGgaaaatagaggaaggtttatactattataataataaatagaaGCTAGTAGGTAAAAACTTT\nATAGAGAAAATACTATTTTTTAATAGGGCTACTCTGTAGCctatctaccttactctctatatatatttaatatatttata\ngtattacttaggttatcttctagttcctataagttcttatggggaatcctttctacttaactaggtataatatcttttct\nttttaaGTCTATTTAGCTAGTATATTCTTTAATTATTCCTCTATTAATATCTTTAAGGGACTAGGGTTCTAATAGtatct\ntattattatcctAAGGGTTTTACTTCTtaagaagggatatataaAAGATAAGATaaatttatctatataattttagtatA\nGAGTATATACCTACTTTctaataatagtagtaaTTTCAAACcgcctattatatttattattaagtttcttatatagGCAC\nTATAGTTTTAAGTTTCTTAAGAAGAGTAGCACCTTTTCTCTAATAAAAAATAAGATATTTCTAGTACTTTTAGTTATAAT\nACCTAGCATATATCTTCTAAGTAGTTACTAGGTGACCTTATAATTAATTTTATAGATCCTAGAGTATTTTAAGGCATTGC\ntctactactattactaataaattattagtaatattagtaagctaGAATAATACTAGAATATAGCTATAAGTAGCTTTAAA\nGGGTATTATTTTTATCAAGGTATAtaggatatagctatagatagTAAAGTTTATCTATAACTATATggcctaattatctt\nactaataattaataaatatTCTTAGGTATATCTCTAATATCTAATTCTAGTattctatctagctatttatttataaataa\naaaGCTATACTTATCTATTTTTTAATTCTTAGATATCTATAGAAATTAGTCTAGTAatagctagtaaatagtatatccct\nattagttactaggcttcttaatatactatattaaaGGATAATTCTCTCTAAGAATATATTAGCAAGTTCCTCTGTAGTAA\nTTATCTTcttaataaagatatatttagctatttttatatattaattaataatcACTAgaattatattaaatatcttTTAT\nTTTTATACTAACAGGGGCAGTAAggtaataaaatctataaataAGTTACTCTATAGCTTTTTTAGTATAGGTAGTGTAGT\nAACCTTCCCTATAGTAAGATACTAAGAAGTCTTAGTCTTctaatatactatatatatagctatattctttaatattatat\ntatatcccCTTCTATTAGAACttatactagagtaaggctaCTATCTTATTAACTCCTAGATATCTAGTAtaggggttatt\natagtatattcTTAGAAGCTCTTATCttattataatattagtaaaTATATAAATAGTATTCTTAAATTATAGCAGACCCT\nTACTATTAATTTTCTAGTTTCTTATCTAAGTATTACTAGAGATCCCCTATAGCTACTTagaaactagttatatataggta\ntcCCTATGCTATAgttcctagataatattctttaAGTAAGTTATTAGTAAGTCTTAAGTAATCTTAAGGCTAGTAGCctc\ntattactagtaatCTTAATAATAGGAATTCTAAACTTCCTACTCTTATTAGcatactataaatactaaggttatttatac\ntagtctagtttaactataagccTAGTATTAGCTTCCTCTAAAGTATAGGTAAGAGTCTAGCTCTAGTTTCTATATTCTTA\nATTAGCCTATAATCTTAGTATCTTAatagtatatctataggGTTAAGCTAtcctagtttatatttaattataaa\n'



The following function allows us to detect ORFs in a DNA sequence:


```python
def startstop_codon(dna, frame):
    dna = dna.upper()
    for i in range(frame, len(dna), 3):
        codon1 = dna[i:i+3]
        if codon1 == 'ATG':
            position1 = i
            for j in range(position1, len(dna), 3):
                codon2 = dna[j:j+3]
                if codon2 in ['TAA', 'TAG', 'TGA']:
                    position2 = j
                    yield (position2-position1+3, dna[position1:position2+3])
                    break


for orflen, orf in startstop_codon('ATGCGAATGTCGATCCGTAACtaaCGATTAA', 0):
    print(orflen, orf)

```

    24 ATGCGAATGTCGATCCGTAACTAA
    18 ATGTCGATCCGTAACTAA
    
