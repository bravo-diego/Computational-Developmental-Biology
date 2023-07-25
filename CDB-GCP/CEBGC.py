# Computational Evolutionary Biology Project   
   
   # Loading libraries
   
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from array import *
from random import *
from matplotlib.figure import Figure

    # Defining codons in a multidimensional array

codons = [[['u'], ['u'], ['u']], [['u'], ['u'], ['c']], [['u'], ['u'], ['a']], [['u'], ['u'], ['g']], [['c'], ['u'], ['u']], [['c'], ['u'], ['c']], [['c'], ['u'], ['a']], [['c'], ['u'], ['g']], [['a'], ['u'], ['u']], [['a'], ['u'], ['c']], [['a'], ['u'], ['a']], [['a'], ['u'], ['g']], [['g'], ['u'], ['u']], [['g'], ['u'], ['c']], [['g'], ['u'], ['a']], [['g'], ['u'], ['g']], [['u'], ['c'], ['u']], [['u'], ['c'], ['c']], [['u'], ['c'], ['a']], [['u'], ['c'], ['g']], [['c'], ['c'], ['u’']], [['c'], ['c'], ['c']], [['c'], ['c'], ['a']], [['c'], ['c'], ['g']], [['a'], ['c'], ['u']], [['a'], ['c'], ['c']], [['a'], ['c'], ['a']], [['a'], ['c'], ['g']], [['g'], ['c'], ['u']], [['g'], ['c'], ['c']], [['g'], ['c'], ['a']], [['g'], ['c'], ['g']], [['u'], ['a'], ['u']], [['u'], ['a'], ['c']], [['u'], ['a'], ['a']], [['u'], ['a'], ['g']], [['c'], ['a'], ['u']], [['c'], ['a'], ['c']], [['c'], ['a'], ['a']], [['c'], ['a'], ['g']], [['a'], ['a'], ['u']], [['a'], ['a'], ['c']], [['a'], ['a'], ['a']], [['a'], ['a'], ['g']], [['g'], ['a'], ['u']], [['g'], ['a'], ['c']], [['g'], ['a'], ['a']], [['g'], ['a'], ['g']], [['u'], ['g'], ['u']], [['u'], ['g'], ['c']], [['u'], ['g'], ['a']], [['u'], ['g'], ['g']], [['c'], ['g'], ['u']], [['c'], ['g'], ['c']], [['c'], ['g'], ['a']], [['c'], ['g'], ['g']], [['a'], ['g'], ['u']], [['a'], ['g'], ['c']], [['a'], ['g'], ['a']], [['a'], ['g'], ['g']], [['g'], ['g'], ['u']], [['g'], ['g'], ['c']], [['g'], ['g'], ['a']], [['g'], ['g'], ['g']]]

    # Defining amino acids in an array

amino_acids = ['Phe', 'Phe', 'Leu', 'Leu', 'Leu', 'Leu', 'Leu', 'Leu', 'Ile', 'Ile', 'Ile', 'Met', 'Val', 'Val', 'Val', 'Val', 'Ser', 'Ser', 'Ser', 'Ser', 'Pro', 'Pro', 'Pro', 'Pro', 'Thr', 'Thr', 'Thr', 'Thr', 'Ala', 'Ala', 'Ala', 'Ala', 'Tyr', 'Tyr', 'End', 'End', 'His', 'His', 'Gln', 'Gln', 'Asn', 'Asn', 'Lys', 'Lys', 'Asp', 'Asp', 'Glu', 'Glu', 'Cys', 'Cys', 'End', 'Trp', 'Arg', 'Arg', 'Arg', 'Arg', 'Ser', 'Ser', 'Arg', 'Arg', 'Gly', 'Gly', 'Gly', 'Gly']

print("There are {} codons, {} amino acids and {} stop codon.".format(len(codons), len(np.unique(np.array(amino_acids))) - 1, len(np.unique(np.array(amino_acids))) - 20))

genetic_code = [] # original genetic code
for i in range(len(codons)):
    genetic_code.append([[], []])
    genetic_code[i][0].append(amino_acids[i])
    for k in range(len(codons[0])):
        genetic_code[i][1].append(codons[i][k])

    # Function to generate random codes

def random_code(empty_code): 
    for i in range(len(codons)):
        empty_code.append([[], []])
        for k in range(len(codons[0])):
            empty_code[i][1].append(codons[i][k])
    amino_acids_copy = amino_acids[:] # amino acids array (copy)
    total_codons = len(codons)
    for i in range(len(codons)):
        random_position = randrange(total_codons) 
        empty_code[i][0].append(amino_acids_copy[random_position])
        del amino_acids_copy[random_position] # delete amino acid in array (copy)
        total_codons -= 1
    return empty_code # random genetic code
    
    # Functions to count mutations around genetic code *by sections
    
def uu_section(genetic_code):
    mutations = 0 # number of mutations
    x = 3
    y = 3
    z = 3
    i = 0 # start position
    k = i + 1 # switch amino acid in same section
    l = i + 4 # switch amino acid in horizontal section
    m = i + 16 # switch amino acid in vertical section
    for A in range(4): # amino acids in section (i.e. Phe, Phe, Leu, Leu)
        amino_acid = genetic_code[i][0]
        for B in range(x): # compare amino acids (same section)
            if amino_acid != genetic_code[k][0]:
                mutations += 1
                k += 1
            else:
                k += 1
        x -= 1
        for C in range(y): # compare amino acids (horizontal sections; 4 out of 4)
            if amino_acid != genetic_code[l][0]:
                mutations += 1
                l += 4
            else:
                l += 4
        for D in range(z): # compare amino acids (vertical sections; 4 out of 4)
            if amino_acid != genetic_code[m][0]:
                mutations += 1
                m += 16
            else:
                m += 16
        i += 1
        k = i + 1
        l = i + 4
        m = i + 16
    return mutations
    
def uc_section(genetic_code):
    mutations = 0 # number of mutations
    x = 3
    y = 2
    z = 3
    i = 4 # start position
    k = i + 1 # switch amino acid in same section
    l = i + 4 # switch amino acid in horizontal section
    m = i + 16 # switch amino acid in vertical section
    for A in range(4): # amino acids in section (i.e. Leu, Leu, Leu, Leu)
        amino_acid = genetic_code[i][0]
        for B in range(x): # compare amino acids (same section)
            if amino_acid != genetic_code[k][0]:
                mutations += 1
                k += 1
            else:
                k += 1
        x -= 1
        for C in range(y): # compare amino acids (horizontal sections; 3 out of 4)
            if amino_acid != genetic_code[l][0]:
                mutations += 1
                l += 4
            else:
                l += 4
        for D in range(z): # compare amino acids (vertical sections; 4 out of 4)
            if amino_acid != genetic_code[m][0]:
                mutations += 1
                m += 16
            else:
                m += 16
        i += 1
        k = i + 1
        l = i + 4
        m = i + 16
    return mutations
    
def ua_section(genetic_code):
    mutations = 0 # number of mutations
    x = 3
    y = 1
    z = 3
    i = 8 # start position
    k = i + 1 # switch amino acid in same section
    l = i + 4 # switch amino acid in horizontal section
    m = i + 16 # switch amino acid in vertical section
    for A in range(4): # amino acids in section (i.e. Ile, Ile, Ile, Met)
        amino_acid = genetic_code[i][0]
        for B in range(x): # compare amino acids (same section)
            if amino_acid != genetic_code[k][0]:
                mutations += 1
                k += 1
            else:
                k += 1
        x -= 1
        for C in range(y): # compare amino acids (horizontal sections; 2 out of 4)
            if amino_acid != genetic_code[l][0]:
                mutations += 1
                l += 4
            else:
                l += 4
        for D in range(z): # compare amino acids (vertical sections; 4 out of 4)
            if amino_acid != genetic_code[m][0]:
                mutations += 1
                m += 16
            else:
                m += 16
        i += 1
        k = i + 1
        l = i + 4
        m = i + 16
    return mutations
    
def ug_section(genetic_code):
    mutations = 0 # number of mutations
    x = 3
    y = 0
    z = 3
    i = 12 # start position
    k = i + 1 # switch amino acid in same section
    l = i + 4 # switch amino acid in horizontal section
    m = i + 16 # switch amino acid in vertical section
    for A in range(4): # amino acids in section (i.e. Val, Val, Val, Val)
        amino_acid = genetic_code[i][0]
        for B in range(x): # compare amino acids (same section)
            if amino_acid != genetic_code[k][0]:
                mutations += 1
                k += 1
            else:
                k += 1
        x -= 1
        for C in range(y): # compare amino acids (horizontal sections; 1 out of 4)
            if amino_acid != genetic_code[l][0]:
                mutations += 1
                l += 4
            else:
                l += 4
        for D in range(z): # compare amino acids (vertical sections; 4 out of 4)
            if amino_acid != genetic_code[m][0]:
                mutations += 1
                m += 16
            else:
                m += 16
        i += 1
        k = i + 1
        l = i + 4
        m = i + 16
    return mutations
    
def cu_section(genetic_code):
    mutations = 0 # number of mutations
    x = 3
    y = 3
    z = 2
    i = 16 # start position
    k = i + 1 # switch amino acid in same section
    l = i + 4 # switch amino acid in horizontal section
    m = i + 16 # switch amino acid in vertical section
    for A in range(4): # amino acids in section (i.e. Ser, Ser, Ser, Ser)
        amino_acid = genetic_code[i][0]
        for B in range(x): # compare amino acids (same section)
            if amino_acid != genetic_code[k][0]:
                mutations += 1
                k += 1
            else:
                k += 1
        x -= 1
        for C in range(y): # compare amino acids (horizontal sections; 4 out of 4)
            if amino_acid != genetic_code[l][0]:
                mutations += 1
                l += 4
            else:
                l += 4
        for D in range(z): # compare amino acids (vertical sections; 3 out of 4)
            if amino_acid != genetic_code[m][0]:
                mutations += 1
                m += 16
            else:
                m += 16
        i += 1
        k = i + 1
        l = i + 4
        m = i + 16
    return mutations
    
def cc_section(genetic_code):
    mutations = 0 # number of mutations
    x = 3
    y = 2
    z = 2
    i = 20 # start position
    k = i + 1 # switch amino acid in same section
    l = i + 4 # switch amino acid in horizontal section
    m = i + 16 # switch amino acid in vertical section
    for A in range(4): # amino acids in section (i.e. Pro, Pro, Pro, Pro)
        amino_acid = genetic_code[i][0]
        for B in range(x): # compare amino acids (same section)
            if amino_acid != genetic_code[k][0]:
                mutations += 1
                k += 1
            else:
                k += 1
        x -= 1
        for C in range(y): # compare amino acids (horizontal sections; 3 out of 4)
            if amino_acid != genetic_code[l][0]:
                mutations += 1
                l += 4
            else:
                l += 4
        for D in range(z): # compare amino acids (vertical sections; 3 out of 4)
            if amino_acid != genetic_code[m][0]:
                mutations += 1
                m += 16
            else:
                m += 16
        i += 1
        k = i + 1
        l = i + 4
        m = i + 16
    return mutations
    
def ca_section(genetic_code):
    mutations = 0 # number of mutations
    x = 3
    y = 1
    z = 2
    i = 24 # start position
    k = i + 1 # switch amino acid in same section
    l = i + 4 # switch amino acid in horizontal section
    m = i + 16 # switch amino acid in vertical section
    for A in range(4): # amino acids in section (i.e. Thr, Thr, Thr, The)
        amino_acid = genetic_code[i][0]
        for B in range(x): # compare amino acids (same section)
            if amino_acid != genetic_code[k][0]:
                mutations += 1
                k += 1
            else:
                k += 1
        x -= 1
        for C in range(y): # compare amino acids (horizontal sections; 2 out of 4)
            if amino_acid != genetic_code[l][0]:
                mutations += 1
                l += 4
            else:
                l += 4
        for D in range(z): # compare amino acids (vertical sections; 3 out of 4)
            if amino_acid != genetic_code[m][0]:
                mutations += 1
                m += 16
            else:
                m += 16
        i += 1
        k = i + 1
        l = i + 4
        m = i + 16
    return mutations    
    
def cg_section(genetic_code):
    mutations = 0 # number of mutations
    x = 3
    y = 0
    z = 2
    i = 28 # start position
    k = i + 1 # switch amino acid in same section
    l = i + 4 # switch amino acid in horizontal section
    m = i + 16 # switch amino acid in vertical section
    for A in range(4): # amino acids in section (i.e. Ala, Ala, Ala, Ala)
        amino_acid = genetic_code[i][0]
        for B in range(x): # compare amino acids (same section)
            if amino_acid != genetic_code[k][0]:
                mutations += 1
                k += 1
            else:
                k += 1
        x -= 1
        for C in range(y): # compare amino acids (horizontal sections; 1 out of 4)
            if amino_acid != genetic_code[l][0]:
                mutations += 1
                l += 4
            else:
                l += 4
        for D in range(z): # compare amino acids (vertical sections; 3 out of 4)
            if amino_acid != genetic_code[m][0]:
                mutations += 1
                m += 16
            else:
                m += 16
        i += 1
        k = i + 1
        l = i + 4
        m = i + 16
    return mutations 
    
def au_section(genetic_code):
    mutations = 0 # number of mutations
    x = 3
    y = 3
    z = 1
    i = 32 # start position
    k = i + 1 # switch amino acid in same section
    l = i + 4 # switch amino acid in horizontal section
    m = i + 16 # switch amino acid in vertical section
    for A in range(4): # amino acids in section (i.e. Tyr, Tyr, End, End)
        amino_acid = genetic_code[i][0]
        for B in range(x): # compare amino acids (same section)
            if amino_acid != genetic_code[k][0]:
                mutations += 1
                k += 1
            else:
                k += 1
        x -= 1
        for C in range(y): # compare amino acids (horizontal sections; 4 out of 4)
            if amino_acid != genetic_code[l][0]:
                mutations += 1
                l += 4
            else:
                l += 4
        for D in range(z): # compare amino acids (vertical sections; 2 out of 4)
            if amino_acid != genetic_code[m][0]:
                mutations += 1
                m += 16
            else:
                m += 16
        i += 1
        k = i + 1
        l = i + 4
        m = i + 16
    return mutations 
    
def ac_section(genetic_code):
    mutations = 0 # number of mutations
    x = 3
    y = 2
    z = 1
    i = 36 # start position
    k = i + 1 # switch amino acid in same section
    l = i + 4 # switch amino acid in horizontal section
    m = i + 16 # switch amino acid in vertical section
    for A in range(4): # amino acids in section (i.e. His, His, Gln, Gln)
        amino_acid = genetic_code[i][0]
        for B in range(x): # compare amino acids (same section)
            if amino_acid != genetic_code[k][0]:
                mutations += 1
                k += 1
            else:
                k += 1
        x -= 1
        for C in range(y): # compare amino acids (horizontal sections; 3 out of 4)
            if amino_acid != genetic_code[l][0]:
                mutations += 1
                l += 4
            else:
                l += 4
        for D in range(z): # compare amino acids (vertical sections; 2 out of 4)
            if amino_acid != genetic_code[m][0]:
                mutations += 1
                m += 16
            else:
                m += 16
        i += 1
        k = i + 1
        l = i + 4
        m = i + 16
    return mutations 
    
def aa_section(genetic_code):
    mutations = 0 # number of mutations
    x = 3
    y = 1
    z = 1
    i = 40 # start position
    k = i + 1 # switch amino acid in same section
    l = i + 4 # switch amino acid in horizontal section
    m = i + 16 # switch amino acid in vertical section
    for A in range(4): # amino acids in section (i.e. Asn, Asn, Lys, Lys)
        amino_acid = genetic_code[i][0]
        for B in range(x): # compare amino acids (same section)
            if amino_acid != genetic_code[k][0]:
                mutations += 1
                k += 1
            else:
                k += 1
        x -= 1
        for C in range(y): # compare amino acids (horizontal sections; 2 out of 4)
            if amino_acid != genetic_code[l][0]:
                mutations += 1
                l += 4
            else:
                l += 4
        for D in range(z): # compare amino acids (vertical sections; 2 out of 4)
            if amino_acid != genetic_code[m][0]:
                mutations += 1
                m += 16
            else:
                m += 16
        i += 1
        k = i + 1
        l = i + 4
        m = i + 16
    return mutations
    
def ag_section(genetic_code):
    mutations = 0 # number of mutations
    x = 3
    y = 0
    z = 1
    i = 44 # start position
    k = i + 1 # switch amino acid in same section
    l = i + 4 # switch amino acid in horizontal section
    m = i + 16 # switch amino acid in vertical section
    for A in range(4): # amino acids in section (i.e. Asp, Asp, Glu Glu)
        amino_acid = genetic_code[i][0]
        for B in range(x): # compare amino acids (same section)
            if amino_acid != genetic_code[k][0]:
                mutations += 1
                k += 1
            else:
                k += 1
        x -= 1
        for C in range(y): # compare amino acids (horizontal sections; 1 out of 4)
            if amino_acid != genetic_code[l][0]:
                mutations += 1
                l += 4
            else:
                l += 4
        for D in range(z): # compare amino acids (vertical sections; 2 out of 4)
            if amino_acid != genetic_code[m][0]:
                mutations += 1
                m += 16
            else:
                m += 16
        i += 1
        k = i + 1
        l = i + 4
        m = i + 16
    return mutations
    
def gu_section(genetic_code):
    mutations = 0 # number of mutations
    x = 3
    y = 3
    z = 0
    i = 48 # start position
    k = i + 1 # switch amino acid in same section
    l = i + 4 # switch amino acid in horizontal section
    m = i + 16 # switch amino acid in vertical section
    for A in range(4): # amino acids in section (i.e. Cys, Cys, End, Trp)
        amino_acid = genetic_code[i][0]
        for B in range(x): # compare amino acids (same section)
            if amino_acid != genetic_code[k][0]:
                mutations += 1
                k += 1
            else:
                k += 1
        x -= 1
        for C in range(y): # compare amino acids (horizontal sections; 4 out of 4)
            if amino_acid != genetic_code[l][0]:
                mutations += 1
                l += 4
            else:
                l += 4
        for D in range(z): # compare amino acids (vertical sections; 1 out of 4)
            if amino_acid != genetic_code[m][0]:
                mutations += 1
                m += 16
            else:
                m += 16
        i += 1
        k = i + 1
        l = i + 4
        m = i + 16
    return mutations
    
def gc_section(genetic_code):
    mutations = 0 # number of mutations
    x = 3
    y = 2
    z = 0
    i = 52 # start position
    k = i + 1 # switch amino acid in same section
    l = i + 4 # switch amino acid in horizontal section
    m = i + 16 # switch amino acid in vertical sgenetic_codeection
    for A in range(4): # amino acids in section (i.e. Arg, Arg, Arg, Arg)
        amino_acid = genetic_code[i][0]
        for B in range(x): # compare amino acids (same section)
            if amino_acid != genetic_code[k][0]:
                mutations += 1
                k += 1
            else:
                k += 1
        x -= 1
        for C in range(y): # compare amino acids (horizontal sections; 3 out of 4)
            if amino_acid != genetic_code[l][0]:
                mutations += 1
                l += 4
            else:
                l += 4
        for D in range(z): # compare amino acids (vertical sections; 1 out of 4)
            if amino_acid != genetic_code[m][0]:
                mutations += 1
                m += 16
            else:
                m += 16
        i += 1
        k = i + 1
        l = i + 4
        m = i + 16
    return mutations
    
def ga_section(genetic_code):
    mutations = 0 # number of mutations
    x = 3
    y = 1
    z = 0
    i = 56 # start position
    k = i + 1 # switch amino acid in same section
    l = i + 4 # switch amino acid in horizontal section
    m = i + 16 # switch amino acid in vertical section
    for A in range(4): # amino acids in section (i.e. Ser, Ser, Arg, Arg)
        amino_acid = genetic_code[i][0]
        for B in range(x): # compare amino acids (same section)
            if amino_acid != genetic_code[k][0]:
                mutations += 1
                k += 1
            else:
                k += 1
        x -= 1
        for C in range(y): # compare amino acids (horizontal sections; 2 out of 4)
            if amino_acid != genetic_code[l][0]:
                mutations += 1
                l += 4
            else:
                l += 4
        for D in range(z): # compare amino acids (vertical sections; 1 out of 4)
            if amino_acid != genetic_code[m][0]:
                mutations += 1
                m += 16
            else:
                m += 16
        i += 1
        k = i + 1
        l = i + 4
        m = i + 16
    return mutations
    
def gg_section(genetic_code):
    mutations = 0 # number of mutations
    x = 3
    y = 0
    z = 0
    i = 60 # start position
    k = i + 1 # switch amino acid in same section
    l = i + 4 # switch amino acid in horizontal section
    m = i + 16 # switch amino acid in vertical section
    for A in range(4): # amino acids in section (i.e. Gly, Gly, Gly, Gly)
        amino_acid = genetic_code[i][0]
        for B in range(x): # compare amino acids (same section)
            if amino_acid != genetic_code[k][0]:
                mutations += 1
                k += 1
            else:
                k += 1
        x -= 1
        for C in range(y): # compare amino acids (horizontal sections; 1 out of 4)
            if amino_acid != genetic_code[l][0]:
                mutations += 1
                l += 4
            else:
                l += 4
        for D in range(z): # compare amino acids (vertical sections; 1 out of 4)
            if amino_acid != genetic_code[m][0]:
                mutations += 1
                m += 16
            else:
                m += 16
        i += 1
        k = i + 1
        l = i + 4
        m = i + 16
    return mutations
    
    # Function to sum up all the mutations around genetic code
    
def total_mutations(genetic_code):
    total = uu_section(genetic_code) + uc_section(genetic_code) + ua_section(genetic_code) + ug_section(genetic_code) + cu_section(genetic_code) + cc_section(genetic_code) + ca_section(genetic_code) + cg_section(genetic_code) + au_section(genetic_code) + ac_section(genetic_code) + aa_section(genetic_code) + ag_section(genetic_code) + gu_section(genetic_code) + gc_section(genetic_code) + ga_section(genetic_code) + gg_section(genetic_code)
    return total
    
    # Function to generate multiple genetic codes and sum up all the mutations around them

mutations_frequencies = []

def frequencies(results):
    for i in range(len(codons)*10000):
        empty_code = []
        random_code(empty_code)
        results.append(total_mutations(empty_code))
    return results

results = frequencies(mutations_frequencies)

print("{} genetic codes were generated and evaluated, with a maximum number of {} mutations and a minimum of {} mutations".format(len(results), max(results), min(results)))

results_dataframe = pd.DataFrame(results, columns=['# Mutations'])
results_dataframe.describe()

    # Histogram

plt.figure()

plt.title("Distribution of mutations of random codes", fontsize = 12)
plt.xlabel("Number of mutations per code", fontsize = 8)
plt.ylabel("Frequency", fontsize = 8)

ax = plt.gca()
ax.spines[['top', 'right']].set_visible(False)
ax.hist(results, bins = 27, facecolor = '#ADD8E6', 
        linewidth = 0.5)
ax.set_xticks(np.arange(min(results), max(results), 5))
ax.set_yticks([0, 20000, 40000, 60000, 80000, 100000])
