import argparse
import re
from dataclasses import dataclass

@dataclass
class MutationTracker:
    guide: int
    pam: int
    mutation: [[]]    # list of all mutations X to Y
    dna: str
# set up the argument parser so we can accept commandline arguments
argParser = argparse.ArgumentParser()
argParser.add_argument("-i", "--input", help="Input File in FSA format")
argParser.add_argument("-o", "--output", help="Output Filename base in FSA format")
args = argParser.parse_args()

in_file = args.input
out_base = args.output

#set up codon lookup table
codons = dict()

leu = ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG']
phe = ['UUU', 'UUC']
ile = ['AAU', 'AUC', 'AUA']
met = ['AUG']
val = ['GUU', 'GUC', 'GUA', 'GUG']
ser = ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC']
pro = ['CCU', 'CCC','CCA', 'CCG']
thr = ['ACU', 'ACC', 'ACA', 'ACG']
ala = ['GCU', 'GCC', 'GCA',  'GCG']
tyr = ['UAU', 'UAC']
his = ['CAU', 'CAC']
gln = ['CAA', 'CAG']
asn = ['AAU', 'AAC']
lys = ['AAA', 'AAG']
asp = ['GAU', 'GAC']
glu = ['GAA', 'GAG']
cys = ['UGU', 'UGC']
trp = ['UGG']
arg = ['CGU', 'CGC', 'CGA', 'CGG', 'AGA','AGG']
gly = ['GGU', 'GGC', 'GGA', 'GGG']

codons['leu'] = leu
codons['phe'] = phe
codons['ile'] = ile
codons['met'] = met
codons['val'] = val
codons['ser'] = ser
codons['pro'] = pro
codons['thr'] = thr
codons['ala'] = ala
codons['tyr'] = tyr
codons['his'] = his
codons['gln'] = gln
codons['asn'] = asn
codons['lys'] = lys
codons['asp'] = asp
codons['glu'] = glu
codons['cys'] = cys
codons['trp'] = trp
codons['arg'] = arg
codons['gly'] = gly

# fill this out when matt/talia get me preferences
preferred_mutations = dict()

def get_dna():
    # open the input file
    input_data = open(in_file, 'r')

    #read data, separate the first line from the DNA and remove all the linebreaks from the DNA string
    all_data = input_data.read()
    input_data.close()
    frontmatter = all_data.partition('\n')[0]
    dna = all_data.partition('\n')[2]
    dna = dna.replace('\n', '')

    return frontmatter, dna

def get_locations(dna):
    #find all of the locations of the NGG or CCN triples
    gg_locs = [loc.start() for loc in re.finditer('(?=GG)', dna)]
    cc_locs = [loc.start() for loc in re.finditer('(?=CC)', dna)]

    return (gg_locs + cc_locs).sort()

def create_guides(dna, loc):
    # grab the 20 pairs before
    #TODO: do we just need the location?

#THis method will return a full dna string for each mutation as part of a MutationTracker type
def create_mutations(guide, mutant):
    # TODO: take the 100 upstream or downstream.  mutate pam and introduce additional mutation
    # TODO:  Ensure that no other GG or CC is introduced
    # TODO:  We don't want to mess with the seed if possible (the 10 before NGG)
    # TODO:  Track the decisions we made in this method so we can output them
    # TODO:  Fix naming of the file and frontmatter

def write_results(frontmatter, results):
    for mutation in results:
        out_file = out_base + '-' + str(mutation.pam) + '-' + mutation.mutation   #create a new file name
        f = open(out_file, 'w')
        new_frontmatter = frontmatter + ', mutated '
        for m in mutation.mutations:    # rewrite the comment to include what we did
            new_frontmatter += str(m[0]) + str(m[1])
        f.write(new_frontmatter)
        f.write('\n')
        f.write(results.dna)
        f.close()

frontmatter, dna = get_dna()
dna_locs = get_locations(dna)
all_mutations = []
for loc in dna_locs:
    guides = create_guides(dna, loc)
    for g in guides:   # for each guide do each mutation
        for m in preferred_mutations:
            mutated_dna = create_mutations(g, m)
            all_mutations.append(frontmatter, mutated_dna)


# at this point, we have everything we need to output the results
write_results(all_mutations)



