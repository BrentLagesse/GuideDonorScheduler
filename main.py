import argparse
import re
from dataclasses import dataclass

@dataclass
class MutationTracker:
    guide: int
    pam: int
    mutation: [[]]    # list of all mutations X to Y
    dna: str

GENE_START_BUFFER = 1000
GENE_END_BUFFER = 1000
# set up the argument parser so we can accept commandline arguments
argParser = argparse.ArgumentParser()
argParser.add_argument("-i", "--input", help="Input File in FSA format")
argParser.add_argument("-o", "--output", help="Output Filename base in FSA format")
args = argParser.parse_args()

in_file = args.input
out_base = args.output

#set up codon lookup table
codons = dict()

leu = ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG']
phe = ['TTT', 'TTC']
ile = ['AAT', 'ATC', 'ATA']
met = ['ATG']
val = ['GTT', 'GTC', 'GTA', 'GTG']
ser = ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC']
pro = ['CCT', 'CCC','CCA', 'CCG']
thr = ['ACT', 'ACC', 'ACA', 'ACG']
ala = ['GCT', 'GCC', 'GCA',  'GCG']
tyr = ['TAT', 'TAC']
his = ['CAT', 'CAC']
gln = ['CAA', 'CAG']
asn = ['AAT', 'AAC']
lys = ['AAA', 'AAG']
asp = ['GAT', 'GAC']
glu = ['GAA', 'GAG']
cys = ['TGT', 'TGC']
trp = ['TGG']
arg = ['CGT', 'CGC', 'CGA', 'CGG', 'AGA','AGG']
gly = ['GGT', 'GGC', 'GGA', 'GGG']

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

string_to_acid = dict()

for acid, strings in codons.items():
    for s in strings:
        string_to_acid[s] = acid


# TODO: fill this out when matt/talia get me preferences
preferred_mutations = dict()
#fake mutations for testing
preferred_mutations['leu'] = 'phe'
preferred_mutations['ala'] = 'tyr'
preferred_mutations['glu'] = 'arg'
preferred_mutations['cys'] = 'thr'
preferred_mutations['val'] = 'ile'


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
    gene_only = dna[GENE_START_BUFFER:len(dna)-GENE_END_BUFFER]
    gg_locs = [loc.start()+GENE_START_BUFFER - 1 for loc in re.finditer('(?=GG)', gene_only)]   # minus one accounts for the N of NGG

    return gg_locs

    #TODO:  add CC later after we get the code working on GG
    #cc_locs = [loc.start()+GENE_START_BUFFER for loc in re.finditer('(?=CC)', gene_only)]

    #return (gg_locs + cc_locs).sort()

def create_guides(dna, loc):
    pass
    # grab the 20 pairs before
    #TODO: do we just need the location?


def perform_mutation(candidate_dna, i, mutant):
    amino_acid_str = candidate_dna[i:3]  # TODO:  We don't want to mess with the seed if possible (the 10 before NGG)
    if amino_acid_str in string_to_acid:   # if this is something that isn't an amino acid, just quit
        amino_acid = string_to_acid[amino_acid_str]
    else:
        return False

    if mutant[0] == amino_acid:  # we found our target, lets make the swap!
        valid_mutations = codons[mutant[1]]  # get a list of valid mutations
        for mutation in valid_mutations:
            if mutation == candidate_dna[i:i+3]:  #This is what we already have, so it isn't a mutation
                continue
            if mutation[0] == 'G' and candidate_dna[i - 1] == 'G':  # this would introduce a GG on the front end
                continue
            if mutation[2] == 'G' and candidate_dna[i - 4] == 'G':  # this would introduce a GG on the back end
                continue
            # we are safe to make a swap here
            candidate_dna = candidate_dna[:i] + mutation + candidate_dna[i+3:]

            return True
    return False

#THis method will return a full dna string for each mutation as part of a MutationTracker type
#pam is the location of the first character of the pam
#mutant is key-val pair of mutant source to mutant destination [0] is key, [1] is value
def create_mutations(dna, pam, mutant):
    UPSTREAM = 100
    DOWNSTREAM = 100
    # TODO: take the 100 upstream and downstream.  mutate pam and introduce additional mutation

    # introduce mutation
    candidate_start = pam - UPSTREAM
    candidate_end = pam + 3 + DOWNSTREAM
    first_amino_acid_loc = int()
    for i in range(0, 2):    # We want to start on the first amino acid that is within our upstream range
        if (candidate_start - GENE_START_BUFFER + i) % 3 == 0:
            first_amino_acid_loc = candidate_start + i
    candidate_dna = dna[first_amino_acid_loc:candidate_end]   #grab starting from the first full amino acid
    # TODO:  Question for Talia/Matt -- Do we prefer to be far from the pam, close, middle?

    for i in range(0, UPSTREAM, 3):    # check upstream, then check downstream
        mutation_successful = perform_mutation(candidate_dna, i, mutant)
        if mutation_successful:
            break
    if not mutation_successful:
        print('Failed to find a valid place to mutate ' + mutant[0] + ' into ' + mutant[1])

    #TODO:  Implement downstream
    #for i in range(0, DOWNSTREAM, 3):
    #    pass

    # mutate pam

    # figure out the pam amino acid situation (does it split, and if so where)
    pam_case = (pam - GENE_START_BUFFER) % 3

    if pam_case == 0: # we only have a single amino acid
        pam_string = dna[pam:pam+3]
        pam_acid = string_to_acid[pam_string]
        pam_mutant=[pam_acid, pam_acid]
        mutation_successful = perform_mutation(candidate_dna, pam, pam_mutant)
        if not mutation_successful:
            print('Failed to find a valid replacement for the pam')
    elif pam_case == 1:  # the N is in one acid and the GG is in another
        # TODO: use a breadth first search, try downstream, then upstream, widening the number of changes
        pam_string_1 = dna[pam-2:pam+1]
        pam_string_2 = dna[pam+1:pam+4]
    elif pam_case == 2: # the NG is in one acid and the G is in another
        pam_string_1 = dna[pam-1:pam+2]
        pam_string_2 = dna[pam+2:pam+5]


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
    # I don't think we need to actually find the guides here since it is just pam - 20
    #guides = create_guides(dna, loc)
    #for g in guides:   # for each guide do each mutation
    for m in preferred_mutations.items():
        mutated_dna = create_mutations(dna, loc, m)
        all_mutations.append(mutated_dna)


# at this point, we have everything we need to output the results
write_results(frontmatter, all_mutations)



