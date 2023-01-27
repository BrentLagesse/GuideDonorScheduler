import argparse
import re
from dataclasses import dataclass
import xlwt
from xlwt import Workbook

@dataclass
class MutationTracker:
    guide: int
    pam: int
    mutation: []
    mutation_loc: int
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
preferred_mutations['ile'] = 'phe'
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


def perform_mutation(candidate_dna, i, mutant, keep_trying=False):
    amino_acid_str = candidate_dna[i:i+3]  # TODO:  We don't want to mess with the seed if possible (the 10 before NGG)
    if amino_acid_str in string_to_acid:   # if this is something that isn't an amino acid, just quit
        amino_acid = string_to_acid[amino_acid_str]
    else:
        return False, None

    if mutant[0] == amino_acid:  # we found our target, lets make the swap!
        valid_mutations = codons[mutant[1]]  # get a list of valid mutations
        for mutation in valid_mutations:
            if mutation == candidate_dna[i:i+3]:  #This is what we already have, so it isn't a mutation
                continue
            if mutation[0] == 'G' and candidate_dna[i - 1] == 'G':  # this would introduce a GG on the front end
                if not keep_trying:
                    continue
            if mutation[2] == 'G' and candidate_dna[i - 4] == 'G':  # this would introduce a GG on the back end
                if not keep_trying:
                    continue

            # we are safe to make a swap here
            candidate_dna = candidate_dna[:i] + mutation + candidate_dna[i+3:]

            return True, candidate_dna
    return False, None

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
    for i in range(0, 3):    # We want to start on the first amino acid that is within our upstream range
        if (candidate_start - GENE_START_BUFFER + i) % 3 == 0:   # it is ok if this is negative for the early pam sites
            first_amino_acid_loc = candidate_start + i
    candidate_dna = dna[first_amino_acid_loc:candidate_end]   #grab starting from the first full amino acid
    # TODO:  Question for Talia/Matt -- Do we prefer to be far from the pam, close, middle?
    mutation_location = -1
    for i in range(0, UPSTREAM, 3):    # check upstream, then check downstream
        mutation_successful, temp_candidate_dna = perform_mutation(candidate_dna, i, mutant)
        if mutation_successful:
            candidate_dna = temp_candidate_dna
            mutation_location = i
            break
    if not mutation_successful:
        print('Failed to find a valid place to mutate ' + mutant[0] + ' into ' + mutant[1])
        return None

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
        mutation_successful, candidate_dna = perform_mutation(candidate_dna, pam - first_amino_acid_loc, pam_mutant)
        if not mutation_successful:
            print('Failed to find a valid replacement for the pam')
            return None
    else:
        if pam_case == 1:  # the N is in one acid and the GG is in another
            # TODO: use a breadth first search, try downstream, then upstream, widening the number of changes
            pam_string_up = dna[pam-2:pam+1]
            pam_string_down = dna[pam+1:pam+4]
        elif pam_case == 2:  # the NG is in one acid and the G is in another
            pam_string_up = dna[pam - 1:pam + 2]
            pam_string_down = dna[pam + 2:pam + 5]
        replaceable_pam = False
        if pam_string_up in string_to_acid:
            pam_acid_up = string_to_acid[pam_string_up]
            replaceable_pam = True
            pam_mutant_up = [pam_acid_up, pam_acid_up]
        elif pam_string_down in string_to_acid:    # this is an if and get rid of return None.  Only for testing
            return None
            pam_acid_down = string_to_acid[pam_string_down]
            replaceable_pam = True
            pam_mutant_down = [pam_acid_down, pam_acid_down]
        if not replaceable_pam:
            return None

        mutation_successful, candidate_dna = perform_mutation(candidate_dna, pam - first_amino_acid_loc - (3-pam_case), pam_mutant_up)

    if mutation_successful:
        result = MutationTracker(pam - first_amino_acid_loc, pam - first_amino_acid_loc - 20, mutant, mutation_location, candidate_dna)
        return result
    else:
        return None

    # TODO:  Track the decisions we made in this method so we can output them
    # TODO:  Fix naming of the file and frontmatter

def write_results(frontmatter, results):



    wb = Workbook()
    sheet1 = wb.add_sheet('Sheet 1')
    sheet1.write(0, 0, 'PAM Location')
    sheet1.write(0, 1, 'Mutation From')
    sheet1.write(0, 2, 'Mutation To')
    sheet1.write(0, 3, 'Mutation Location')
    sheet1.write(0, 4, 'Result')
    for i,mutation in enumerate(results):
        #out_file = out_base + '-' + str(mutation.pam) + '-' + mutation.mutation   #create a new file name
        #f = open(out_file, 'w')
        #new_frontmatter = frontmatter + ', mutated '
        sheet1.write(i + 1, 0, mutation.pam)
        sheet1.write(i + 1, 1, mutation.mutation[0])
        sheet1.write(i + 1, 2, mutation.mutation[1])
        sheet1.write(i + 1, 3, mutation.mutation_loc)
        sheet1.write(i+1, 4, mutation.dna)
        #sheet1.write(1,i,str(mutation.mutation))
        #sheet1.write(2,i,str(m[1]))
        #f.write(new_frontmatter)
        #f.write('\n')
        #f.write(results.dna)
        #f.close()
    wb.save('out_base.xls')

frontmatter, dna = get_dna()
dna_locs = get_locations(dna)
all_mutations = []
for loc in dna_locs:
    # I don't think we need to actually find the guides here since it is just pam - 20
    #guides = create_guides(dna, loc)
    #for g in guides:   # for each guide do each mutation
    for m in preferred_mutations.items():
        mutated_dna = create_mutations(dna, loc, m)
        if mutated_dna is not None:
            all_mutations.append(mutated_dna)


# at this point, we have everything we need to output the results
write_results(frontmatter, all_mutations)



