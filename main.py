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

@dataclass
class GlobalStats:
    failed_due_to_mutate: int
    failed_due_to_pam: int
    succeeded: int

gs = GlobalStats(0,0,0)



GENE_START_BUFFER = 1000
GENE_END_BUFFER = 1000
UP_ACIDS = 6
DOWN_ACIDS = 4
# set up the argument parser so we can accept commandline arguments
argParser = argparse.ArgumentParser()
argParser.add_argument("-i", "--input", help="Input File in FSA format")
argParser.add_argument("-o", "--output", help="Output Filename base in FSA format")
args = argParser.parse_args()

in_file = args.input
out_base = args.output

#set up codon lookup table
codons = dict()

# The one we prefer most is placed at the front of the list
leu = ['TTG', 'TTA', 'CTT', 'CTC', 'CTA', 'CTG']
phe = ['TTT', 'TTC']
ile = ['ATT', 'ATC', 'ATA']
met = ['ATG']
val = ['GTT', 'GTC', 'GTA', 'GTG']
ser = ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC']
pro = ['CCA', 'CCC', 'CCT', 'CCG']
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
arg = ['AGA', 'CGC', 'CGA', 'CGG', 'CGT', 'AGG']
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

def insert_extra_sequence(candidate_dna, guide):
    first = 'GACCGTGCGACTGGGCGTCTCGGATC'
    second = 'GTTTGAAGAGCATACGCTCTTCTTCT'
    third = 'ACATCGAGACGTGTCCCTGCCTTGCG'
    return first + guide + second + candidate_dna + third


def perform_mutation(candidate_dna, first_amino_acid_loc, i, mutant, keep_trying=False):
    amino_loc = first_amino_acid_loc - GENE_START_BUFFER + i
    amino_acid_str = candidate_dna[amino_loc: amino_loc+3]
    if amino_acid_str in string_to_acid:   # if this is something that isn't an amino acid, just quit
        amino_acid = string_to_acid[amino_acid_str]
    else:
        return False, None

    if mutant[0] == amino_acid:  # we found our target, lets make the swap!
        valid_mutations = codons[mutant[1]]  # get a list of valid mutations
        for mutation in valid_mutations:
            if mutation == candidate_dna[amino_loc:amino_loc+3]:  #This is what we already have, so it isn't a mutation
                continue

                #     if possible, I would like to avoid mutating the PAM to NAG.
            #    2. If you can’t mutate the PAM then mutate at least one location in the “seed” region (10 bases upstream of the PAM). The closer the silent mutation is to the PAM the better it works. We may decide that we want to do two silent mutations if we find that one in the seed region isn’t enough to prevent re-cutting.
    #3. Oftentimes the mutation we intend to make will be in the seed.
     #   a. If it is still possible to make a silent PAM mutation then that would be good (Although we are currently testing this and this parameter might change).
      #  b. If there is no way to make a silent PAM mutation and the mutation is more than 5 bases away from the PAM then the next best thing would be to make a silent mutation within the 5 bases upstream of the PAM.
       # c. If the mutation is within 5 bases from the PAM and you can’t make a silent PAM mutation, then I wouldn’t make any additional mutation.


                # This is from when I thought we couldn't introduce a new GG
            #if mutation[0] == 'G' and candidate_dna[amino_loc - 1] == 'G':  # this would introduce a GG on the front end
            #    if not keep_trying:
            #        continue
            #if mutation[2] == 'G' and candidate_dna[amino_loc - 4] == 'G':  # this would introduce a GG on the back end
            #    if not keep_trying:
            #        continue

            # we are safe to make a swap here
            candidate_dna = candidate_dna[:amino_loc] + mutation + candidate_dna[amino_loc+3:]


            return True, candidate_dna
    return False, None

#THis method will return a full dna string for each mutation as part of a MutationTracker type
#pam is the location of the first character of the pam
#mutant is key-val pair of mutant source to mutant destination [0] is key, [1] is value
def create_mutations(dna, pam, mutant):
    global gs
    # it seems like we are only looking at the 6 upstream and 4 downstream amino acids
    UPSTREAM = UP_ACIDS * 3
    DOWNSTREAM = DOWN_ACIDS * 3

    guide = dna[pam-20:pam]

    # introduce mutation
    # TODO:  Grab 132 base pairs surrounding the halfway point between the pam and the mutation
    # TODO:  Or we can grab 132 from the center of the guide
    # PAM + 10 is the center of the guide + 66 on each side of it
    candidate_start = pam - 10 - 66 #pam - UPSTREAM
    candidate_end = pam - 10 + 66 #pam + 3 + DOWNSTREAM
    first_amino_acid_loc = int()
    for i in range(0, 3):    # We want to start on the first amino acid that is within our upstream range
        if (candidate_start - GENE_START_BUFFER + i) % 3 == 0:
            first_amino_acid_loc = candidate_start + i
    while first_amino_acid_loc < GENE_START_BUFFER:   #ignore acids outside the gene
        first_amino_acid_loc += 3
    if first_amino_acid_loc > pam:   # if we can't get anything upstream, just start at the beginning of the gene
        first_amino_acid_loc = GENE_START_BUFFER

    #candidate_dna = dna[first_amino_acid_loc:candidate_end]   #grab starting from the first full amino acid
    candidate_dna = dna[candidate_start:candidate_end]
    mutation_successful = False
    mutation_location = -1
    if first_amino_acid_loc >= GENE_START_BUFFER and first_amino_acid_loc + 3 < pam:   # only do upstream if we are still in the gene
        for i in range(UPSTREAM - 3, -1, -3):    # check upstream, then check downstream
            if i + first_amino_acid_loc + 3 >= pam:   # don't go into the pam (TODO:  I think this is true)
                continue
            mutation_successful, temp_candidate_dna = perform_mutation(candidate_dna, first_amino_acid_loc, i, mutant)
            if mutation_successful:
                candidate_dna = temp_candidate_dna
                mutation_location = i
                break

    if candidate_end > (len(dna) - GENE_END_BUFFER):   # only do downstream if we are still in the gene
        pass
    #TODO:  Implement downstream
    #for i in range(0, DOWNSTREAM, 3):
    #    pass


    if not mutation_successful:
        print('Failed to find a valid place to mutate ' + mutant[0] + ' into ' + mutant[1])
        gs.failed_due_to_mutate += 1
        return None

    # mutate pam
    # figure out the pam amino acid situation (does it split, and if so where)
    pam_case = (pam - GENE_START_BUFFER) % 3
    mutation_successful = False
    if pam_case == 0: # we only have a single amino acid
        pam_string = dna[pam:pam+3]
        pam_acid = string_to_acid[pam_string]
        pam_mutant=[pam_acid, pam_acid]
        mutation_successful, candidate_dna = perform_mutation(candidate_dna, pam - first_amino_acid_loc, pam_mutant)
        if not mutation_successful:
            print('Failed to find a valid replacement for the pam')
            return None
    else:
        if pam_case == 1:  # the N is in one acid and the GG is in another so we can only replace down

            pam_string_up = None #dna[pam - 1:pam + 2]
            pam_string_down = dna[pam + 2:pam + 5]

        elif pam_case == 2:  # the NG is in one acid and the G is in another
            pam_string_up = dna[pam-2:pam+1]
            pam_string_down = dna[pam + 1:pam + 4]

        replaceable_pam = False
        pam_mutant_up = None
        if pam_string_up is not None and pam_string_up in string_to_acid:
            pam_acid_up = string_to_acid[pam_string_up]
            replaceable_pam = True
            pam_mutant_up = [pam_acid_up, pam_acid_up]
        elif pam_string_down in string_to_acid:
            pam_acid_down = string_to_acid[pam_string_down]
            replaceable_pam = True
            pam_mutant_down = [pam_acid_down, pam_acid_down]
        if not replaceable_pam:
            return None
        if pam_mutant_up is not None:
            mutation_successful, candidate_dna = perform_mutation(candidate_dna, pam - first_amino_acid_loc - pam_case, pam_mutant_up)
        if not mutation_successful:   #try downstream if upstream didn't work
            mutation_successful, candidate_dna = perform_mutation(candidate_dna, pam - first_amino_acid_loc - pam_case, pam_mutant_down)

    if mutation_successful:
        candidate_dna = insert_extra_sequence(candidate_dna, guide)
        #guide pam mutation mutationloc dna
        result = MutationTracker(pam - first_amino_acid_loc - 20, pam - first_amino_acid_loc, mutant, mutation_location, candidate_dna)
        gs.succeeded += 1
        return result
    else:
        gs.failed_due_to_pam += 1
        return None

    # TODO:  Track the decisions we made in this method so we can output them
    # TODO:  Fix naming of the file and frontmatter

def write_results(frontmatter, results):

    global gs

    print('\nfailed due to mutate: ' + str(gs.failed_due_to_mutate))
    print('failed due to pam: ' + str(gs.failed_due_to_pam))
    print('succeeded: ' + str(gs.succeeded))

    wb = Workbook()
    sheet1 = wb.add_sheet('Sheet 1')
    sheet1.write(0, 0, 'Guide Location')
    sheet1.write(0, 1, 'PAM Location')
    sheet1.write(0, 2, 'Mutation From')
    sheet1.write(0, 3, 'Mutation To')
    sheet1.write(0, 4, 'Mutation Location')
    sheet1.write(0, 5, 'Result')
    for i,mutation in enumerate(results):
        #out_file = out_base + '-' + str(mutation.pam) + '-' + mutation.mutation   #create a new file name
        #f = open(out_file, 'w')
        #new_frontmatter = frontmatter + ', mutated '
        sheet1.write(i + 1, 0, mutation.guide)
        sheet1.write(i + 1, 1, mutation.pam)
        sheet1.write(i + 1, 2, mutation.mutation[0])
        sheet1.write(i + 1, 3, mutation.mutation[1])
        sheet1.write(i + 1, 4, mutation.mutation_loc)
        if mutation.mutation_loc < mutation.pam:
            #TODO:  is it ok if the mutation location is in the guide?
            guide_bonus = 0  # in case we need to expand to include the guide
            if mutation.pam - mutation.mutation_loc < 20:
                guide_bonus = 20
            sheet1.write(i+1, 5, mutation.dna[mutation.mutation_loc-guide_bonus:mutation.pam+3])
        else:
            sheet1.write(i + 1, 5, mutation.dna[mutation.pam - 20 :mutation.mutation_loc])   # make sure to include the upstream guide if we go downstream
        #sheet1.write(1,i,str(mutation.mutation))
        #sheet1.write(2,i,str(m[1]))
        #f.write(new_frontmatter)
        #f.write('\n')
        #f.write(results.dna)
        #f.close()
    wb.save(out_base + '.xls')

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



