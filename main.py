import argparse
import re
from dataclasses import dataclass
import sys

import fastaparser
# fastaparser, xlrd, xlutils and xlwt modules need to be installed.
# Can be done via pip install <module>

#import fastaparser
import xlwt
from xlwt import Workbook
import xlrd
import xlutils

# Other modules
import config

@dataclass
class MutationTracker:
    guide: str
    pam: int
    mutation: []
    mutation_loc: int
    dna: str
    complement: bool
    pam_location_in_gene: int
    distance_from_pam: int
    original_pam: str
    justification: str

@dataclass
class GlobalStats:
    failed_due_to_mutate: int
    failed_due_to_pam: int
    failed_due_to_guide_library: int
    failed_due_to_rank: int
    guides_used_from_rank_file: int
    succeeded: int


@dataclass
class PrebuiltGuide:
    gene: str
    id: str
    guide: str
    mutation_loc: int
    priority: float
    compliment: bool


gs = GlobalStats(0, 0, 0, 0, 0, 0)
guide_lib = []

# set up the argument parser so we can accept commandline arguments
argParser = argparse.ArgumentParser()
argParser.add_argument("-i", "--input", help="Input File in FSA format")
argParser.add_argument("-o", "--output", help="Output Filename base in FSA format")
args = argParser.parse_args()

in_files = args.input
out_base = args.output

# set up codon lookup table
codons = dict()

# all codons are defined in config.py in the CODON SETUP section
codons['leu'] = config.leu
codons['phe'] = config.phe
codons['ile'] = config.ile
codons['met'] = config.met
codons['val'] = config.val
codons['ser'] = config.ser
codons['pro'] = config.pro
codons['thr'] = config.thr
codons['ala'] = config.ala
codons['tyr'] = config.tyr
codons['his'] = config.his
codons['gln'] = config.gln
codons['asn'] = config.asn
codons['lys'] = config.lys
codons['asp'] = config.asp
codons['glu'] = config.glu
codons['cys'] = config.cys
codons['trp'] = config.trp
codons['arg'] = config.arg
codons['gly'] = config.gly
codons['nns'] = config.nns
codons['stop'] = config.stop

invert_mapping = dict()
invert_mapping['T'] = 'A'
invert_mapping['A'] = 'T'
invert_mapping['G'] = 'C'
invert_mapping['C'] = 'G'
invert_mapping['N'] = 'R'
invert_mapping['S'] = 'B'
invert_mapping['Z'] = 'Z'  # FOR DEBUG

string_to_acid = dict()

for acid, strings in codons.items():
    for s in strings:
        string_to_acid[s] = acid


def get_dna():
    # Checks if default input/output files should be used
    confirm_input_args()
    input_data = []

    # If the input file fails to load, the program exits with an error message
    try:
        # open the input files and insert into input_data
        frontmatter = []
        for line in in_files:
            if config.KILL_MODE:
                fasta_file = open(line, 'r')
                r = fastaparser.Reader(fasta_file)
                for seq in r:
                    # seq is a FastaSequence object
                    # print('ID:', seq.id)
                    # print('Description:', seq.description)
                    # print('Sequence:', seq.sequence_as_string())
                    # print()
                    input_data.append(seq.sequence_as_string())
                    frontmatter.append(seq.id)
                fasta_file.close()

                return frontmatter, input_data, None
            else:
                input_data.append(open(line, 'r'))
    except:
        if (config.QUIT_ON_NO_DATA):
            print("Opening DNA file failed, program exiting.")
            sys.exit()
    print(input_data)

    # read data, separate the first line from the DNA and remove all the linebreaks from the DNA string
    all_data = []

    for data in input_data:
        all_data.append(data.read())
        data.close()

    frontmatter = []
    dna = []

    for data in all_data:
        frontmatter.append(data.partition('\n')[0])
        temp = data.partition('\n')[2]
        dna.append(temp.replace('\n', ''))

    # PROCESS EXITS IF GUIDE LIBRARY IS NOT IN USE
    if (config.USE_GUIDE_LIBRARY != True):
        return frontmatter, dna, None

    workbook = xlrd.open_workbook(config.GUIDE_LIBRARY_INPUT_FILE + ".xls")
    worksheet = workbook.sheet_by_index(0)

    running = True
    i = 1

    guide_data = []

    while running:
        i += 1
        cur = worksheet.cell_value(i, 0)

        if (cur != ""):
            entry = PrebuiltGuide(cur,
                                  worksheet.cell_value(i, 1),
                                  worksheet.cell_value(i, 2),
                                  worksheet.cell_value(i, 3),
                                  worksheet.cell_value(i, 5),
                                  False)
            entry_inv = PrebuiltGuide(cur,
                                      worksheet.cell_value(i, 1),
                                      worksheet.cell_value(i, 4),
                                      worksheet.cell_value(i, 3),
                                      worksheet.cell_value(i, 6),
                                      True)

            guide_data.append(entry)

        if (cur == config.GUIDE_LIBRARY_EOF):
            running = False

    return frontmatter, dna, guide_data


# Returns locations of NGG and NCC triples as a set of arrays : [NGGs, NCCs]
def get_locations(dna):
    # find all of the locations of the NGG or CCN triples
    gene_only = dna[config.GENE_START_BUFFER:len(dna) - config.GENE_END_BUFFER]

    gg_locs = [loc.start() + config.GENE_START_BUFFER - 1 for loc in
               re.finditer('(?=GG)', gene_only)]  # minus one accounts for the N of NGG
    cc_locs = []
    if not config.KILL_MODE:
        cc_locs = [loc.start() + config.GENE_START_BUFFER for loc in
                   re.finditer('(?=CC)', gene_only)]
    # Make sure its using nGG pams for reverse compliment, look into this and see why nCCs are used

    return [gg_locs, cc_locs]


# Returns the 20 base-pairs before the PAM location
def create_guides(dna, loc, complement):
    if complement:
        return invert_dna(dna[loc+3:loc+23])
    else:
        return dna[loc - config.GUIDE_LENGTH:loc]


# Returns a list of guides subtracting all duplicates, leaving guides of either highest priority, or the first instance of a gene
def filter_guides(_guide_list):
    # TODO
    _return_list = []

    print("Filtering guide library guides.")
    print("Original Count: " + str(len(_guide_list)))
    for guide_lib_object in _guide_list:

        ind = -1
        for i in range(len(_return_list)):
            glo = _return_list[i]
            if glo.guide == guide_lib_object.guide and glo.gene == guide_lib_object.gene:
                ind = i

                continue

        if ind != -1:
            glo = _return_list[ind]
            if guide_lib_object.priority > glo.priority:
                _return_list[ind] = guide_lib_object
        else:
            _return_list.append(guide_lib_object)
    print("New Count: " + str(len(_return_list)))
    return _return_list


# Returns the earliest possible guide modified with an end codon to kill the protein
def create_kill_guide(mutation_tracker_to_modify):
    m = mutation_tracker_to_modify
    m2 = MutationTracker(m.guide, m.pam, m.mutation, m.mutation_loc, m.dna, m.complement, m.pam_location_in_gene, m.distance_from_pam, m.justification, '')

    m2.mutation = config.KILL_MUTATION
    m2.dna = str(m.dna[:m.mutation_loc]) + str(config.KILL_MUTATION[1]) + str(m.dna[m.mutation_loc + 3:])

    return m2


# Returns if the given guide is held within the guide library
def is_guide_in_library(guide, guide_library):
    for guide_entry in guide_library:
        if guide_entry.guide == guide:
            return True
    return False


# Returns if the given mutation tracker is allowed
def is_mutation_permitted(guide, tracker, guide_library):
    for guide_entry in guide_library:
        if guide_entry.guide == guide:
            if int(guide_entry.mutation_loc) == int(tracker.mutation_loc):
                return True
    return False


# Combines the candidate dna with the guide, as well as the extra sequences defined within config
def insert_extra_sequence(candidate_dna, guide):
    first = config.first_sequence
    second = config.second_sequence
    third = config.third_sequence
    return first + guide + second + candidate_dna + third


def perform_mutation(candidate_dna, first_amino_acid_loc, pam_case, mutant, decision_path, keep_trying=False, distance_from_pam=0,
                     mutation_location=-1, complement=False, down=False, total_mutations=0):

    #global amino_acid_number
    if complement:
        pam_char_one = 'C'
        pam_char_two = 'T'
        avoid_one = 'TC'
        avoid_two = 'CT'
        pam_indicator = 'CC'
        offset = 3
    else:
        pam_char_one = 'G'
        pam_char_two = 'A'
        avoid_one = 'GA'
        avoid_two = 'AG'
        pam_indicator = 'GG'
        offset = -3

    # config for silent mutations
    try_until = 8  # this is our default for now, but I'm still doing the checks in case we change this later
    if config.SILENT_MUTATION_MODE == 'seed':
        try_until = 8
    elif config.SILENT_MUTATION_MODE == 'guide':
        try_until = config.GUIDE_LENGTH


    actual_mutation = [mutant[0], mutant[1]]
    amino_acid_str = candidate_dna[first_amino_acid_loc: first_amino_acid_loc + 3]
    if first_amino_acid_loc == mutation_location:  # we would be replacing the mutation, so skip it
        #if distance_from_pam <= 5:  # don't mutate
        if config.VERBOSE_EXECUTION:
            print('we did not mutate the pam because the mutation was within 5 base pairs of pam')
        mutant[0] = '*'  # these two *s force any silent mutation
        mutant[1] = '*'
        decision_path += "Did not Mutate " + amino_acid_str + " because it is the real mutation. "
        return perform_mutation(candidate_dna, first_amino_acid_loc + offset, pam_case, mutant, decision_path, mutation_location=mutation_location,
                                distance_from_pam=distance_from_pam + 3, down=down, complement=complement, total_mutations=total_mutations)

        #return False, None, None, actual_mutation, decision_path
    if config.PRINT_MUTATION_CHECKS:
        print(amino_acid_str)
    if amino_acid_str in string_to_acid:  # if this is something that isn't an amino acid, just quit
        amino_acid = string_to_acid[amino_acid_str]
    else:
        if amino_acid_str in config.stop:
            if config.VERBOSE_EXECUTION:
                print('Ran into a stop codon: ' + amino_acid_str)
            return False, None, None, actual_mutation, decision_path
        else:
            if config.VERBOSE_EXECUTION:
                print('Somehow we ran into something that was not an amino acid: ' + amino_acid_str)
            return False, None, None, actual_mutation, decision_path

    if config.PRINT_MUTATION_CHECKS:
        print("Currently checking " + str(amino_acid) + " for " + str(mutant[0]) + ".")

    if mutant[0] == amino_acid or mutant[0] == '*' or mutant[0] == mutant[1]:  # we found our target, lets make the swap!
        if mutant[0] == '*' and mutant[1] == '*' or mutant[0] == mutant[1]:
            valid_mutations = codons[amino_acid]  # choose a silent mutation for whatever we are looking at
        else:
            valid_mutations = codons[mutant[1]]  # get a list of valid mutations
        replaceable = True
        for mutation in valid_mutations:
            if mutant[0] == mutant[1] or pam_case != 0:  # this is only true for PAM or seed-option for pam
                replaceable = True

                # identify the cases where we can't replace the PAM
                if distance_from_pam == 0:
                    if pam_case == 0 and not complement and (mutation.endswith(avoid_one) or mutation.endswith(avoid_two) or mutation.endswith(pam_indicator)):
                        replaceable = False # This will avoid making NAG/NGA mutations on PAM
                    if pam_case == 0 and complement and (mutation.startswith(pam_indicator) or mutation.startswith(avoid_one) or mutation.startswith(avoid_two)):
                        replaceable = False # Avoid mutating a PAM to a PAM (ex: CCT to CCA)
                    if pam_case == 1 and not down and (mutation[2] == pam_char_one or mutation[2] == pam_char_two):  # this would replace GG with GG
                        replaceable = False
                    if pam_case == 1 and down and (mutation[0] == pam_char_one or mutation[0] == pam_char_two):  # this would replace GG with GG
                        replaceable = False
                    if pam_case == 2 and not down:
                        temp = candidate_dna[first_amino_acid_loc + 3:first_amino_acid_loc + 5]
                        if temp == pam_indicator or temp == avoid_one or temp == avoid_two:
                            replaceable = False # This will not change any of the GGs so no point mutating it
                    if pam_case == 2 and down and (mutation.startswith(pam_indicator) or mutation.startswith(avoid_one) or mutation.startswith(avoid_two)):  # this would replace GG with GG
                        replaceable = False
                    if pam_case == 3 and not down and (mutation[0] == pam_char_one or mutation[0] == pam_char_two):  # this would replace CC with CC
                        replaceable = False
                    if pam_case == 3 and down and (mutation[2] == pam_char_one or mutation[2] == pam_char_two):  # this would replace CC with CC
                        replaceable = False
                    if pam_case == 4 and down and (mutation.endswith(pam_indicator) or mutation.endswith(avoid_one) or mutation.endswith(avoid_two)):  # this would replace CC with CC
                        replaceable = False
                    if pam_case == 4 and not down:
                        temp = candidate_dna[first_amino_acid_loc - 2:first_amino_acid_loc]
                        if temp == pam_indicator or temp == avoid_one or temp == avoid_two:
                            replaceable = False

                if mutation == candidate_dna[first_amino_acid_loc:first_amino_acid_loc + 3]:  # can't replace pam
                    replaceable = False

                # TODO: don't overwrite the original mutation

            elif mutation == candidate_dna[
                             first_amino_acid_loc:first_amino_acid_loc + 3]:  # This is what we already have, so it isn't a mutation
                continue
                # TODO:  If we are a pam, don't replace with G in the same place

                #     if possible, I would like to avoid mutating the PAM to NAG.
            #    2. If you can’t mutate the PAM then mutate at least one location in the “seed” region (10 bases upstream of the PAM). The closer the silent mutation is to the PAM the better it works. We may decide that we want to do two silent mutations if we find that one in the seed region isn’t enough to prevent re-cutting.

            # 3. Oftentimes the mutation we intend to make will be in the seed.
            #   a. If it is still possible to make a silent PAM mutation then that would be good (Although we are currently testing this and this parameter might change).
            #  b. If there is no way to make a silent PAM mutation and the mutation is more than 5 bases away from the PAM then the next best thing would be to make a silent mutation within the 5 bases upstream of the PAM.
            # c. If the mutation is within 5 bases from the PAM and you can’t make a silent PAM mutation, then I wouldn’t make any additional mutation.

            # This is from when I thought we couldn't introduce a new GG


            # if we are replacing a silent mutation outside of our range
            if (mutant[0] == '*' and mutant[1] == '*' or mutant[0] == mutant[1]) and distance_from_pam > try_until:
                replaceable = False

            # we are safe to make a swap here




            if replaceable:

                actual_mutation[0] = string_to_acid[candidate_dna[first_amino_acid_loc:first_amino_acid_loc + 3]]
                actual_mutation[1] = string_to_acid[mutation]
                if not config.USE_DEBUG_MUTATION:
                    candidate_dna = candidate_dna[:first_amino_acid_loc] + mutation + candidate_dna[first_amino_acid_loc + 3:]
                    decision_path += "Mutated " + amino_acid_str + " to " + mutation + ". "
                else:
                    candidate_dna = candidate_dna[:first_amino_acid_loc] + 'ZZZ' + candidate_dna[
                                                                                   first_amino_acid_loc + 3:]
                # keep going if we are doing silent mutations
                if (mutant[0] == '*' and mutant[1] == '*' or mutant[0] == mutant[1]) and (total_mutations + 1) < config.SILENT_MUTATION_MAXIMUM:
                    return perform_mutation(candidate_dna, first_amino_acid_loc + offset, pam_case, mutant,
                                            decision_path, mutation_location=mutation_location,
                                            distance_from_pam=distance_from_pam + 3, down=down, complement=complement, total_mutations=total_mutations+1)
                elif mutant[0] == mutant[1] and total_mutations < config.SILENT_MUTATION_MINIMUM:  # if we made a mutation, but it wasn't enough to reach the minimum, fail this one
                    return False, None, None, actual_mutation, decision_path
                else:
                    return True, candidate_dna, distance_from_pam, actual_mutation, decision_path
        if not replaceable and not down:
            if distance_from_pam > try_until:  # we ran out of space to look
                #if config.VERBOSE_EXECUTION:
                #    print('Could not find a replacement in the ' + config.SILENT_MUTATION_MODE)
                if total_mutations < config.SILENT_MUTATION_MINIMUM:  # We didn't find enough mutations
                    if config.VERBOSE_EXECUTION:
                        print('Could only find ' + str(total_mutations) + ' mutations and minimum is ' + str(config.SILENT_MUTATION_MINIMUM))
                    decision_path += 'Could only find ' + str(total_mutations) + ' mutations and minimum is ' + str(config.SILENT_MUTATION_MINIMUM) + '. '  # this doesn't actually show up because it is a failed mutation
                    return False, None, None, actual_mutation, decision_path
                else:
                    # we got enough mutations.  Carry on.
                    return True, candidate_dna, distance_from_pam, actual_mutation, decision_path
                #return False, None, None, actual_mutation, decision_path

            mutant[0] = '*'  # these two *s force any silent mutation
            mutant[1] = '*'
            decision_path += "Couldn't Mutate " + amino_acid_str + ". "
            return perform_mutation(candidate_dna, first_amino_acid_loc + offset, pam_case, mutant, decision_path, mutation_location=mutation_location,
                                        distance_from_pam=distance_from_pam + 3, down=down, complement=complement, total_mutations=total_mutations)

    if config.VERBOSE_EXECUTION:
        print('Mutant was not desirable')

    decision_path += "Couldn't Mutate " + amino_acid_str + ". "
    return False, None, None, actual_mutation, decision_path


# THis method will return a full dna string for each mutation as part of a MutationTracker type
# pam is the location of the first character of the pam
# mutant is key-val pair of mutant source to mutant destination [0] is key, [1] is value

def create_mutations(dna, pam, mutant, complement=False, only_once=False):
    global gs
    global guide_lib

    if pam == 1126:
        pass

    if pam == config.GENE_START_BUFFER + 8:    # was 1008
        pass
    # it seems like we are only looking at the 6 upstream and 4 downstream amino acids
    UPSTREAM = (config.UP_ACIDS) * 3
    DOWNSTREAM = config.DOWN_ACIDS * 3
    order = 0
    if complement:
        temp = UPSTREAM
        UPSTREAM = DOWNSTREAM
        DOWNSTREAM = temp
        order = 1

    # 1c) Take the 20 bases upstream of the NGG and that is the guide.
    guide = create_guides(dna, pam, complement)

    if (config.USE_RANK):
        rank = int(get_rank(guide)) # get the rank of the current guide

        # if rank is less than the threshold we skip this mutation
        # if rank is -1 (meaning the guide wasn't found in the file) we go ahead and do the mutation
        if (rank != -1 and rank < config.RANK_THRESHOLD):
            gs.failed_due_to_rank += 1
            return None
        else:
            gs.guides_used_from_rank_file += 1

    # If using the guide library, automatically reject any guide not present in the library
    if (config.USE_GUIDE_LIBRARY and not (is_guide_in_library(guide, guide_lib))):
        if config.VERBOSE_EXECUTION:
            print("Failed to find guide within guide library")
        gs.failed_due_to_guide_library += 1
        return None

    # 2a a. To make the donor, take 132 base pairs surrounding the mutations (either centered around both the main mutation and the PAM mutation, or if easier could just center all the donors for a given guide around the PAM).
    # PAM + 10 is the center of the guide + 66 on each side of it

    half = int(config.BP_LENGTH / 2)
    candidate_start = pam - half  # pam - UPSTREAM
    candidate_end = pam + half  # pam + 3 + DOWNSTREAM
    candidate_dna = dna[candidate_start:candidate_end]  # this is the 132 base pairs surrounding the middle of the guide

    # 2)  Find the amino acid you want to mutate


    # figure out the pam amino acid situation (does it split, and if so where)
    # 0 --> NGG
    # 1 --> XNG GXX
    # 2 --> XXN GGX

    pam_case = (pam - config.GENE_START_BUFFER) % 3
    if not complement:
        first_amino_acid_loc = pam - UPSTREAM - pam_case
    else:
        first_amino_acid_loc = pam - pam_case

    while first_amino_acid_loc < config.GENE_START_BUFFER:  # ignore acids outside the gene
        first_amino_acid_loc += 3
    while first_amino_acid_loc >= len(dna) - config.GENE_END_BUFFER:  # ignore acids outside the gene
        first_amino_acid_loc -= 3

    mutation_successful = False
    mutation_location = -1
    if (config.PRINT_MUTATION_CHECKS):
        print("Checking " + str(candidate_dna) + " for " + str(mutant[0]) + ".")

    successful_mutations = []
    candidate_dnas = []
    mutation_locations = []
    mutation_acids = []
    comments = []
    actual_mutation = []

    # ISSUE 25 - Only mutate in the PAM/Seed if "NULL" mutation is active.
    # If "NULL" mutation is active we skip this section that does the original mutations
    if (mutant[0] == 'NULL' and mutant[1] == 'NULL'):
        null_active = True;
        candidate_dnas.append(candidate_dna)
        actual_mutation = ['NULL', 'NULL']
    else:
        null_active = False

    if (not null_active):
        # Setting up for multiple potentials

        for ordering in range(2):
            if (ordering+order) % 2 == 0:  # only do upstream if we are still in the gene
                if complement:   # we have fo fix up the first amino acid location if we are on the reverse
                    first_amino_acid_loc = pam - UPSTREAM - pam_case

                for i in range(UPSTREAM - 3, -1, -3):  # check upstream, then check downstream
                    decision_path = ''
                    candidate_first_amino_acid_loc = first_amino_acid_loc + i - candidate_start
                    temp = candidate_first_amino_acid_loc + candidate_start

                    if not complement and (temp + 3 > pam):  # don't go into the pam (TODO:  I think this is true)
                        continue

                    if temp < config.GENE_START_BUFFER:
                        continue

                    if (config.PRINT_MUTATION_CHECKS):
                        print(candidate_dna[:candidate_first_amino_acid_loc] + " | " + candidate_dna[
                                                                                       candidate_first_amino_acid_loc:candidate_first_amino_acid_loc + 3] + " | " + candidate_dna[
                                                                                                                                                                    candidate_first_amino_acid_loc + 3:])
                    # 2)  Actually perform the mutation
                    mutation_successful, temp_candidate_dna, d_pam, actual_mutation, decision_path = perform_mutation(candidate_dna,
                                                                                      candidate_first_amino_acid_loc, 0, mutant, decision_path, complement=complement)
                    if config.TRACE_CANDIDATE_DNA_GENERATION:
                        print("Candidate DNA:")
                        print(temp_candidate_dna)
                    if mutation_successful:
                        # candidate_dna = temp_candidate_dna
                        # mutation_location = candidate_first_amino_acid_loc
                        # break

                        # Instead of setting them once, am now pushing them onto the list
                        candidate_dnas.append(temp_candidate_dna)
                        mutation_locations.append(candidate_first_amino_acid_loc)
                        mutation_acids.append(actual_mutation)
                        comments.append(decision_path)
                        if only_once:
                            break

            # if candidate_end > (len(dna) - config.GENE_END_BUFFER):   # Original version before my tweak
            # only do downstream if we are still in the gene
            if (ordering+order) % 2 == 1:

                if not complement:
                    first_amino_acid_loc = pam - pam_case

                for i in range(0, DOWNSTREAM, 3):  # check upstream, then check downstream
                    decision_path = ''
                    # convert first_amino_acid_loc from global dna to candidate dna
                    # candidate_first_amino_acid_loc = first_amino_acid_loc - candidate_start

                    candidate_first_amino_acid_loc = first_amino_acid_loc + i - candidate_start # this used to add UPSTREAM and I don't know why
                    temp = candidate_first_amino_acid_loc + candidate_start

                    if temp + 3 > (len(dna) - config.GENE_END_BUFFER):
                        continue

                    if (config.PRINT_MUTATION_CHECKS):
                        print(candidate_dna[:candidate_first_amino_acid_loc] + " | " + candidate_dna[
                                                                                       candidate_first_amino_acid_loc:candidate_first_amino_acid_loc + 3] + " | " + candidate_dna[
                                                                                                                                                                    candidate_first_amino_acid_loc + 3:])
                    # 2)  Actually perform the mutation
                    mutation_successful, temp_candidate_dna, d_pam, actual_mutation, decision_path = perform_mutation(candidate_dna,
                                                                                      candidate_first_amino_acid_loc, 0, mutant, decision_path, complement=complement)
                    if config.TRACE_CANDIDATE_DNA_GENERATION:
                        print("Reverse:")
                        print(temp_candidate_dna)
                    if mutation_successful:
                        candidate_dnas.append(temp_candidate_dna)
                        mutation_locations.append(candidate_first_amino_acid_loc)
                        mutation_acids.append(actual_mutation)
                        comments.append(decision_path)
                        if only_once:
                            break

        # if not mutation_successful:
        if len(candidate_dnas) == 0:
            if config.VERBOSE_EXECUTION:
                print('Failed to find a valid place to mutate ' + mutant[0] + ' into ' + mutant[1])
            gs.failed_due_to_mutate += 1
            return None

    for i in range(len(candidate_dnas)):
        candidate_dna = candidate_dnas[i]
        offset = 0
        decision_path = ''
        d_pam = None
        if (not null_active):
            mutation_location = mutation_locations[i]
            actual_mutation = mutation_acids[i]
            decision_path = comments[i]

        pam_loc_in_candidate = pam - candidate_start
        pam_string = dna[pam:pam + 3]
        mutation_successful = False
        make_pam_mutation = False
        decision_path += "Done Making OG Mutation. "

        # Decide if we need to make pam/seed mutation. OG mutation might have not disrupted the pam
        # or might have turned the pam into NGA/NAG (top strand) and NTC/NCT (bottom strand)
        if complement:
            if pam_case == 0:
                if (pam_string.startswith('CC')) or (pam_string.startswith('TC')) or (pam_string.startswith('CT')):
                    make_pam_mutation = True
                else:
                    make_pam_mutation = False
            elif pam_case == 1:
                pam_acid = candidate_dna[pam_loc_in_candidate - 1:pam_loc_in_candidate + 2]  # XCC
                if (pam_acid.endswith('CC')) or (pam_acid.endswith('TC')) or (pam_acid.endswith('CT')):
                    make_pam_mutation = True
                else:
                    make_pam_mutation = False
            else:
                pam_acid_one = candidate_dna[pam_loc_in_candidate - 2:pam_loc_in_candidate + 1]  # XXC
                pam_acid_two = candidate_dna[pam_loc_in_candidate + 1:pam_loc_in_candidate + 4]  # CNX
                if ((pam_acid_one.endswith('C') and pam_acid_two.startswith('C')) or
                        (pam_acid_one.endswith('C') and pam_acid_two.startswith('T')) or
                        (pam_acid_one.endswith('T') and pam_acid_two.startswith('C'))):
                    make_pam_mutation = True
                else:
                    make_pam_mutation = False
        else:
            if pam_case == 0:
                if (pam_string.endswith('GG')) or (pam_string.endswith('AG')) or (pam_string.endswith('GA')):
                    make_pam_mutation = True
                else:
                    make_pam_mutation = False
            elif pam_case == 1:
                pam_acid_one = candidate_dna[pam_loc_in_candidate - 1:pam_loc_in_candidate + 2] # XNG
                pam_acid_two = candidate_dna[pam_loc_in_candidate + 2:pam_loc_in_candidate + 5] # GXX
                if ((pam_acid_one.endswith('G') and pam_acid_two.startswith('G')) or
                        (pam_acid_one.endswith('G') and pam_acid_two.startswith('A')) or
                        (pam_acid_one.endswith('A') and pam_acid_two.startswith('G'))):
                    make_pam_mutation = True
                else:
                    make_pam_mutation = False
            else:
                pam_acid = candidate_dna[pam_loc_in_candidate + 1:pam_loc_in_candidate + 4]  # GGX
                if (pam_acid.startswith('GG')) or (pam_acid.startswith('AG')) or (pam_acid.startswith('GA')):
                    make_pam_mutation = True
                else:
                    make_pam_mutation = False

        if make_pam_mutation:
            decision_path += "Have to make Pam/seed mutation since OG mutation did not disrupt the PAM. "
            # 2)  mutate pam
            if pam_case == 0:  # we only have a single amino acid
            # go upstream only
                decision_path += "Pam is not split. "
                pam_acid = string_to_acid[pam_string]
                pam_mutant = [pam_acid, pam_acid]
                mutation_successful, candidate_dna, d_pam, pam_mutation, decision_path = perform_mutation(candidate_dna, pam_loc_in_candidate, 0,
                                                                             pam_mutant, decision_path,
                                                                             mutation_location=mutation_location, complement=complement)
                if config.TRACE_CANDIDATE_DNA_GENERATION:
                    print("PAM Candidate DNA:")
                    print(temp_candidate_dna)
                if not mutation_successful:
                    if config.VERBOSE_EXECUTION:
                        print('Failed to find a valid replacement for the pam')
                    gs.failed_due_to_pam += 1
                    return None
            else:
                if complement and pam_case == 1:
                    pam_case = 4
                elif complement and pam_case == 2:
                    pam_case = 3

                if pam_case == 2:  # the N is in one acid and the GG is in another so we can only replace down
                    # mutate down once (ONLY MUTATE ONE OF THE G)
                    # if failed then go upstream
                    # TODO:  make sure we don't overwrite the mutation by removing the pam
                    pam_string_up = dna[pam - 2:pam + 1]
                    pam_string_down = dna[pam + 1:pam + 4]
                    offset = 1
                    decision_path += "Pam is split, XXN GGX. "

                elif pam_case == 1:  # the NG is in one acid and the G is in another
                    # Check down stream first then go upstream
                    pam_string_up = dna[pam - 1:pam + 2]
                    pam_string_down = dna[pam + 2:pam + 5]
                    offset = 2
                    decision_path += "Pam is split, XNG GXX. "

                elif pam_case == 3:  # the NC is in second acid and the C is in first
                    # TODO:  make sure we don't overwrite the mutation by removing the pam
                    pam_string_up = dna[pam + 1:pam + 4]
                    pam_string_down = dna[pam - 2:pam + 1]
                    offset = -2
                    decision_path += "Pam is split, XXC CNX. "

                elif pam_case == 4:   # the N is in second acid and the CC is in first so we can only replace up
                    pam_string_up = dna[pam + 2:pam + 5]
                    pam_string_down = dna[pam - 1:pam + 2]
                    offset = -1
                    decision_path += "Pam is split, XCC NXX. "

                replaceable_pam = False
                pam_mutant_up = None
                pam_mutant_down = None
                if pam_string_up is not None and pam_string_up in string_to_acid:
                    pam_acid_up = string_to_acid[pam_string_up]
                    replaceable_pam = True
                    pam_mutant_up = [pam_acid_up, pam_acid_up]
                if pam_string_down is not None and pam_string_down in string_to_acid:
                    pam_acid_down = string_to_acid[pam_string_down]
                    replaceable_pam = True
                    pam_mutant_down = [pam_acid_down, pam_acid_down]
                if not replaceable_pam:
                    return None

                if pam_string_down is not None: # Try downstream first

                    mutation_successful, temp_candidate_dna, d_pam, pam_mutation, decision_path = perform_mutation(candidate_dna,
                        pam_loc_in_candidate + offset, pam_case, pam_mutant_down, decision_path,
                        mutation_location=mutation_location, complement=complement, down=True)

                    if config.TRACE_CANDIDATE_DNA_GENERATION:
                        print("PAM Candidate DNA 2:")
                        print(temp_candidate_dna)
                    if mutation_successful:
                        candidate_dna = temp_candidate_dna
                if not mutation_successful and pam_mutant_up is not None:  # try upstream if downstream didn't work
                    if config.TRACE_CANDIDATE_DNA_GENERATION:
                        print("PAM Candidate DNA 3:")
                        print(temp_candidate_dna)

                    if complement:
                        if pam_case == 3:
                            offset = 1
                        else:
                            offset = 2
                    else:
                        offset = -pam_case

                    mutation_successful, temp_candidate_dna, d_pam, pam_mutation, decision_path = perform_mutation(candidate_dna,
                        pam_loc_in_candidate + offset, pam_case, pam_mutant_up, decision_path,
                        mutation_location=mutation_location, complement=complement)
                    if mutation_successful:
                        candidate_dna = temp_candidate_dna
        else:
            decision_path += "OG mutation did disrupt the PAM. No Need to make pam/seed mutation. "
            mutation_successful = True  # we mutated the PAM already, so we didn't need to do another mutation

        if mutation_successful:

            if d_pam == None: # mutation was in the pam
                d_pam = 0

            if complement:
                pam_mutation_location = pam_loc_in_candidate + d_pam + offset
            else:
                pam_mutation_location = pam_loc_in_candidate - d_pam + offset

            d_pam = pam_mutation_location - pam_loc_in_candidate
            pam_loc = pam_loc_in_candidate
#            if complement:  # if we are on the reverse complement, invert it back before we add the other stuff
#                mutation_location = len(candidate_dna) - mutation_location - 3
#                pam_loc = len(candidate_dna) - pam_loc - 3

            size_of_first_sequence = len(config.first_sequence)
            size_of_second_sequence = len(config.second_sequence)
            total_to_add = size_of_first_sequence + size_of_second_sequence + config.GUIDE_LENGTH
            candidate_dna = insert_extra_sequence(candidate_dna, guide)
            # we just added 52 + 20 (guide) basepairs
            # guide pam mutation mutationloc dna
            # if we used a wildcard, replace it with what we actually mutated
            #
            # class MutationTracker:
            #     guide: str
            #     pam: int
            #     mutation: []
            #     mutation_loc: int
            #     dna: str
            #     complement: bool
            #     pam_location_in_gene: int
            #     distance_from_pam: int
            #     original_pam: str
            #     justification: str

            result = MutationTracker(guide, pam_loc + total_to_add, actual_mutation, mutation_location + total_to_add, candidate_dna, complement, pam,
                                     d_pam, pam_string, decision_path)

            # If using the guide library and only allowing one mutation per guide, automatically reject any guide not present in the library
            if (config.USE_GUIDE_LIBRARY and config.ONE_MUTATION_PER_GUIDE and not (
                    is_mutation_permitted(guide, result, guide_lib))):
                if config.VERBOSE_EXECUTION:
                    print("Mutation not permitted from guide library")
                gs.failed_due_to_guide_library += 1
            else:
                gs.succeeded += 1
                #prevent dups
                if result not in successful_mutations:
                    successful_mutations.append(result)
        else:
            if config.VERBOSE_EXECUTION:
                print('Mutation failed due to pam')  # TODO:  output why
            gs.failed_due_to_pam += 1


    if len(successful_mutations) == 0:
        return None

    else:
        return successful_mutations
    # TODO:  Track the decisions we made in this method so we can output them
    # TODO:  Fix naming of the file and frontmatter

# This method inverts the DNA so we can get the reverse complement strand
def invert_dna(dna):
    backwards_dna = dna[::-1]
    inv_dna = str()
    for base in backwards_dna:
        inv_dna += invert_mapping[base]
    return inv_dna


def write_results(frontmatter_list, results_list, dna_list, use_output_file=True):
    # use_output_file determines if the output in config will be used, or
    # if a title will automatically be generated from the frontmatter

    global gs

    if use_output_file:
        output_file = out_base;
    else:
        output_file = config.MULTI_OUTPUT_PREFIX + (str(frontmatter_list[0]).partition(' ')[0])[1:]

    if (config.PRINT_MUTATION_RESULTS):

        wb = Workbook()

        if config.PRINT_MUTATION_SUCCESS_COUNTS:
            print('\nfailed due to mutate: ' + str(gs.failed_due_to_mutate))
            print('failed due to pam: ' + str(gs.failed_due_to_pam))
            print('failed due to rank: ' + str(gs.failed_due_to_rank))
            if (config.USE_GUIDE_LIBRARY):
                print('failed due to guide library: ' + str(gs.failed_due_to_guide_library))
            print('succeeded: ' + str(gs.succeeded))
            print('\nGuides used from the rank file: ' + str(gs.guides_used_from_rank_file))
            print('Total guides in the file: ' + str(get_total_guides_from_file()))

        column_pos = 0

        sheet1 = wb.add_sheet('Mutation Results')
        sheet1.write(column_pos, 0, 'ID')
        sheet1.write(column_pos, 1, 'Mutation From')
        sheet1.write(column_pos, 2, 'Mutation To')
        sheet1.write(column_pos, 3, 'Mutation Offset from Start of Gene')
        sheet1.write(column_pos, 4, 'Reverse Complement')
        sheet1.write(column_pos, 5, 'Mutation Distance from Cut Site')
        sheet1.write(column_pos, 6, 'Original PAM')
        sheet1.write(column_pos, 7, 'Seed Mutation Distance From PAM')
        sheet1.write(column_pos, 8, 'Guides')
        sheet1.write(column_pos, 9, 'Result')
        sheet1.write(column_pos, 10, 'Comments')

        column_pos += 2

        extra_font = xlwt.easyfont('color_index gray50')
        mutation_font = xlwt.easyfont('color_index red')
        guide_font = xlwt.easyfont('color_index blue')
        pam_font = xlwt.easyfont('color_index green')
        dna_font = xlwt.easyfont('color_index black')
        pam_mut_font = xlwt.easyfont('color_index orange')

        first = config.first_sequence
        second = config.second_sequence
        third = config.third_sequence

        for g in range(len(results_list)):

            frontmatter = frontmatter_list[g]
            results = results_list[g]
            dna = dna_list[g]

            cur_id = (str(frontmatter).partition(' ')[0])[1:]
            i = 0
            for i, mutation in enumerate(results):

                if config.TRACE_CANDIDATE_DNA_GENERATION:
                    print("DNA during writing:")
                    print(mutation.dna)

                if (i == len(results) - 1):
                    sheet1.write(i + column_pos, 0, cur_id + "_" + config.KILL_MUTATION_ID_SUFFIX)
                else:
                    sheet1.write(i + column_pos, 0, cur_id + "_" + str(i))
                # We need the actual mutation location within the gene.
                sheet1.write(i + column_pos, 1, mutation.mutation[0] + str(int((mutation.pam_location_in_gene + (
                            mutation.mutation_loc - mutation.pam) - config.GENE_START_BUFFER) / 3) + 1))
                sheet1.write(i + column_pos, 2, mutation.mutation[1])
                sheet1.write(i + column_pos, 3, mutation.pam_location_in_gene + (
                            mutation.mutation_loc - mutation.pam) - config.GENE_START_BUFFER)
                sheet1.write(i + column_pos, 4, mutation.complement)

                if mutation.complement:
                    sheet1.write(i + column_pos, 5, str(abs(mutation.mutation_loc - (
                            mutation.pam + 6))))  # 3 bp downstream of pam + the length of the mutation (3 bp)
                else:
                    sheet1.write(i + column_pos, 5, str(abs(mutation.mutation_loc - (
                            mutation.pam - 6))))  # 3 bp upstream of pam + the length of the mutation (3 bp)

                sheet1.write(i + column_pos, 6, mutation.original_pam)
                sheet1.write(i + column_pos, 7, abs(mutation.distance_from_pam))
                sheet1.write(i + column_pos, 8, mutation.guide)
                sheet1.write(i + column_pos, 10, mutation.justification)

                if (mutation.complement):
                    pass

                # Get the extra fonts
                seg_first = (mutation.dna[0:len(first)], extra_font)
                seg_guide = (mutation.dna[len(first):len(first) + config.GUIDE_LENGTH], guide_font)
                seg_second = (mutation.dna[len(first) + config.GUIDE_LENGTH:len(first) + config.GUIDE_LENGTH + len(second)],
                    extra_font)
                seg_third = (mutation.dna[len(mutation.dna) - len(third):], extra_font)
                blank = ('', pam_font)

                # If mutation is NULL, we don't need to show og
                if mutation.mutation[0] == 'NULL':
                    og_mutation = blank
                else:
                    og_mutation = (mutation.dna[mutation.mutation_loc: mutation.mutation_loc + 3], mutation_font)

                mode = mutation.pam_location_in_gene % 3  # this is which "mode" we were in, it gives us the offset from pam start
                seed_mutation = None

                start_of_pam_mutation = mutation.pam + mutation.distance_from_pam
                seed_mutation = (mutation.dna[start_of_pam_mutation:start_of_pam_mutation + 3], pam_mut_font)
                og_mutation_distance_from_pam = abs(mutation.pam - mutation.mutation_loc)
                distance_from_pam_positive = abs(mutation.distance_from_pam)

                # Cases where we don't need to show the seed mutation
                if mutation.distance_from_pam == 0 and (og_mutation_distance_from_pam == 1 or og_mutation_distance_from_pam == 2):
                    seed_mutation = blank
                elif start_of_pam_mutation == mutation.mutation_loc:
                    seed_mutation = blank

                # identify how much, if any of the pam we mutated
                # 1 -> NGG
                # 0 -> XXN GGX
                # 2 -> XNG GXX
                if mode == 1 and (distance_from_pam_positive == 0 or og_mutation_distance_from_pam == 0):  # we mutated the entire pam
                    seg_pam = blank  # we don't need to show the pam
                elif mode == 0 and ((og_mutation_distance_from_pam == 2 and distance_from_pam_positive == 1) or (og_mutation_distance_from_pam == 1 and distance_from_pam_positive == 2)):
                    seg_pam = blank # we don't need to show the pam
                elif mode == 2 and ((og_mutation_distance_from_pam == 1 and distance_from_pam_positive == 2) or (og_mutation_distance_from_pam == 2 and distance_from_pam_positive == 1)):
                    seg_pam = blank # we don't need to show the pam
                elif mode == 0 and (distance_from_pam_positive == 2 or og_mutation_distance_from_pam == 2):  # we mutated the first char
                    seg_pam = (mutation.dna[mutation.pam + 1: mutation.pam + 3], pam_font)
                elif mode == 0 and (distance_from_pam_positive == 1 or og_mutation_distance_from_pam == 1):  # we mutated the last two characters
                    seg_pam = (mutation.dna[mutation.pam: mutation.pam + 1], pam_font)
                elif mode == 2 and (distance_from_pam_positive == 2 or og_mutation_distance_from_pam == 2):  # we mutated the last char
                    seg_pam = (mutation.dna[mutation.pam: mutation.pam + 2], pam_font)
                elif mode == 2 and (distance_from_pam_positive == 1 or og_mutation_distance_from_pam == 1):  # we mutated the first two chars
                    seg_pam = (mutation.dna[mutation.pam + 2: mutation.pam + 3], pam_font)
                else:  # we mutated the seed #DONE
                    seg_pam = (mutation.dna[mutation.pam: mutation.pam + 3], pam_font)

                if mutation.mutation_loc <= mutation.pam: # upstream mutation
                    if start_of_pam_mutation > mutation.pam and start_of_pam_mutation > mutation.mutation_loc: # seed mutation is after the pam and after og
                        seg_dna1 = (mutation.dna[len(first) + config.GUIDE_LENGTH + len(second):mutation.mutation_loc], dna_font)
                        seg_dna2 = (mutation.dna[mutation.mutation_loc + 3:mutation.pam], dna_font)
                        seg_dna3 = (mutation.dna[mutation.pam + 3: start_of_pam_mutation], dna_font)
                        seg_dna4 = (mutation.dna[start_of_pam_mutation + 3:len(mutation.dna) - len(third)], dna_font)
                        sheet1.write_rich_text(i + column_pos, 9, (
                            seg_first, seg_guide, seg_second, seg_dna1, og_mutation, seg_dna2, seg_pam,
                            seg_dna3, seed_mutation, seg_dna4, seg_third))
                    elif start_of_pam_mutation > mutation.mutation_loc and start_of_pam_mutation <= mutation.pam: # seed mutation is between og and pam
                        seg_dna1 = (mutation.dna[len(first) + config.GUIDE_LENGTH + len(second):mutation.mutation_loc], dna_font)
                        seg_dna2 = (mutation.dna[mutation.mutation_loc + 3:start_of_pam_mutation], dna_font)
                        seg_dna3 = (mutation.dna[start_of_pam_mutation + 3: mutation.pam], dna_font)
                        seg_dna4 = (mutation.dna[mutation.pam + 3 :len(mutation.dna) - len(third)], dna_font)
                        sheet1.write_rich_text(i + column_pos, 9, (
                            seg_first, seg_guide, seg_second, seg_dna1, og_mutation, seg_dna2, seed_mutation,
                            seg_dna3, seg_pam, seg_dna4, seg_third))
                    elif start_of_pam_mutation < mutation.mutation_loc and start_of_pam_mutation < mutation.pam: # seed mutation is before the og and pam
                        seg_dna1 = (mutation.dna[len(first) + config.GUIDE_LENGTH + len(second):start_of_pam_mutation], dna_font)
                        seg_dna2 = (mutation.dna[start_of_pam_mutation + 3:mutation.mutation_loc], dna_font)
                        seg_dna3 = (mutation.dna[mutation.mutation_loc + 3:mutation.pam], dna_font)
                        seg_dna4 = (mutation.dna[mutation.pam + 3:len(mutation.dna) - len(third)], dna_font)
                        sheet1.write_rich_text(i + column_pos, 9, (
                            seg_first, seg_guide, seg_second, seg_dna1, seed_mutation, seg_dna2, og_mutation,
                            seg_dna3, seg_pam, seg_dna4, seg_third))
                else: # downstream mutation
                    if (start_of_pam_mutation <= mutation.pam): # seed mutation is before the pam or completely overlap
                        seg_dna1 = (mutation.dna[len(first) + config.GUIDE_LENGTH + len(second):start_of_pam_mutation], dna_font)
                        seg_dna2 = (mutation.dna[start_of_pam_mutation + 3:mutation.pam], dna_font)
                        seg_dna3 = (mutation.dna[mutation.pam + 3:mutation.mutation_loc], dna_font)
                        seg_dna4 = (mutation.dna[mutation.mutation_loc + 3:len(mutation.dna) - len(third)], dna_font)
                        sheet1.write_rich_text(i + column_pos, 9, (
                            seg_first, seg_guide, seg_second, seg_dna1, seed_mutation, seg_dna2, seg_pam,
                            seg_dna3, og_mutation, seg_dna4, seg_third))
                    elif start_of_pam_mutation > mutation.pam and start_of_pam_mutation < mutation.mutation_loc: # seed mutation is after the pam but before og
                        seg_dna1 = (mutation.dna[len(first) + config.GUIDE_LENGTH + len(second):mutation.pam], dna_font)
                        seg_dna2 = (mutation.dna[mutation.pam + 3:start_of_pam_mutation], dna_font)
                        seg_dna3 = (mutation.dna[start_of_pam_mutation + 3:mutation.mutation_loc], dna_font)
                        seg_dna4 = (mutation.dna[mutation.mutation_loc + 3:len(mutation.dna) - len(third)], dna_font)
                        sheet1.write_rich_text(i + column_pos, 9, (
                            seg_first, seg_guide, seg_second, seg_dna1, seg_pam, seg_dna2, seed_mutation,
                            seg_dna3, og_mutation, seg_dna4, seg_third))
                    elif start_of_pam_mutation > mutation.pam and start_of_pam_mutation >= mutation.mutation_loc: # seed mutation is after the pam and after og
                        seg_dna1 = (mutation.dna[len(first) + config.GUIDE_LENGTH + len(second):mutation.pam], dna_font)
                        seg_dna2 = (mutation.dna[mutation.pam + 3:mutation.mutation_loc], dna_font)
                        seg_dna3 = (mutation.dna[mutation.mutation_loc + 3:start_of_pam_mutation], dna_font)
                        seg_dna4 = (mutation.dna[start_of_pam_mutation + 3:len(mutation.dna) - len(third)], dna_font)
                        sheet1.write_rich_text(i + column_pos, 9, (
                            seg_first, seg_guide, seg_second, seg_dna1, seg_pam, seg_dna2, og_mutation,
                            seg_dna3, seed_mutation, seg_dna4, seg_third))

                if (mutation.complement):
                    mutation.dna = invert_dna(mutation.dna)

            column_pos += i + 2

            sheet1.write(column_pos, 0, "Original Sequence")
            column_pos += 1
            sheet1.write(column_pos, 0, dna)
            column_pos += 1
            sheet1.write(column_pos, 0, "Original Gene")
            column_pos += 1
            sheet1.write(column_pos, 0, dna[config.GENE_START_BUFFER:-config.GENE_END_BUFFER])
            # NOTE // Please check in and make sure that crops correctly
            column_pos += 2
            # Print both full sequence, as well as just the gene - the thousand surrounding pairs

            # Saves the file
            wb.save(output_file + '.xls')

    if (config.PRINT_GUIDE_LIBRARY):

        wb = Workbook()

        column_pos = 0

        sheet2 = wb.add_sheet('Guide Library')
        sheet2.write(column_pos, 0, 'GENE ID')
        sheet2.write(column_pos, 1, 'GUIDE ID')
        sheet2.write(column_pos, 2, 'GUIDE')
        sheet2.write(column_pos, 3, 'MUTATION POSITION')
        sheet2.write(column_pos, 4, 'INVERSE COMPLIMENT')
        sheet2.write(column_pos, 5, 'GUIDE PRIORITY')
        sheet2.write(column_pos, 6, 'INVERSE COMPLIMENT PRIORITY')
        sheet2.write(column_pos, 7, 'GUIDE WITH HEADER')
        sheet2.write(column_pos, 8, 'INVERSE COMPLIMENT WITH HEADER')

        column_pos += 2

        extra_font = xlwt.easyfont('color_index gray50')
        guide_font = xlwt.easyfont('color_index blue')

        for g in range(len(results_list)):

            frontmatter = frontmatter_list[g]
            results = results_list[g]
            dna = dna_list[g]

            cur_id = (str(frontmatter).partition(' ')[0])[1:]

            guides = []
            inv_guides = []

            i = 0
            for i, mutation in enumerate(results):
                sheet2.write(i + column_pos, 0, cur_id)
                sheet2.write(i + column_pos, 1, cur_id + "_" + str(i))
                sheet2.write(i + column_pos, 3, mutation.mutation_loc)

                # Temp priority system for testing
                sheet2.write(i + column_pos, 5, i)

                guide = (mutation.dna[len(first):len(first) + config.GUIDE_LENGTH], guide_font)
                inv_guide = (invert_dna(mutation.dna[len(first):len(first) + config.GUIDE_LENGTH]), guide_font)

                prefix = config.GUIDE_LIBRARY_STRAND_PREFIX
                inv_prefix = config.GUIDE_LIBRARY_INVERSE_PREFIX
                _blank = ("", dna_font)
                if mutation.mutation_loc < mutation.pam:  # upstream mutation
                    _prefix = (prefix, dna_font)
                    sheet2.write_rich_text(i + column_pos, 2, (_blank, guide))
                    sheet2.write_rich_text(i + column_pos, 7, (_prefix, guide))
                    guides.append(_prefix[0] + guide[0])

                    _prefix = (inv_prefix, dna_font)
                    sheet2.write_rich_text(i + column_pos, 4, (_blank, inv_guide))
                    sheet2.write_rich_text(i + column_pos, 8, (_prefix, inv_guide))
                    inv_guides.append(_prefix[0] + inv_guide[0])
                else:  # downstream mutation
                    _prefix = (inv_prefix, dna_font)
                    sheet2.write_rich_text(i + column_pos, 2, (_blank, guide))
                    sheet2.write_rich_text(i + column_pos, 7, (_prefix, guide))
                    inv_guides.append(_prefix[0] + guide[0])

                    _prefix = (prefix, dna_font)
                    sheet2.write_rich_text(i + column_pos, 4, (_blank, inv_guide))
                    sheet2.write_rich_text(i + column_pos, 8, (_prefix, inv_guide))
                    guides.append(_prefix[0] + inv_guide[0])

            column_pos += i + 2

        # END OF FILE LINE FOR READER
        sheet2.write(column_pos, 0, config.GUIDE_LIBRARY_EOF)

        # PROBABLY REDO THIS
        # Dont use json, figure out how to read the excel back in directly
        # Formatting the output data for readability

        wb.save(config.GUIDE_LIBRARY_OUTPUT_FILE + '.xls')

        # out_data = ["Guides", guides, "Guides inverted", inv_guides]
        # with open(config.GUIDE_LIBRARY_OUTPUT_FILE+".json", "w") as file:
        #    json.dump(out_data, file)

# Checks if no arguments are given, will set input and output files to the default
# if none are present.
# Prints warnings if doing so, or if no defaults are configured.
def confirm_input_args():
    global in_files
    global out_base

    # Checks if argument inputs are present, otherwise defaults to what is in config

    # Checking in_files
    if (in_files == None):
        print("No input file detected, using default.")
        in_files = config.DEFAULT_IN_FILES
        if (in_files == None):
            print("No default input file present.")

    # Checking out_file
    if (out_base == None):
        print("No output file detected, using default.")
        out_base = config.DEFAULT_OUT_FILE
        if (config.DEFAULT_OUT_FILE == None):
            print("No default output file present.")


def get_all_mutations(_dna_locs, _inv_dna_locs, _dna, _inv_dna, only_once=False):
    _only_once = only_once

    mutations_output = []
    mutation_count = 0
    if not config.DEBUG_INVERSE:

        for loc in _dna_locs:
            # I don't think we need to actually find the guides here since it is just pam - 20
            # guides = create_guides(dna, loc)
            # for g in guides:   # for each guide do each mutation

            for m in config.mutations_to_attempt.items():
                mutated_dna = create_mutations(_dna, loc, list(m), only_once=_only_once)
                if mutated_dna is not None:
                    for md in mutated_dna:
                        mutations_output.append(md)
                        mutation_count += 1
            if config.KILL_MODE and mutation_count >= config.MAX_PAMS:  # Talia only wants 3 pams to get a stop codon
                break


    # class MutationTracker:
    #    guide: int
    #    pam: int
    #    mutation: []
    #    mutation_loc: int
    #    dna: str

    for loc in _inv_dna_locs:
        for m in config.mutations_to_attempt.items():
            mutated_dna = create_mutations(_dna, loc, list(m), only_once=_only_once,
                                           complement=True)  # returns a mutation tracker
            if mutated_dna is not None:
                for md in mutated_dna:
                    mutation_count += 1
                    mutations_output.append(md)

    return mutations_output


# Quickly moved the execution to its own function, for later use with a seperate driver (potentially, easy enough to revert if not)

def execute_program():
    global guide_lib

    # 1)  Choose the gene to mutate
    frontmatter, dna_list, guide_lib = get_dna()
    # At this stage, the dna is the full dna, buffer still included.

    if config.USE_GUIDE_LIBRARY:
        guide_lib = filter_guides(guide_lib)
    combined_mutation_page = []

    for i in range(0, len(dna_list)):
        dna = dna_list[i]
        inv_dna_full = invert_dna(dna)
        # 1b. Search for NGG (where N is any base, A, T, C, or G) (aka the PAM)
        #  This function finds all the pams

        pams = get_locations(dna)  # Returns pams from both regular strand and inverse compliment
        dna_locs = pams[0]
        inv_dna_locs = pams[1]

        if (len(dna_locs) > 0):  # Make sure there actually are pams to use
            candidate_start = int(dna_locs[0]) - 10 - 66  # pam - UPSTREAM
            candidate_end = int(dna_locs[0]) - 10 + 66  # pam + 3 + DOWNSTREAM
            candidate_dna = dna[candidate_start:candidate_end]
        if config.KILL_MODE:
            all_mutations = get_all_mutations(dna_locs, [], dna, [], only_once=True)  # we don't need to do inverse because we can get 3 on the ,main side
        else:
            all_mutations = get_all_mutations(dna_locs, inv_dna_locs, dna, inv_dna_full)

        if not config.KILL_MODE:
            # NOTE // Kill guide is inserted as the very last one
            if (len(all_mutations) > 0):
                all_mutations.append(create_kill_guide(all_mutations[0]))

        combined_mutation_page.append(all_mutations)
        # at this point, we have everything we need to output the results
        if not config.OUTPUT_TO_ONE_FILE:
            write_results([frontmatter[i]], [all_mutations], [dna], False)

    if config.OUTPUT_TO_ONE_FILE:
        write_results(frontmatter, combined_mutation_page, dna_list)

# This method gets the total number of guides in a file
def get_total_guides_from_file():
    workbook = xlrd.open_workbook(config.RANK_FILE + ".xls")
    worksheet = workbook.sheet_by_index(0)

    # Get the index of guide column first
    guide_index = -1
    reading = True
    i = 0
    while reading:
        column_name = worksheet.cell_value(0, i)

        if column_name != "":
            if column_name == config.GUIDE_COLUMN_IN_RANK_FILE:
                guide_index = i

            i += 1

            # we found the indicies, stop searching
            if guide_index != -1:
                reading = False

    # now count the guides
    reading = True
    i = 1
    while (reading):
        guide_from_file = worksheet.cell_value(i, guide_index)

        # We're at the end
        if (guide_from_file == "END"):
            return i - 1

        i += 1

# This method takes in a guide and finds that guide's rank in the rank file
# and returns it. If the guide is not in the file, returns -1.
def get_rank(guide):
    workbook = xlrd.open_workbook(config.RANK_FILE + ".xls")
    worksheet = workbook.sheet_by_index(0)

    # rank files can be different format so we want to find the index values of "Guide" and "Rank" column
    # get the index values of guide and rank
    guide_index = -1
    rank_index = -1
    reading = True
    i = 0
    while reading:
        column_name = worksheet.cell_value(0, i)

        if column_name != "":
            if column_name == config.GUIDE_COLUMN_IN_RANK_FILE:
                guide_index = i
            elif column_name == config.RANK_COLUMN_IN_RANK_FILE:
                rank_index = i

            i += 1

            # we found the indicies, stop searching
            if guide_index != -1 and rank_index != -1:
                reading = False

    # now look for the guide
    reading = True
    i = 1
    while (reading):
        guide_from_file = worksheet.cell_value(i, guide_index)

        if (guide_from_file != ""):
            if (guide_from_file == guide):
                return worksheet.cell_value(i, rank_index) # return the rank

            if (guide_from_file == config.GUIDE_LIBRARY_EOF):
                return -1 # guide wasn't in the file

        i += 1

def test_execution():
    global guide_lib

    # 1)  Choose the gene to mutate
    frontmatter, dna_list, guide_lib = get_dna()
    # At this stage, the dna is the full dna, buffer still included.

    if config.USE_GUIDE_LIBRARY:
        guide_lib = filter_guides(guide_lib)
    combined_mutation_page = []

    dna = dna_list[0]
    inv_dna_full = invert_dna(dna)
    print(inv_dna_full)
    # 1b. Search for NGG (where N is any base, A, T, C, or G) (aka the PAM)
    #  This function finds all the pams

    pams = get_locations(dna)  # Returns pams from both regular strand and inverse compliment
    dna_locs = pams[0][0]
    inv_dna_locs = pams[1][0]

    candidate_start = dna_locs - 10 - 66  # pam - UPSTREAM
    candidate_end = dna_locs - 10 + 66  # pam + 3 + DOWNSTREAM
    candidate_dna = dna[candidate_start:candidate_end]
    print(dna_locs)
    all_mutations = get_all_mutations([dna_locs], [inv_dna_locs], dna, inv_dna_full, only_once=False)

    # at this point, we have everything we need to output the results
    write_results([frontmatter[0]], [all_mutations], [dna], False)

if config.RUN_IN_EXECUTION_TESTING_MODE:
    test_execution()
else:
    execute_program()
