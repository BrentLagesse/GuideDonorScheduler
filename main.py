import argparse
import re
from dataclasses import dataclass
import sys
import json

# xlrd, xlutils and xlwt modules need to be installed.  
# Can be done via pip install <module>

import xlwt
from xlwt import Workbook
import xlrd
import xlutils

# Other modules
import config


@dataclass
class MutationTracker:
    guide: int
    pam: int
    mutation: []
    mutation_loc: int
    dna: str
    complement: bool
    pam_location_in_gene: int
    distance_from_pam: int
    original_pam: str


@dataclass
class GlobalStats:
    failed_due_to_mutate: int
    failed_due_to_pam: int
    failed_due_to_guide_library: int
    succeeded: int


@dataclass
class PrebuiltGuide:
    gene: str
    id: str
    guide: str
    mutation_loc: int
    priority: float
    compliment: bool


gs = GlobalStats(0, 0, 0, 0)
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

invert_mapping = dict()
invert_mapping['T'] = 'A'
invert_mapping['A'] = 'T'
invert_mapping['G'] = 'C'
invert_mapping['C'] = 'G'
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
        for line in in_files:
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
    # AMINO_ACID_IGNORE = 1 * 3 # // Ignores the first amino acid in the sequence
    AMINO_ACID_IGNORE = 3 - config.GENE_START_BUFFER % 3  # // Ignores the first amino acid in the sequence that has a piece in the buffer
    gene_only = dna[config.GENE_START_BUFFER + AMINO_ACID_IGNORE:len(dna) - config.GENE_END_BUFFER]
    regular = gene_only[AMINO_ACID_IGNORE:]
    gg_locs = [loc.start() + config.GENE_START_BUFFER - 1 + AMINO_ACID_IGNORE for loc in
               re.finditer('(?=GG)', gene_only)]  # minus one accounts for the N of NGG
    inverted = invert_dna(gene_only[AMINO_ACID_IGNORE:])
    cc_locs = [loc.start() + config.GENE_START_BUFFER - 1 for loc in
               re.finditer('(?=GG)', inverted)]  # Double check if this also needs a -1?
    # Make sure its using nGG pams for reverse compliment, look into this and see why nCCs are used

    return [gg_locs, cc_locs]


# Returns the 20 base-pairs before the PAM location
def create_guides(dna, loc):
    return dna[loc - 20:loc]


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
    m2 = MutationTracker(m.guide, m.pam, m.mutation, m.mutation_loc, m.dna, m.complement, m.pam_location_in_gene, 0, '')

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


def perform_mutation(candidate_dna, first_amino_acid_loc, pam_case, mutant, keep_trying=False, distance_from_pam=0,
                     mutation_location=-1):
    # if first_amino_acid_loc > 76, we are replacing downstream
    if first_amino_acid_loc == mutation_location:  # we would be replacing the mutation
        if distance_from_pam <= 5:  # don't mutate
            if config.VERBOSE_EXECUTION:
                print('we did not mutate the pam because the mutation was withing 5 base pairs of pam')
            return True, candidate_dna, distance_from_pam
        return False, None, None
    amino_acid_str = candidate_dna[first_amino_acid_loc: first_amino_acid_loc + 3]
    if config.PRINT_MUTATION_CHECKS:
        print(amino_acid_str)
    if amino_acid_str in string_to_acid:  # if this is something that isn't an amino acid, just quit
        amino_acid = string_to_acid[amino_acid_str]
    else:
        if amino_acid_str in config.stop:
            if config.VERBOSE_EXECUTION:
                print('Ran into a stop codon: ' + amino_acid_str)
            return False, None, None
        else:
            if config.VERBOSE_EXECUTION:
                print('Somehow we ran into something that was not an amino acid: ' + amino_acid_str)
            return False, None, None

    if config.PRINT_MUTATION_CHECKS:
        print("Currently checking " + str(amino_acid) + " for " + str(mutant[0]) + ".")

    if mutant[0] == amino_acid or mutant[0] == '*':  # we found our target, lets make the swap!
        if mutant[0] == '*' and mutant[1] == '*':
            valid_mutations = codons[amino_acid]  # choose a silent mutation for whatever we are looking at
        else:
            valid_mutations = codons[mutant[1]]  # get a list of valid mutations
        replaceable = True
        for mutation in valid_mutations:
            if mutant[0] == mutant[1] or pam_case != 0:  # this is only true for PAM or seed-option for pam
                replaceable = True
                # identify the cases where we can't replace the PAM
                if pam_case == 1 and mutation[0] == 'G' and mutation[1] == 'G':  # this would replace GG with GG
                    replaceable = False

                # if mutation[2] == 'G' and amino_acid_str[2] == 'G':  # this would introduce a GG on the back end

                if mutation == candidate_dna[first_amino_acid_loc:first_amino_acid_loc + 3]:  # can't replace pam
                    replaceable = False

                # TODO: don't overwrite the original mutation

                # we couldn't replace it, so lets try option #2, the seed (10 bases upstream)


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

            # we are safe to make a swap here
            if replaceable:
                if not config.USE_DEBUG_MUTATION:
                    candidate_dna = candidate_dna[:first_amino_acid_loc] + mutation + candidate_dna[
                                                                                      first_amino_acid_loc + 3:]
                else:
                    candidate_dna = candidate_dna[:first_amino_acid_loc] + 'ZZZ' + candidate_dna[
                                                                                   first_amino_acid_loc + 3:]
                return True, candidate_dna, distance_from_pam
        if not replaceable:
            if distance_from_pam > 8:  # Couldn't find anything in the seed, so quit -- we would be 11 from pam on next run
                if config.VERBOSE_EXECUTION:
                    print('Could not find a replacement in the seed')
                return False, None, None
            mutant[0] = '*'  # these two *s force any silent mutation
            mutant[1] = '*'
            return perform_mutation(candidate_dna, first_amino_acid_loc - 3, 3, mutant,
                                    distance_from_pam=distance_from_pam + 3)

    if config.VERBOSE_EXECUTION:
        print('Mutant was not desireable')
    return False, None, None


# THis method will return a full dna string for each mutation as part of a MutationTracker type
# pam is the location of the first character of the pam
# mutant is key-val pair of mutant source to mutant destination [0] is key, [1] is value

# only_once is used for debug, only runs through one mutation
def create_mutations(dna, pam, mutant, complement=False, only_once=False):
    global gs
    global guide_lib
    # it seems like we are only looking at the 6 upstream and 4 downstream amino acids
    UPSTREAM = (config.UP_ACIDS + 2) * 3
    DOWNSTREAM = config.DOWN_ACIDS * 3

    # 1c) Take the 20 bases upstream of the NGG and that is the guide.
    guide = create_guides(dna, pam)

    # If using the guide library, automatically reject any guide not present in the library
    if (config.USE_GUIDE_LIBRARY and not (is_guide_in_library(guide, guide_lib))):
        if config.VERBOSE_EXECUTION:
            print("Failed to find guide within guide library")
        gs.failed_due_to_guide_library += 1
        return None

    # If we are on the reverse compliment, invert the guide // NOTE this might not work yet
    #    if (complement):
    #        guide = invert_dna(guide)

    # introduce mutation

    # 2a a. To make the donor, take 132 base pairs surrounding the mutations (either centered around both the main mutation and the PAM mutation, or if easier could just center all the donors for a given guide around the PAM).
    # PAM + 10 is the center of the guide + 66 on each side of it
    candidate_start = pam - 10 - 66  # pam - UPSTREAM
    candidate_end = pam - 10 + 66  # pam + 3 + DOWNSTREAM
    candidate_dna = dna[candidate_start:candidate_end]  # this is the 132 base pairs surrounding the middle of the guide

    # 2)  Find the amino acid you want to mutate

    first_amino_acid_loc = int()
    origin = first_amino_acid_loc = pam - UPSTREAM
    for i in range(0, 3):  # We want to start on the first amino acid that is within our upstream range
        if (pam - config.GENE_START_BUFFER + i) % 3 == 0:
            first_amino_acid_loc = pam - UPSTREAM + i
    while first_amino_acid_loc < config.GENE_START_BUFFER:  # ignore acids outside the gene
        first_amino_acid_loc += 3

    pam_modifier = origin - first_amino_acid_loc
    # pam -= pam_modifier

    # print(dna[pam-20:pam+20])
    # print(dna[pam-20:pam] + " " + dna[pam:pam+3] + " " + dna[pam+3:pam+20] )

    # first_amino_acid_loc += 3 # Adding 3 to ignore the start codon

    # candidate_dna = dna[first_amino_acid_loc:candidate_end]   #grab starting from the first full amino acid
    mutation_successful = False
    mutation_location = -1
    if (config.PRINT_MUTATION_CHECKS):
        print("Checking " + str(candidate_dna) + " for " + str(mutant[0]) + ".")

    # Setting up for multiple potentials

    candidate_dnas = []
    mutation_locations = []

    if first_amino_acid_loc >= config.GENE_START_BUFFER and first_amino_acid_loc + 6 < pam:  # only do upstream if we are still in the gene
        #        for i in range(UPSTREAM - 3, -1, -3):    # check upstream, then check downstream
        for i in range(UPSTREAM - 3, -1, -3):  # check upstream, then check downstream
            if i + first_amino_acid_loc + 3 >= pam:  # don't go into the pam (TODO:  I think this is true)
                continue
            # convert first_amino_acid_loc from global dna to candidate dna
            # candidate_first_amino_acid_loc = first_amino_acid_loc - candidate_start
            candidate_first_amino_acid_loc = first_amino_acid_loc + i - candidate_start
            if (config.PRINT_MUTATION_CHECKS):
                print(candidate_dna[:candidate_first_amino_acid_loc] + " | " + candidate_dna[
                                                                               candidate_first_amino_acid_loc:candidate_first_amino_acid_loc + 3] + " | " + candidate_dna[
                                                                                                                                                            candidate_first_amino_acid_loc + 3:])
            # 2)  Actually perform the mutation
            mutation_successful, temp_candidate_dna, d_pam = perform_mutation(candidate_dna,
                                                                              candidate_first_amino_acid_loc, 0, mutant)
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
                if only_once:
                    break

    # if candidate_end > (len(dna) - config.GENE_END_BUFFER):   # Original version before my tweak
    if candidate_end < (len(dna) - config.GENE_END_BUFFER):  # only do downstream if we are still in the gene
        for i in range(0, DOWNSTREAM, 3):  # check upstream, then check downstream
            # convert first_amino_acid_loc from global dna to candidate dna
            # candidate_first_amino_acid_loc = first_amino_acid_loc - candidate_start
            candidate_first_amino_acid_loc = first_amino_acid_loc + i - candidate_start  # this used to add UPSTREAM and I don't know why
            if (config.PRINT_MUTATION_CHECKS):
                print(candidate_dna[:candidate_first_amino_acid_loc] + " | " + candidate_dna[
                                                                               candidate_first_amino_acid_loc:candidate_first_amino_acid_loc + 3] + " | " + candidate_dna[
                                                                                                                                                            candidate_first_amino_acid_loc + 3:])
            # 2)  Actually perform the mutation
            mutation_successful, temp_candidate_dna, d_pam = perform_mutation(candidate_dna,
                                                                              candidate_first_amino_acid_loc, 0, mutant)
            if config.TRACE_CANDIDATE_DNA_GENERATION:
                print("Reverse:")
                print(temp_candidate_dna)
            if mutation_successful:
                candidate_dnas.append(temp_candidate_dna)
                mutation_locations.append(candidate_first_amino_acid_loc)
                if only_once:
                    break

    # if not mutation_successful:
    if len(candidate_dnas) == 0:
        if config.VERBOSE_EXECUTION:
            print('Failed to find a valid place to mutate ' + mutant[0] + ' into ' + mutant[1])
        gs.failed_due_to_mutate += 1
        return None

    # if we wrote over the pam already, we are fine, I think

    successful_mutations = []

    for i in range(len(candidate_dnas)):
        candidate_dna = candidate_dnas[i]
        mutation_location = mutation_locations[i]

        pam_loc_in_candidate = 76  # this is always true
        # pam_string = candidate_dna[pam_loc_in_candidate:pam_loc_in_candidate + 3]
        pam_string = dna[pam:pam + 3]
        mutation_successful = False
        if 'GG' in (candidate_dna[pam_loc_in_candidate:pam_loc_in_candidate + 3]):

            # 2)  mutate pam
            # figure out the pam amino acid situation (does it split, and if so where)
            pam_case = (pam - config.GENE_START_BUFFER) % 3

            if pam_case == 0:  # we only have a single amino acid

                pam_acid = string_to_acid[pam_string]
                pam_mutant = [pam_acid, pam_acid]
                mutation_successful, candidate_dna, d_pam = perform_mutation(candidate_dna, pam_loc_in_candidate, 0,
                                                                             pam_mutant,
                                                                             mutation_location=mutation_location)
                if config.TRACE_CANDIDATE_DNA_GENERATION:
                    print("PAM Candidate DNA:")
                    print(temp_candidate_dna)
                if not mutation_successful:
                    if config.VERBOSE_EXECUTION:
                        print('Failed to find a valid replacement for the pam')
                    gs.failed_due_to_pam += 1
                    return None
            else:
                if pam_case == 1:  # the N is in one acid and the GG is in another so we can only replace down
                    # TODO:  make sure we don't overwrite the mutation by removing the pam
                    pam_string_up = None  # dna[pam - 1:pam + 2]
                    pam_string_down = dna[pam + 1:pam + 4]

                elif pam_case == 2:  # the NG is in one acid and the G is in another
                    pam_string_up = dna[pam - 1:pam + 2]
                    pam_string_down = dna[pam + 2:pam + 5]

                replaceable_pam = False
                pam_mutant_up = None
                if pam_string_up is not None and pam_string_up in string_to_acid:
                    pam_acid_up = string_to_acid[pam_string_up]
                    replaceable_pam = True
                    pam_mutant_up = [pam_acid_up, pam_acid_up]
                if pam_string_down in string_to_acid:
                    pam_acid_down = string_to_acid[pam_string_down]
                    replaceable_pam = True
                    pam_mutant_down = [pam_acid_down, pam_acid_down]
                if not replaceable_pam:
                    return None
                if pam_mutant_up is not None:
                    mutation_successful, temp_candidate_dna, d_pam = perform_mutation(candidate_dna,
                                                                                      pam_loc_in_candidate - 1,
                                                                                      pam_case, pam_mutant_up,
                                                                                      mutation_location=mutation_location)
                    if config.TRACE_CANDIDATE_DNA_GENERATION:
                        print("PAM Candidate DNA 2:")
                        print(temp_candidate_dna)
                    if mutation_successful:
                        candidate_dna = temp_candidate_dna
                if not mutation_successful:  # try downstream if upstream didn't work
                    if config.TRACE_CANDIDATE_DNA_GENERATION:
                        print("PAM Candidate DNA 3:")
                        print(temp_candidate_dna)
                    mutation_successful, temp_candidate_dna, d_pam = perform_mutation(candidate_dna,
                                                                                      pam_loc_in_candidate + pam_case,
                                                                                      pam_case, pam_mutant_down,
                                                                                      mutation_location=mutation_location)
                    if mutation_successful:
                        candidate_dna = temp_candidate_dna
        else:
            mutation_successful = True  # we mutated the PAM already, so we didn't need to do another mutation
        if mutation_successful:
            pam_loc = pam_loc_in_candidate
            if complement:  # if we are on the reverse complement, invert it back before we add the other stuff
                candidate_dna = invert_dna(candidate_dna)
                guide = invert_dna(guide)
                mutation_location = len(candidate_dna) - mutation_location - 3
                pam_loc = len(candidate_dna) - pam_loc - 3

            candidate_dna = insert_extra_sequence(candidate_dna, guide)
            # we just added 52 + 20 (guide) basepairs
            # guide pam mutation mutationloc dna
            result = MutationTracker(0, pam_loc + 72, mutant, mutation_location + 72, candidate_dna, complement, pam,
                                     d_pam, pam_string)

            # If using the guide library and only allowing one mutation per guide, automatically reject any guide not present in the library
            if (config.USE_GUIDE_LIBRARY and config.ONE_MUTATION_PER_GUIDE and not (
            is_mutation_permitted(guide, result, guide_lib))):
                if config.VERBOSE_EXECUTION:
                    print("Mutation not permitted from guide library")
                gs.failed_due_to_guide_library += 1
            else:
                gs.succeeded += 1
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


# print(invert_dna('GACCGTGCGACTGGGCGTCTCGGATCTAAGCTTTTGAATATTCCCTGTTTGAAGAGCATACGCTCTTCTTCTAACTTGATAAAATAAATATCCAGTCTGATAAATTGACAAGCTCAATTAAATCCAGAAAGCTGAAAGCTGAGGGAATATTCAAAAGCTTACTGGATACGTTGAGGCAATACGATTCGTCGATACAAAATTTAAACATCGAGACGTGTCCCTGCCTTGCG'))

# Now expect all of these to be arrays.
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
            if (config.USE_GUIDE_LIBRARY):
                print('failed due to guide library: ' + str(gs.failed_due_to_guide_library))
            print('succeeded: ' + str(gs.succeeded))

        column_pos = 0

        sheet1 = wb.add_sheet('Mutation Results')
        sheet1.write(column_pos, 0, 'ID')
        sheet1.write(column_pos, 1, 'Mutation From')
        sheet1.write(column_pos, 2, 'Mutation To')
        sheet1.write(column_pos, 3, 'Mutation Location')
        sheet1.write(column_pos, 4, 'Reverse Complement')
        sheet1.write(column_pos, 5, 'Mutation Distance Start of PAM')
        sheet1.write(column_pos, 6, 'Original PAM')
        sheet1.write(column_pos, 7, 'Seed Mutation Distance From PAM')
        sheet1.write(column_pos, 8, 'Result')

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
                sheet1.write(i + column_pos, 1, mutation.mutation[0] + str(int((mutation.pam_location_in_gene + (mutation.mutation_loc - mutation.pam ) - config.GENE_START_BUFFER) / 3) + 1))
                sheet1.write(i + column_pos, 2, mutation.mutation[1])
                sheet1.write(i + column_pos, 3, mutation.pam_location_in_gene + (mutation.mutation_loc - mutation.pam ) - config.GENE_START_BUFFER)
                sheet1.write(i + column_pos, 4, mutation.complement)
                sheet1.write(i + column_pos, 5, str(abs(mutation.mutation_loc - mutation.pam)))
                sheet1.write(i + column_pos, 6, mutation.original_pam)
                sheet1.write(i + column_pos, 7, mutation.distance_from_pam)

                if (mutation.complement):
                    pass
                    # mutation.dna = invert_dna(mutation.dna)
                # Potentially reading positions (and donor) from non inverted gene, when it should be read from reverse compliment?
                # Not fully sure, still gotta look into this // NOTE \\

                seg_first = (mutation.dna[0:len(first)], extra_font)
                seg_guide = (mutation.dna[len(first):len(first) + config.GUIDE_LENGTH], guide_font)
                seg_second = (
                mutation.dna[len(first) + config.GUIDE_LENGTH:len(first) + config.GUIDE_LENGTH + len(second)],
                extra_font)

                if mutation.complement:
                    seg_mutation = (
                    invert_dna(mutation.dna[mutation.mutation_loc: mutation.mutation_loc + 3]), mutation_font)
                else:
                    seg_mutation = (mutation.dna[mutation.mutation_loc: mutation.mutation_loc + 3], mutation_font)

                mutation.distance_from_pam  # this is how far away from the pam we did the mutation
                mode = mutation.pam_location_in_gene % 3  # this is which "mode" we were in, it gives us the offset from pam start
                pam_mut_seg = None
                if mutation.distance_from_pam is not None:  # None means the mutation took care of the pam for us
                    start_of_pam_mutation = mutation.pam + mode - mutation.distance_from_pam
                    pam_mut_seg = (mutation.dna[start_of_pam_mutation:start_of_pam_mutation + 3], pam_mut_font)
                mod_dna2 = 0  # how much we need to remove from dna 2 because of pam mut
                mod_dna3 = 0  # how much we need to remove from dna 3 because of pam mut
                # this is a little hacky, but figure out the ordering of if the mut or pam goes first
                blank = ('', pam_font)
                if mutation.pam < mutation.mutation_loc < mutation.pam + 3:  # we mutated the pam with the OG mutation.
                    seg_pam = [
                        (mutation.dna[mutation.pam: mutation.pam + (mutation.mutation_loc - mutation.pam)], pam_font),
                        blank]
                else:
                    # identify how much, if any of the pam we mutated
                    if mode == 0 and mutation.distance_from_pam == 0:  # we mutated the entire pam
                        seg_pam = [pam_mut_seg, blank]  # we don't need to show the pam
                    elif mode == 1 and mutation.distance_from_pam == 3:  # we mutated the first char
                        seg_pam = [pam_mut_seg, (mutation.dna[mutation.pam + 1: mutation.pam + 3], pam_font)]
                        mod_dna2 = 2
                    elif mode == 1 and mutation.distance_from_pam == 0:  # we mutated the last two characters
                        seg_pam = [pam_mut_seg, (mutation.dna[mutation.pam: mutation.pam + 1], pam_font)]
                        mod_dna3 = 1
                    elif mode == 2 and mutation.distance_from_pam == 0:  # we mutated the last char
                        seg_pam = [(mutation.dna[mutation.pam: mutation.pam + 2], pam_font), pam_mut_seg]
                        mod_dna3 = 2
                    elif mode == 2 and mutation.distance_from_pam == 3:  # we mutated the first two chars
                        seg_pam = [pam_mut_seg, (mutation.dna[mutation.pam + 2: mutation.pam + 3], pam_font)]
                        mod_dna2 = 1
                    else:  # we mutated the seed
                        seg_pam = [(mutation.dna[mutation.pam: mutation.pam + 3], pam_font), blank]
                seg_third = (mutation.dna[len(mutation.dna) - len(third):], extra_font)

                # if we modified the seed instead of the pam, this is what we need to do
                update_dna_2 = False
                if mutation.distance_from_pam is not None and not ((mode == 1 and mutation.distance_from_pam == 3) or (
                        mode == 2 and mutation.distance_from_pam == 3)) and mutation.distance_from_pam > 0:
                    update_dna_2 = True

                if mutation.mutation_loc < mutation.pam:  # upstream mutation
                    seg_dna1 = (
                    mutation.dna[len(first) + config.GUIDE_LENGTH + len(second):mutation.mutation_loc], dna_font)
                    if update_dna_2:
                        seg_dna2 = [(mutation.dna[
                                     mutation.mutation_loc + 3:mutation.pam - mutation.distance_from_pam + mode],
                                     dna_font), pam_mut_seg, (
                                    mutation.dna[mutation.pam - mutation.distance_from_pam + mode + 3:mutation.pam],
                                    dna_font)]
                    else:
                        seg_dna2 = [(mutation.dna[mutation.mutation_loc + 3:mutation.pam - mod_dna2], dna_font), blank,
                                    blank]
                    seg_dna3 = (mutation.dna[mutation.pam + 3 + mod_dna3:len(mutation.dna) - len(third)], dna_font)
                    sheet1.write_rich_text(i + column_pos, 8, (
                    seg_first, seg_guide, seg_second, seg_dna1, seg_mutation, seg_dna2[0], seg_dna2[1], seg_dna2[2],
                    seg_pam[0], seg_pam[1], seg_dna3, seg_third))
                else:  # downstream mutation
                    if update_dna_2:
                        seg_dna1 = [(mutation.dna[
                                     len(first) + config.GUIDE_LENGTH + len(second):mutation.pam - mod_dna2], dna_font),
                                    blank, blank]
                    else:
                        seg_dna1 = [(mutation.dna[
                                     len(first) + config.GUIDE_LENGTH + len(second):mutation.pam - mod_dna2], dna_font),
                                    blank, blank]
                    seg_dna2 = (mutation.dna[mutation.pam + 3 + mod_dna3:mutation.mutation_loc], dna_font)
                    seg_dna3 = (mutation.dna[mutation.mutation_loc + 3:len(mutation.dna) - len(third)], dna_font)
                    sheet1.write_rich_text(i + column_pos, 8, (
                        seg_first, seg_guide, seg_second, seg_dna1[0], seg_dna1[1], seg_dna1[2], seg_pam[0], seg_pam[1],
                        seg_dna2, seg_mutation, seg_dna3, seg_third))

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

    for loc in _dna_locs:
        # I don't think we need to actually find the guides here since it is just pam - 20
        # guides = create_guides(dna, loc)
        # for g in guides:   # for each guide do each mutation

        for m in config.mutations_to_attempt.items():
            mutated_dna = create_mutations(_dna, loc, m, only_once=_only_once)
            if mutated_dna is not None:
                for md in mutated_dna:
                    mutations_output.append(md)

    # class MutationTracker:
    #    guide: int
    #    pam: int
    #    mutation: []
    #    mutation_loc: int
    #    dna: str

    for loc in _inv_dna_locs:
        for m in config.mutations_to_attempt.items():
            mutated_dna = create_mutations(_inv_dna, loc, m, only_once=_only_once,
                                           complement=True)  # returns a mutation tracker
            if mutated_dna is not None:
                # revert to original  -- Maybe we don't need to do this?
                # inv_guide = mutated_dna.guide
                # inv_pam = len(mutated_dna.dna) - mutated_dna.pam
                # inv_mutation = []
                # inv_mutation_loc = len(mutated_dna.dna) - mutated_dna.mutation_loc
                # inv_dna = invert_dna(mutated_dna.dna)
                # inv_mutated_dna = MutationTracker(inv_guide, inv_pam, mutated_dna.mutation, inv_mutation_loc, inv_dna)
                # all_mutations.append(inv_mutated_dna)

                for md in mutated_dna:
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

        all_mutations = get_all_mutations(dna_locs, inv_dna_locs, dna, inv_dna_full)

        # NOTE // Kill guide is inserted as the very last one
        if (len(all_mutations) > 0):
            all_mutations.append(create_kill_guide(all_mutations[0]))

        combined_mutation_page.append(all_mutations)
        # at this point, we have everything we need to output the results
        if not config.OUTPUT_TO_ONE_FILE:
            write_results([frontmatter[i]], [all_mutations], [dna], False)

    if config.OUTPUT_TO_ONE_FILE:
        write_results(frontmatter, combined_mutation_page, dna_list)


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
    all_mutations = get_all_mutations([dna_locs], [inv_dna_locs], dna, inv_dna_full, only_once=True)

    # at this point, we have everything we need to output the results
    write_results([frontmatter[0]], [all_mutations], [dna], False)


if config.RUN_IN_EXECUTION_TESTING_MODE:
    test_execution()
else:
    execute_program()
