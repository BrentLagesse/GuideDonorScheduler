# ----------------------------------------------------
#
#                   MUTATIONS CONFIG
#
# ----------------------------------------------------

mutations_to_attempt = dict()




#fake mutations for testing
mutations_to_attempt['met'] = 'val'
mutations_to_attempt['ala'] = 'tyr'
#mutations_to_attempt['*'] = 'arg'
mutations_to_attempt['cys'] = 'thr'
mutations_to_attempt['val'] = 'ile'

first_sequence = 'GACCGTGCGACTGGGCGTCTCGGATC'
second_sequence = 'GTTTGAAGAGCATACGCTCTTCTTCT'
third_sequence = 'ACATCGAGACGTGTCCCTGCCTTGCG'

# ----------------------------------------------------
#
#                   FILE SETUP
#
# ----------------------------------------------------

GUIDE_LIBRARY_OUTPUT_FILE = 'guide_library'
GUIDE_LIBRARY_INPUT_FILE = 'guide_library'
USE_GUIDE_LIBRARY = False

DEFAULT_IN_FILE = 'S288C_YCL032W_STE50_flanking.fsa'
DEFAULT_OUT_FILE = 'mutation_output'

QUIT_ON_NO_DATA = True

MUTATION_RESULT_ID_PREFIX = 'GD'
GUIDE_LIBRARY_ID_PREFIX = 'GL'

# ----------------------------------------------------
#
#                   OUTPUTS
#
# ----------------------------------------------------

PRINT_MUTATION_RESULTS = True
PRINT_GUIDE_LIBRARY = True



# ----------------------------------------------------
#
#                   PARAMETER CONFIG
#
# ----------------------------------------------------


GENE_START_BUFFER = 1000
GENE_END_BUFFER = 1000
UP_ACIDS = 6
DOWN_ACIDS = 4

GUIDE_LENGTH = 20

# ----------------------------------------------------
#
#                   CODON SETUP
#
# ----------------------------------------------------


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
stop = ['TAA', 'TAG', 'TGA']





# FINISHED ISSUES ( Ignore, this is for me )

# Skip mutating first acid complete
# Guide library ( Though I should add an export that could then be read back in )
# Config ( Will likely be adding more to this as time goes on however )
# Improved output ( Still need to confirm IDs, finish that * thing?? )
# Reverse strand needs reverse compliment ( The fix seemed too easy, probably not right, gotta look further into it + the additional tasks in that issue )