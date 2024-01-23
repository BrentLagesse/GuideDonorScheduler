# ----------------------------------------------------
#
#                   MUTATIONS CONFIG
#
# ----------------------------------------------------

mutations_to_attempt = dict()

# TESTING

#fake mutations for testing
#mutations_to_attempt['glu'] = 'arg'
mutations_to_attempt['asp'] = 'arg'
#mutations_to_attempt['glu'] = 'lys'
#mutations_to_attempt['gly'] = 'phe'
#mutations_to_attempt['glu'] = 'thr'
#mutations_to_attempt['asn'] = 'arg'



# ----------------------------------------------------
#
#                   FILE SETUP
#
# ----------------------------------------------------

# WARNING !! Guide Library ouput files will be overwritten when output,
# be sure to move any modifications to input to avoid losing data.

GUIDE_LIBRARY_OUTPUT_FILE = 'Guide_Library_Output'
GUIDE_LIBRARY_INPUT_FILE = 'Guide_Library_Input'
USE_GUIDE_LIBRARY = False
ONE_MUTATION_PER_GUIDE = True

GUIDE_LIBRARY_EOF = "END"

DEFAULT_IN_FILES = ['S288C_YCL032W_STE50_flanking.fsa',
                    'S288C_YIL144W_NDC80_flanking.fsa',
                    'NUF2_1000flanking.fsa',
                    'TRP1_1000flanking.fsa'] # Set this up
DEFAULT_OUT_FILE = "Mutation_Results_Combined_Output"



OUTPUT_TO_ONE_FILE = True
# Choose to output to multiple outputs or just one

# If false, each gene will be output to an individual file following the
# nameing format of:
# MULTI_OUTPUT_PREFIX + gene
# the output prefix is defined below

MULTI_OUTPUT_PREFIX = "Mutation_Results_"
# \\ NOTE! Only used if outputting to more than one file

QUIT_ON_NO_DATA = True

# Program exits if no files are given


# ----------------------------------------------------
#
#                   OUTPUTS
#
# ----------------------------------------------------

PRINT_MUTATION_RESULTS = True
PRINT_GUIDE_LIBRARY = True

KILL_MUTATION_ID_SUFFIX = 'KILL_GUIDE'
# Add collumn to guide library output of just guide, just one,
# prefer GATC-guide
# Add priority collumn, automatically rank the guides in based on their priority

# collumn 1 = id
# 2 = guide
# 3 = priority // NOTE \\ 

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
#                   DNA SETUP
#
# ----------------------------------------------------

first_sequence = 'GACCGTGCGACTGGGCGTCTCGGATC'
second_sequence = 'GTTTGAAGAGCATACGCTCTTCTTCT'
third_sequence = 'ACATCGAGACGTGTCCCTGCCTTGCG'

GUIDE_LIBRARY_STRAND_PREFIX = 'GATC'
GUIDE_LIBRARY_INVERSE_PREFIX = 'AAAC'

# STOP CODON
# Find a guide that is both early and good, replace some acid with a stop codon
# Must be within first half of gene, earlier is better

# eliminate last half first
# best within first half

# use guide library for this

# In output, colouring may be off for reverse compliment gds // NOTE \\ 

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

# Kill mutation
KILL_MUTATION = ['*', stop[0]]

# ----------------------------------------------------
#
#                   DEBUG AND TESTING
#
# ----------------------------------------------------

RUN_UNIT_TESTS = False
VERBOSE_TESTING = False # Prints progress of unit testing
VERBOSE_EXECUTION = False # Prints progress of main program
PRINT_MUTATION_SUCCESS_COUNTS = True # Print success vs failed mutations after execution
PRINT_MUTATION_CHECKS = False

RUN_IN_EXECUTION_TESTING_MODE = False # Outputs only a few mutations for easier tracing and debug purposes
USE_DEBUG_MUTATION = False # Output ZZZ as mutation for debug purposes
TRACE_CANDIDATE_DNA_GENERATION = False # Print candidate DNA on each step of perform_mutation

DEBUG_INVERSE = False
# NOTES
# DNA - GUIDE ID - MUTATION ID

# 