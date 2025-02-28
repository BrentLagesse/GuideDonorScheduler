# ----------------------------------------------------
#
#                   MUTATIONS CONFIG
#
# ----------------------------------------------------

# This dictionary stores the mutations to attempt
mutations_to_attempt = dict()

# To add mutations refer to template:
# mutations_to_attempt['mutation_from'] = 'mutation_to'
# Use '*' in place of 'mutation_from' to mutate every acid
mutations_to_attempt['ile'] = 'val'

BP_LENGTH = 132 # This changes the length of the donor sequence

# This is the threshold value of the rank.
# Guides that have rank of less than this threshold value will be ignored
# Rank that is greater than or equal to this threshold value will be used to perform the mutation
RANK_THRESHOLD = 0

# This variable controls whether to use ranks to do the mutations or not
USE_RANK = False

# The rank file must have two columns titled "Guide" and "Rank" that will have guide and its rank
# in that column. If the title for guide column is something different change the "GUIDE_COLUMN_IN_RANK_FILE".
# If the title of rank column is something different change the "RANK_COLUMN_IN_RANK_FILE".
# The order of the guide and rank column doesn't have to be consistent... the guide and rank columns can be
# anywhere in the file as long as the title of them is in the top row and also have to make sure
# the last entry in the guide column has "END" to indicate end of the file.
GUIDE_COLUMN_IN_RANK_FILE = 'Guide1'
RANK_COLUMN_IN_RANK_FILE = 'Rank1'

# We found out that a single, silent seed mutation was not particularly effective.  Use this to set how you
# want to silently mutate the seed/guide if it cannot silently mutate the PAM.
# Options are --
# guide -- mutates everything in the seed/guide that it can
# seed -- mutates everything in the seed that it can
# SILENT_MUTATION_MINIMUM tells us when to fail if we cannot do that number of mutations in the configured mode
SILENT_MUTATION_MODE = 'guide'
SILENT_MUTATION_MINIMUM = 2
SILENT_MUTATION_MAXIMUM = 10   # this is for testing purposes currently

# ----------------------------------------------------
#
#                   FILE SETUP
#
# ----------------------------------------------------

# WARNING !! Guide Library output files will be overwritten when output,
# be sure to move any modifications to input to avoid losing data.

KILL_MODE = False
RANK_FILE = 'guide_only_ranked_ordered_mock'
GUIDE_LIBRARY_OUTPUT_FILE = 'Guide_Library_Output'
GUIDE_LIBRARY_INPUT_FILE = 'Guide_Library_Input'
USE_GUIDE_LIBRARY = False
ONE_MUTATION_PER_GUIDE = True

GUIDE_LIBRARY_EOF = "END"

DEFAULT_IN_FILES = ['S288C_YCL032W_STE50_flanking.fsa',
                    'S288C_YIL144W_NDC80_flanking.fsa',
                    'NUF2_1000flanking.fsa',
                    'TRP1_1000flanking.fsa'] # Set this up

# For KILL MODE testing
if KILL_MODE:
    DEFAULT_IN_FILES = ['cds_from_genomic.fna']
    MAX_PAMS = 3

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
# Add column to guide library output of just guide, just one,
# prefer GATC-guide
# Add priority column, automatically rank the guides in based on their priority

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
if KILL_MODE:
    GENE_START_BUFFER = 0
    GENE_END_BUFFER = 0
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
nns = ['NNS']
stop = ['TAA', 'TAG', 'TGA']

# Kill mutation
KILL_MUTATION = ['*', stop[0]]

if KILL_MODE:
    mutations_to_attempt = dict()
    mutations_to_attempt['*'] = 'stop'

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
