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


# ----------------------------------------------------
#
#                   PARAMETER CONFIG
#
# ----------------------------------------------------


GENE_START_BUFFER = 1000
GENE_END_BUFFER = 1000
UP_ACIDS = 6
DOWN_ACIDS = 4

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