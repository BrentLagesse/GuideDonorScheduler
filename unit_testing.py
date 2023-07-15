import config
import main
import numpy as np

# ----------------------------------------------------
#
#                   UNIT TESTING
#
# ----------------------------------------------------

def test_finding_pams():
    
    if (config.VERBOSE_TESTING):
        print("Testing PAM Locator")
    
    frontmatter, dna, guide_lib = main.get_dna();
    inv_dna = main.invert_dna(dna)
    locs = np.array(main.get_locations(dna));
    pams = []
    success = True
    for l in locs[0]:
        p = dna[l:l+3]
        pams.append(p)
        if (not p[1:] == "GG"):
            success = False
    for l in locs[1]:
        p = inv_dna[l:l+3]
        pams.append(p)
        if (not p[:-1] == "CC"):
            success = False
    if (config.VERBOSE_TESTING):
        print("PAMs identified:")
        print(pams)
    if (success):
        if (config.VERBOSE_TESTING):
            print("Pam Locator passed testing")
    else:
        print("Pam Locator failed testing")

def test_guide_library():
    pass

# Execute all unit tests
    
def run_unit_tests():
    test_finding_pams()
    

if (config.RUN_UNIT_TESTS):
    run_unit_tests()