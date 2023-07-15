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
    # Run the program, generate and save a guide library, compare this to confirm validity
    # Next, reload the library back in, remove some entries, and apply it back to the program
    # Ensure only the entries in the library are used and that non entries are discarded


def test_mutation_identification():
    
    # Values to check correctness
    
    locs_to_find = 0 # Number of pams that *should* be returned
    mutations_to_find = 0;
    
    # Test
    
    print("------------------------------ !! WARNING !! ------------------------------")
    print("Mutation locator is not yet configured, values are needed to make this work")
    
    if (config.VERBOSE_TESTING):
        print("Testing mutation Locator")
    
    frontmatter, dna, guide_lib = main.get_dna();
    inv_dna = main.invert_dna(dna)
    locs = np.array(main.get_locations(dna));
    mutations = main.get_all_mutations(locs[0], locs[1], dna, inv_dna)
    
    # Verification

    if (config.VERBOSE_TESTING):
        print("Checking loc count : " + str(len(locs)) + " of " + str(locs_to_find))
        print("Checking mutation count : " + str(len(mutations)) + " of " + str(mutations_to_find))

    if (len(locs) == locs_to_find and len(mutations) == mutations_to_find):
        if (config.VERBOSE_TESTING):
            print("Mutation Locator passed testing")
    else:
        print("Mutation Locator failed testing")
        
    
# Execute all unit tests
    
def run_unit_tests():
    test_finding_pams()
    test_guide_library()
    test_mutation_identification()

if (config.RUN_UNIT_TESTS):
    run_unit_tests()