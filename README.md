# GuideDonorScheduler
Creates Guide+Donor for gene editing.

### **Overview**

This code is designed to facilitate the process of mutating a specific gene and generating corresponding donor sequences for CRISPR/Cas9 gene editing. The mutations include changes to both the target gene and the adjacent Protospacer Adjacent Motif (PAM). The code follows specific rules for selecting guide sequences and designing donor sequences around the mutations.

### **Usage**
1. #### Choose Target Gene and Guide Sequence: 
   * Choose the target gene (e.g., Ndc80) and identify the guide sequence by locating the NGG PAM sequence nearby.
     * ```python main.py -i INPUT_FILE -o OUTPUT_FILE```
       * The input file is in .fsa format
       * The output file is in xls format
       * If either input or output files are omitted, DEFAULT_IN_FILES and DEFAULT_OUT_FILE in config.py will be used
   * Identify the amino acid to be mutated and specify the desired mutation (e.g., Serine 49 to Alanine).
     * Edit the mutations_to_attempt dict in config.py.  (e.g., ```mutations_to_attempt['ser'] = ['ala']``` to mutate Serine to Alanine) 
2. #### Guide and Donor Design:
   * The guide sequence is the 20 bases upstream of the NGG PAM.
   * Design the donor sequence by selecting 132 base pairs surrounding the mutations, centered around both the main mutation and the PAM mutation.
3. #### Reverse Complement Strand:
   * For mutations on the reverse complement strand, reverse the guide sequence and design the donor accordingly.
4. #### Final Guide+Donor:
   * Combine the guide and donor sequences, considering any additional consistent sequence.
5. #### Guide PAM Mutation for Cas9 Avoidance:
   * The code is designed to introduce silent mutations in the PAM (NGG) sequence, preventing Cas9 re-cutting.
   * The preference is to mutate either G in the PAM to any other base.
   * Caution is advised against mutating PAM to NAG.
   * If modifying the PAM is challenging, the code targets at least one location in the "seed" region (10 bases upstream of the PAM).
   * Consider performing two silent mutations in the seed region if a single mutation proves insufficient to prevent re-cutting.
   * If the intended mutation is in the seed, the code prioritizes silent mutations.
   * If the mutation is more than 5 bases away from the PAM, the code aims for a silent mutation within 5 bases upstream of the PAM.
     * This ensures that, even though the mutation itself is not in the immediate vicinity of the PAM, a silent mutation is made in a nearby region to potentially influence the Cas9 activity.
   * If the mutation is within 5 bases of the PAM and a silent PAM mutation is not achievable, the code may refrain from additional mutation.
     * In cases where the mutation is already close to the PAM (within 5 bases) and it is not feasible to perform a silent mutation in the PAM itself, the code may decide not to introduce any additional mutation.
