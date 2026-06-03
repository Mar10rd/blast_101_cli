# blast_101_tests.py
# Automated unit and integration tests for BLAST 101.
# SW tests use BLOSUM62 analytical ground truth.
# BLAST tests validated against NCBI blastp (SwissProt, April 2026).
#
# Author: Mario Antonio Rodriguez Diaz
# MSc Bioinformatics, University of Edinburgh, 2025-2026
#
# Usage:
#   python3 blast_101_tests.py
#   python3 blast_101_cli.py --test

import sys

import programme_settings


#====================================================================
# Test function
#====================================================================
def check(test_name, real, expected):
    '''
    It returns pass/fail message for test cases
    '''
    if real == expected:
        print(f"PASSED {test_name}")
        return True
    else:
        print(f"FAILURE {test_name}")
        print(f"real: {real} | expected: {expected}")
        return False

#======================================================================
# SW tests
#======================================================================
def sw_test():
    '''
    Performs several tests in SW search
    :return:
    '''

    print("\nRunning SW tests")

    # Keep track of SW tests
    SW_register = []

    import smith_waterman_p as SW
    import blosum as bl
    matrix = bl.BLOSUM(62)

    # Test 1 Short sequence test: MARIQ vs MARIQ
    seq1 = "MARIQ"
    score_idem = SW.perform_smith_waterman(seq1,seq1,False,False)
    val_idem = sum([matrix[AA][AA] for AA in seq1])

    SW_register.append(check("SW-1 Short seq test: MARIQ",score_idem,val_idem))

    # Test 2 Missmatch checker: MARIQ vs MARIA
    seq2 = "MARIE"
    score_missmatch = SW.perform_smith_waterman(seq1,seq2,False,False)
    val_missmatch = sum([matrix[seq1[i]][seq2[i]] for i in range(len(seq1))])

    SW_register.append(check("SW-2 Missmatch check: MARIQ vs MARIE",score_missmatch,val_missmatch))

    # Test 3 Finding the best local alignment: QY vs MARIQYMARIA
    val_local = matrix["Q"]["Q"] + matrix["Y"]["Y"]
    score_local = SW.perform_smith_waterman("QY","MARIQYMARIA",False,False)

    SW_register.append(check("SW-3 Best Local Alignment",score_local,val_local))

    return SW_register

#======================================================================
# Blast tests
#======================================================================

def blast_test():
    '''
    Performs several tests in BLAST search
    '''

    print("\nRunning BLAST tests")

    blast_register = []

    # test 1 Find the correct best hit in DB: query AAWVTEGMYSFCY vs uniprot_bit2.fasta
    # query sequence corresponds to CP7A1_HUMAN Cytochrome P450 in uniprot_bit2.fasta
    # blast must find it as the best hit
    import blast_101_search as blast
    import create_seq_word_dict as cwd
    import process_fasta_file as pff

    programme_settings.read()
    programme_settings.settings["DEFAULT"]["database"] = "uniprot_bit2.fasta"

    # Set query
    query = "AAWVTEGMYSFCY"
    blast.qsequence = query
    blast.query_sequence = cwd.create_word_dict(query)
    blast.init_print_timer()
    # Reset res and scores
    pff.res = []
    pff.bestscore = 0

    # run blast
    res = blast.process_fasta_file()

    # checking for at least 1 hit and for the expected hit (accessing header of first element in res)
    blast_register.append(check("BLAST-1: at least one hit returned", len(res) > 0, True))
    blast_register.append(check("BLAST-2: top hit es CYP7A1", "CYP7A1" in res[0][0], True))

    # Test #2 blast_101_search vs blastp ncbi
    # Running both algorithms using the same query should return the same proteins (although scores might vary)
    programme_settings.settings["DEFAULT"]["database"] = "uniprot_sprot.fasta"

    # Set query
    query2 = "QKMRTVFSQAQLCALKDRFQKQKYLSLQQMQELSSILNLSYKQVKTWFQNQRMKCK"
    blast.qsequence = query2
    blast.query_sequence = cwd.create_word_dict(query2)
    blast.init_print_timer()
    # Reset res and scores
    pff.res = []
    pff.bestscore = 0

    # run blast
    res2 = blast.process_fasta_file()

    # Comparing top blast101 hits vs top Blastp NCBI results (08/04/2026)
    # NCBI blastp against SwissProt returns NANOG_MOUSE as top hit for this query
    blast_register.append(check("BLAST-3: at least one hit found in SwissProt", len(res2) > 0, True))

    if len(res2) > 0:
        blast_register.append(check("BLAST-4: top hit es NANOG_MOUSE", "NANOG_MOUSE" in res2[0][0], True))


    return blast_register

#======================================================================
# Main
#======================================================================
def run_tests():
    '''
    Call to all tests and issue a final report
    '''

    results = []

    results += sw_test()
    results += blast_test()

    passed = sum(results)
    total = len(results)

    print("\n" + "=" * 50)
    print(f"  TOTAL: {passed}/{total} tests passed")
    print("=" * 50)


    return


if __name__ == "__main__":
    run_tests()
    sys.exit(0)
