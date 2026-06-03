# blast_101_cli.py
# Command-line interface for the BLAST 101 protein sequence search application.
# Handles argument parsing, input validation, and mode dispatching.
#
# Author: Mario Antonio Rodriguez Diaz
# MSc Bioinformatics, University of Edinburgh, 2025-2026
# Part of Bioinformatics Algorithms coursework (Dr Simon Tomlinson)
#
# Usage:
#   python3 blast_101_cli.py -q MNAGWTL -d uniprot_bit2.fasta
#   python3 blast_101_cli.py -q MNAGWTL -d uniprot_bit2.fasta --mode sw
#   python3 blast_101_cli.py --test

import argparse
import configparser
import os
import sys

# =====================================================================
# Parser
# =====================================================================
# Creating parser
def create_parser():
    '''
    Creates the argument parser
    '''
    parser = argparse.ArgumentParser(prog="blast_101_cli.py", description='BLAST 101 Search')

    # Creating arguments
    parser.add_argument('-q','--query',help='Query sequence',required=False, default=None)
    parser.add_argument('-d','--database',help='fasta file database',required=False, default=None)
    parser.add_argument('-m','--mode',help='Search mode or stats',choices=['blast','sw','stats'], default='blast')
    parser.add_argument('-t','--test',help='Run automated tests for blast 101', action='store_true')

    return parser


#==============================================================================
# Input Validation
#=============================================================================

# Sequence validation
def validate_sequence(query):
    '''
    Checks if the query sequence is valid
    Check if query sequence:
    -is not empty
    -is not DNA
    -do not have invalid characters
    '''

    # Extracting valid residues from settings
    valid = set("GAVLITSMCPFYWHKRDENQ")
    # DNA bases, but also N for unknown nucleotides and Uracilo, just in case
    dna_bases = set("ATCGNU")
    query=query.upper()
    # Check if sequence is empty
    if len(query) < 1:
        print(f"Error: query sequence is empty")
        sys.exit(1)
    # Check if seq is DNA
    # Since A,C,G and C represent both aminoacids and nucleotides, we most look for protein only characters
    protein_only = set(query) - dna_bases
    if not protein_only:
        print(f"Error: query sequence might be DNA\n Check again for protein entry only")
        sys.exit(1)
    # Check if seq do not have invalid characters
    pos = 0
    for amino_acid in query:
        pos += 1
        if amino_acid not in valid:
            print(f"Error: invalid residue '{amino_acid}' at position {pos}")
            sys.exit(1)
    print(f"Query sequence: {query} is valid")
    return

#-----------------------------------------------------------------------------------
# Database validation
def db_validation(database):
    '''
    Checks if the database is valid

    Check if database file:
    -exists
    -is not empty
    -is fasta format

    '''

    # Check if file exists
    if not os.path.exists(database):
        print("Error in arguments \n database does not exist")
        sys.exit(1)
    # Check if file is not empty
    if os.path.getsize(database) == 0:
        print("Error in arguments \n database is empty")
        sys.exit(1)
    # Check if file is in fasta format: first non-empty line must start with '>'
    # modified from BA7b_ICAandObjects_2per.pdf
    with open(database, "r") as f:
        for line in f:
            line = line.strip()
            if line:  # skip empty lines
                if line[0] != ">":
                    print(f"Error in arguments \n database file is not fasta file")
                    sys.exit(1)
                else:
                    print(f"Database file is valid")
                break  # only check first non-empty line

# ==============================================================================
# Update query and db in settings.ini
#===============================================================================

def update_settings(query, database):
    '''
    Overrides settings file with cli parameters
    '''

    config = configparser.ConfigParser()
    settings_path = "settings.ini"

    if os.path.isfile(settings_path):
    # If file exist, we update query and database by overwriting
        config.read("settings.ini")

        config["DEFAULT"]["query_sequence"] = query.upper()
        config["DEFAULT"]["database"]       = database

        with open(settings_path, "w") as f:
            config.write(f)
        print(f"Settings file '{settings_path}' updated.")
    # If not, a warning message is displayed on screen
    else:
        print(f"Error: settings file '{settings_path}' not found.")
        sys.exit(1)


#=====================================================================
# Printing all parameters before searching
#=====================================================================
def print_args(args):
    '''
    Prints CLI arguments and key BLAST parameters from settings.ini
    '''
    config = configparser.ConfigParser()
    config.read("settings.ini")

    print("-" * 50)
    print("Values for input variables")
    print("-" * 50)
    print(f"  Mode         : {args.mode}")
    print(f"  Query        : {args.query.upper()} ({len(args.query)} aa)")
    print(f"  Database     : {args.database}")
    print(f"  BLOSUM       : {config['BLAST']['blosum']}")
    print(f"  Word size    : {config['BLAST']['word_size']}")
    print(f"  Max scores   : {config['BLAST']['max_scores']}")
    print(f"  Gap score    : {config['DEFAULT']['seq_gap']}")
    print("-" * 50)
    return

#=====================================================================
# Launching modules according to selected mode
#=====================================================================
def run_mode(args):
    '''
    Launches the module corresponding to the search mode
    '''

    if args.mode == "blast":

        try:
            import blast_101_search as blast
            blast.blast101_run()
        except Exception as e:
            print(f"Error in launching blast: {e}")
            sys.exit(1)
    elif args.mode == "sw":
        try:
            import smith_waterman_search as smith
            smith.sw_run()
        except Exception as e:
            print(f"Error in launching sw: {e}")
            sys.exit(1)
    elif args.mode == "stats":
        try:
            import calc_bit_and_evalues as cbe
            cbe.build_fit()
        except Exception as e:
            print(f"Error in launching stats: {e}")
            sys.exit(1)

#================================================================================
# Main
#===============================================================================
def main():

    args = create_parser().parse_args()

    # Welcoming user
    print("*"*50)
    print("Welcome to BLAST 101 Command Line Interface")
    print("*"*50)

    if args.test:
        import blast_101_tests
        blast_101_tests.run_tests()
        sys.exit(0)
    if args.mode == "stats":
        run_mode(args)
        return
    if args.query is None or args.database is None:
        print("Error: -q/--query and -d/--database are required for blast and sw modes")
        sys.exit(1)
    else:
        validate_sequence(args.query)
        db_validation(args.database)
        update_settings(args.query, args.database)
        print_args(args)
        run_mode(args)

    return
    
    
if __name__ == "__main__":
    main()
