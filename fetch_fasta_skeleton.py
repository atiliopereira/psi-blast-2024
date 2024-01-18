#!/usr/bin/python

"""
Vrije Universiteit Amsterdam, 2021
Principles of Bioinformatics, X_401094, Assignment 2 Skeleton Script.

This script downloads all sequences that were provided by the user,
putting them into a single file to be used as a database and
storing the fetched fasta sequences into individual files
which will be used as queries for the (PSI-)BLAST search.

example script call:
python3 fetch_fasta_skeleton.py \
    -i data/SCOP_selections.txt \
    -db data/myDatabase.fasta
"""

from argparse import ArgumentParser
from pathlib import Path
from typing import Iterable, Tuple

import httpx
from tqdm import tqdm


FILE_EXTENSIONS: Tuple[str, ...] = ('.tsv', '.txt')
PLOT_EXTENSIONS: Tuple[str, ...] = ('.png', '.jpeg')
UNIPROT_URL: str = "https://rest.uniprot.org/uniprotkb/"


# Helper Functions
def argument_parser() -> ArgumentParser:
    """Provide an interface to parse arguments
    into the script's Namespace.

    :return:    ArgumentParser object with possible
                script arguments.
    """
    parser = ArgumentParser(description="""
    Download FASTA-formatted files for a list 
    of UniProt Accessions.
    """)

    # Input/Output Arguments
    parser.add_argument(
        "-i", "--inputFilename", type=str,
        help="File location for the input data.",
        default="")
    parser.add_argument(
        "-q", "--queriesFolder", type=str,
        help="Output folder for your individual query storage. "
        "Default = ./queries/",
        default="./queries/")
    parser.add_argument(
        "-db", "--databaseFilename", type=str,
        help="File name for the database.",
        default="")

    # Optional Arguments
    parser.add_argument(
        "-s", "--silent",
        help="Mute the progress report.",
        action="store_true")

    return parser


def test_link() -> bool:
    """Test whether the uniprot URL is available.

    :return: True if the url is working.
    """
    url: str = f"{UNIPROT_URL}P12345.fasta"

    return httpx.get(url).status_code == 200


def progress(sequence: Iterable,
             mute: bool) -> Iterable:
    """Check if the user wants to mute the
    download progression.

    :param sequence:    the iterable to follow
    :param mute:        does the user want to mute

    :return: a tracktable iterable if desired
    """
    if not mute:
        sequence = tqdm(sequence)

    return sequence


def parse_scop(scop_file: str) -> Iterable:
    """Parse the file content from the
    SCOP selection into an iterable.

    :param scop_file:   location of the SCOP file.

    :return: all UniProt IDs in a single Iterable.
    """
    lines = list()
    with open(scop_file, "r") as file_content:
        for line in file_content:
            lines.append(line.strip())

    return lines


# Functions that require your input
def fetch_protein_sequences(proteins: Iterable,
                            query_folder: str,
                            database_filename: str,
                            mute: bool) -> None:
    """Fetch the fasta formatted sequence for each
    UniProtID and store the content as a `.fasta`
    inside the queries folder and create a single
    database by appending all individual fastas
    into 1 file.

    :param proteins:            a sequence of uniprot ids
    :param query_folder:        folder containing fasta files
    :param database_filename:   path to the file saving the database
    :param mute:                disable progress tracking

    :return: None
    """
    # Make sure we point to a directory and not a file.
    if not query_folder.endswith("/"): query_folder += "/"
    # Create the folder if it does not exist
    if not Path(query_folder).exists(): Path(query_folder).mkdir()

    with open(database_filename, "a") as database, \
            httpx.Client() as client:
        # For every uniprot ID:
        for protein_id in progress(proteins, mute):

            #####################
            # START CODING HERE #
            #####################

            raise ValueError("I need to be removed!")
            # 1. Fetch the fasta-formatted file for each UniProt ID (protein_id).
            # You can use `fasta = client.get(url).text` to retrieve the
            # file from an URL. Make sure to set the `url` variable correctly!
            # According to: https://www.uniprot.org/help/api_retrieve_entries
            # Hint: take a look at the 'test_link()' function.

            # 2. Store the fasta content as individual fasta files in your **queries**
            # directory. (Available through the `query_folder` variable (string))

            # 3. Append the fasta content to the 'database' (a single file).
            # You have access to the database file through the `database`
            # variable (file object). Note: the database variable has the file
            # opened in append mode!

            #####################
            #  END CODING HERE  #
            #####################

    return


def main() -> None:
    """Run this function through the
    `if name == main` statement."""
    if not test_link():
        raise ValueError(
            "The URL seems to be dysfunctional. Please check with the teacher.")

    # Set script arguments inside a dictionary.
    args = vars(argument_parser().parse_args())

    if ".fasta" not in args['databaseFilename']:
        raise ValueError(
            "Please provide a Database name with a .fasta extension!"
        )

    # Get a collection of proteins.
    proteins = parse_scop(args['inputFilename'])

    if not args['silent']: print("Gathering protein FASTA files from Uniprot!")

    fetch_protein_sequences(
        proteins,
        args['queriesFolder'],
        args['databaseFilename'],
        args['silent']
    )

    if not args['silent']: print(
        f"Done... Downloaded {len(list(Path(args['queriesFolder']).glob('*.fasta')))} samples."
    )

    return


if __name__ == "__main__": main()

