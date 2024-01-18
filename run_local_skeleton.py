#!/usr/bin/python

"""
Vrije Universiteit Amsterdam, 2021
Principles of Bioinformatics, X_401094, Assignment 2 Skeleton Script.

Run (PSI-)BLAST locally by using the blastp and psiblast
terminal commands through this script.

example script call (running BLAST), add `-psi` for PSI-BLAST:
python3 run_local_skeleton.py \
    -i data/SCOP_selections.txt \
    -db data/myDatabase.fasta \
    -o results/blastEvalues.tsv

example script call (plotting):
python3 run_local_skeleton.py \
    -r ./results/blastEvalues.tsv \
    -o ./results/BLAST_histogram.png
"""

from argparse import ArgumentParser
from datetime import datetime
import json
from subprocess import Popen, PIPE, STDOUT
from typing import Iterable, Tuple
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

FILE_EXTENSIONS: Tuple[str, ...] = ('.tsv', '.txt')
PLOT_EXTENSIONS: Tuple[str, ...] = ('.png', '.jpeg')


# Helper Functions
def argument_parser() -> ArgumentParser:
    """Provide an interface to parse arguments
    into the script's Namespace.

    :return:    ArgumentParser object with possible
                script arguments.
    """
    parser = ArgumentParser(description="""
    Run BLAST and PSI-BLAST through a Python script. 
    The generated output is captured and stored into 
    a properly formatted file.
    """)

    # Input/Output arguments
    parser.add_argument(
        "-i", "--inputFilename", type=str,
        help="File location for the input data.",
        default="")
    parser.add_argument(
        "-r", "--algoResult", type=str,
        help="Input data for the plotting function if a name is provided.",
        default="")
    parser.add_argument(
        "-q", "--queriesFolder", type=str,
        help="Output folder for your individual query storage. "
             "Default = ./queries/",
        default="./queries/")
    parser.add_argument(
        "-db", "--databaseFilename", type=str,
        help="File name of the database.",
        default="")
    parser.add_argument(
        "-o", "--outFilename", type=str,
        help=f"Output file name. Use either {' or '.join(FILE_EXTENSIONS)} "
             f"to store the results from running an algorithm and either "
             f"{' or '.join(PLOT_EXTENSIONS)} to store plots.",
        default="")

    # Optional arguments
    parser.add_argument(
        "-psi", "--runPsiblast",
        help="Run PSI-BLAST instead of BLAST.",
        action="store_true")
    parser.add_argument(
        "-t", "--threads", type=str,
        help="The amount of threads that (PSI-)BLAST can use.",
        default="1")
    parser.add_argument(
        "-s", "--silent",
        help="Mute the progress report.",
        action="store_true")

    return parser


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


def unique_cartesian_products(*args: Iterable,
                              repeats: int = 1) -> Iterable:
    """Create a cartesian product as we would
    by using itertools but ensure we only get unique pairs.

    :param args:    iterables to combine.
    :param repeats: amount of time to include a cartesian product.

    :return:    A generator that holds as many unique cartesian
                products as specified by repeats.
    """
    pools = [tuple(pool) for pool in args] * repeats
    result = [[]]
    for pool in pools:
        result = [x + [y] for x in result for y in pool if x != [y]]
    for prod in result:
        yield tuple(prod)


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


def parse_blast(blast_results: str) -> Tuple[list, dict]:
    """Parse the BLAST E-values file into a usable dict.

    :param blast_results:  location of the blast output file

    :return: mapping of pair to E-value
    """
    with open(blast_results, "r") as blast_content:
        metadata = blast_content.readline().split("\t")
        return metadata, {
            tuple(line.split("\t")[:2]): float(line.split("\t")[-2])
            for line in blast_content
        }


def write_output_file(proteins: Iterable,
                      metadata: str,
                      algo_result: dict,
                      output: str) -> None:
    """This function writes the scores of
    all-against-all protein pairs to the specified
    output file.

    :param proteins:    list of all uniprot IDs
    :param metadata:    what created these values
    :param algo_result: mapping of (query, subject) to e-value
    :param output:      file name to store output

    :return: None
    """
    # See the `unique_cartesian_products` docstring for more info.
    combinations = unique_cartesian_products(
        proteins, repeats=2)
    
    # Create dir if it doesnt exist
    dir = Path(output).parent
    if not dir.exists(): dir.mkdir()


    with open(output, "w") as store:
        # For every available pair, write the
        # corresponding e-value: 'protein1\tprotein2\t0.987'.
        store.write(metadata + "\n")
        for pair in combinations:
            if pair in algo_result:
                store.write(
                    '\t'.join(pair + (str(algo_result[pair]), "\n"))
                )

    return


# Functions that require your input
def blast(query: str,
          queries: str,
          database: str,
          psi: bool,
          threads: str) -> dict:
    """This function executes blast or
    psi-blast for the given query and db.

    :param query:		    query protein name
    :param queries:         queries folder
    :param database: 		database file name
    :param psi:		        perform psiblast when True, else do blastp
    :param threads:         amount of threads for the BLAST executable

    :return:    the (psi-)blast results in a dictionary
                with a `(query, subject): E-value` mapping.
    """
    cmd = ""  # Dummy value

    if psi:

        #####################
        # START CODING HERE #
        #####################

        # Create a command with the appropriate arguments
        # for running PSI-BLAST. Store this command inside
        # `cmd = "..."`
        # Take a good look at `psiblast -help` to
        # check the available options.
        cmd = f"psiblast -query {queries}{query} -db {database} -outfmt 15 -out {queries}{query}.json"

        # Ensure the *output* format is a 'Single-file BLAST JSON'!
        # Hint: listed under `-outfmt`!

        #####################
        #  END CODING HERE  #
        #####################

    else:

        #####################
        # START CODING HERE #
        #####################

        # Create a command with the appropriate arguments
        # for running BLAST. Store this command inside
        # `cmd = "..."`
        # Take a good look at `blastp -help` to
        # check the available options.
        cmd = f"blastp -query {queries}{query} -db {database} -outfmt 15 -out {queries}{query}.json"

        # Ensure the *output* format is a 'Single-file BLAST JSON'!
        # Hint: listed under `-outfmt`!

        #####################
        #  END CODING HERE  #
        #####################

    # Run a shell command from the script.
    p = Popen(
        cmd + " -num_threads " + threads,
        shell=True,
        stdin=PIPE,
        stdout=PIPE,
        stderr=STDOUT,
        close_fds=True
    )
    p.wait()  # wait for the process to finish

    # Load the JSON content as a string
    # and parse it as a dictionary.
    try:
        with open(f"{queries}{query}.json", "r") as f:
            content = json.load(f)
    except Exception as exc:
        raise ValueError(
            "Looks like your CMD did not run successfully. "
            f"Encountered: {exc}"
        )

    #####################
    # START CODING HERE #
    #####################

    # Create a dict such that you will produce:
    # result = {
    #   (query, subject1): float(E-value),
    #   (query, subject2): float(E-value),
    #   ...
    # }
    result = dict()
    
    # Where subject1 refers to the 'accession' of a hit.

    # You can use the `content` variable (dict) to retrieve the items.
    # For psi-blast, you can check:
    # `content['BlastOutput2'][0]['report']['results']['iterations'][-1]['search']['hits']`
    if psi:
        for i in range(len(content['BlastOutput2'][0]['report']['results']['iterations'])):
            for j in range(len(content['BlastOutput2'][0]['report']['results']['iterations'][i]['search']['hits'])):
                result[(query, content['BlastOutput2'][0]['report']['results']['iterations'][i]['search']['hits'][j]['description'][0]['accession'])] = float(content['BlastOutput2'][0]['report']['results']['iterations'][i]['search']['hits'][j]['hsps'][0]['evalue'])
    # For blast, you can check:
    # `content['BlastOutput2'][0]['report']['results']['search']['hits']`
    else:
        for i in range(len(content['BlastOutput2'][0]['report']['results']['search']['hits'])):
            result[(query, content['BlastOutput2'][0]['report']['results']['search']['hits'][i]['description'][0]['accession'])] = float(content['BlastOutput2'][0]['report']['results']['search']['hits'][i]['hsps'][0]['evalue'])

    #####################
    #  END CODING HERE  #
    #####################

    return result


def plot_value_distribution(blast_result: tuple,
                            output_name: str,
                            ) -> Tuple[plt.Figure, plt.Axes, int]:
    """This function plots the distribution of the log(e-value).

    :param blast_result:    mapping of (query, subject) to e-value
    :param output_name:     location to save the e-value distribution figure

    :return:    a matplotlib figure object
    """
    metadata, mapping = blast_result
    array = np.array(list(mapping.values()))

    # Create a histogram with all log10(e-values)
    fig, ax = plt.subplots()
    frequency = 0  # Dummy Value
    ax.hist(
        np.log10(
            array + (np.min(array[np.nonzero(array)]) / 1000)
        )
    )

    #####################
    # START CODING HERE #
    #####################

    # Give the plot appropriate labels. Check `metadata`
    # for more info on the algorithm that was used :)
    ax.set_xlabel("log10(E-value)")
    ax.set_ylabel("Frequency")
    ax.set_title(f"Distribution of log10(E-values) for {metadata[0]}")

    # ax.set_xlabel()
    # ax.set_ylabel()
    # ax.set_title()

    # (Required for Question 3.1)
    # frequency = ...
    frequency = np.sum(np.log10(array + (np.min(array[np.nonzero(array)]) / 1000)) < -3)

    # If frequency is greater than 0, it will
    # automatically be printed to the terminal.
    # Hint: You can use `np.sum()` to count
    # the amount of values in an array that match
    # a statement: e.g. np.sum(np.array([0,0,1]) > 0) = 1

    #####################
    #  END CODING HERE  #
    #####################

    if output_name: fig.savefig(output_name)

    return fig, ax, frequency


# No additional coding input required!
def main() -> None:
    """Run this function through the
    `if name == main` statement."""
    args = vars(argument_parser().parse_args())

    if args['inputFilename'] and args['databaseFilename']:
        # Run the specified algorithm.
        if not args['outFilename'].lower().endswith(FILE_EXTENSIONS):
            raise ValueError(
                f"Make sure the results file extension is either {' or '.join(FILE_EXTENSIONS)}!"
            )

        proteins = parse_scop(args['inputFilename'])
        algo_result = dict()

        for query in progress(proteins, args['silent']):
            algo_result.update(
                blast(
                    query,
                    args['queriesFolder'],
                    args['databaseFilename'],
                    args['runPsiblast'],
                    args['threads']
                )
            )

        # Write the results to a file on disk.
        # Store the (PSI-)BLAST results with some metadata.
        write_output_file(
            proteins,
            "\t".join((
                "PSI-BLAST" if args['runPsiblast'] else "BLAST",
                datetime.now().strftime("%m/%d/%Y\t%H:%M:%S")
            )),
            algo_result,
            args['outFilename']
        )

        # Set the plotting values properly.
        # - args.r will be the provided args.o
        # - args.p will be an arbitrary plot name
        args['algoResult'], args['outFilename'] = (
            args['outFilename'],
            f"value_distribution_{'blast' if not args['runPsiblast'] else 'psiblast'}.png"
        )

    if args['algoResult']:  # Plot E-value distribution!
        if not args['outFilename'].lower().endswith(PLOT_EXTENSIONS):
            raise ValueError(
                f"Make sure the plot file extension is either {' or '.join(PLOT_EXTENSIONS)}!"
            )

        _, _, freq = plot_value_distribution(parse_blast(args['algoResult']), args['outFilename'])

        if not args['silent']: print(f"Successfully stored a new plot in '{args['outFilename']}'.")
        if freq and not args['silent']: print(f"E-values *lower* than threshold: {freq}")

    else:
        # In case no arguments are (correctly) provided.
        print("ERROR: Did not receive any valid parameters.")

    return


if __name__ == "__main__": main()

