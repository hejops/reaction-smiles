#!/usr/bin/env python3
# from rdkit import Chem
# from rdkit.Chem.Draw import rdMolDraw2D
# import matplotlib.image as mpimg	    much slower than PIL
# import matplotlib.pyplot as plt
# import os
# import psutil
# import re
# import readline
# import requests
# import urllib
from pathlib import Path
import argparse
import csv
import json
import logging
import pickle
import shutil
import sys
import time

# from drawops import *	        # draw images
# from molops import *	        # apply actions to molecules (e.g. unmap, neutralise); deprecated
from smartjoin import get_full_reaction

logging.basicConfig(stream=sys.stderr, level=logging.ERROR)

dashes = "-" * shutil.get_terminal_size().columns

# with open(pkl, 'rb') as handle:
#     b = pickle.load(handle)


def dump_all(x):
    with open("test.pickle", "wb") as handle:
        pickle.dump(x, handle, protocol=pickle.HIGHEST_PROTOCOL)
    with open("test.json", "w") as handle:
        json.dump(x, handle)


# Previously, reactions were split to obtain unique reactants for manual
# mapping/joining. This has been deprecated in favour of a more robust and
# maintainable approach using SMARTS.

# # flatten, remove header
# reactants = [i.split(".") for i in reactions]
# reactants = [x for sublist in reactants for x in sublist]
# seen = set()
# seen_add = seen.add
# reactants = [x for x in reactants if not (x in seen or seen_add(x))]


def main(db, col, startidx, strict):
    """
    Extract all reactants from a csv file and generate products

    By default, [file].csv will be written to [file]_out.csv
    """

    if not Path(db).is_file():
        print(f"{db} is not a valid file!")
        sys.exit()

    warn = 0
    passed_rxns = []
    failed_rxns = []
    start_time = time.time()

    reader = csv.DictReader(open(db))
    t = sum(1 for row in open(db))

    for i, row in enumerate(reader):
        # col = "SMILES"
        N = row["N"]
        reaction = row[col]

        if ">>" in reaction:
            print(f"The {col} column of {db} appears to contain products")
            sys.exit()

        if i < startidx:
            continue

        print(f"\rReaction {i}: ", end="")
        full_smiles, name, sanit = get_full_reaction(reaction)

        # might not need to return sanit if lowercase check is used
        if strict and "c" in full_smiles:
            sys.exit()

        if sanit:
            print(f"{i} OK")

        else:
            print(f"{i} OK, but not sanitized")
            failed_rxns.append(
                {
                    "Line": i + 1,
                    "Reaction": full_smiles,
                }
            )
            warn += 1

        print(dashes)

        passed_rxns.append(
            {
                "Line": i + 1,
                "N": N,
                "Reactant 1": name[0],
                "Reactant 2": name[1],
                "Sanitized": sanit,
                "Reactants": reaction,
                # keep headers consistent between in and out file
                col: full_smiles,
            }
        )

        # if i == 100:
        #     break

        # TODO: if smarts <8 colons, warn simple; this is not so easy because smarts is no longer returned

    # TODO: because there are duplicate reactions, results always < rows
    # fails = [x for x in reactions if x not in passed_rxns]

    print(
        f"""
    {len(passed_rxns)}/{t} rows processed successfully
    Warning: {warn} products were unsanitized
    Finished in {time.time() - start_time} seconds"""
    )
    # Note: {warn} reactions were obtained with relatively simple rules

    def dict_to_csv(listofdicts, file):
        if not listofdicts:
            return

        with open(file, "w") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=list(listofdicts[0]))
            writer.writeheader()
            for data in listofdicts:
                writer.writerow(data)
        print(f"Wrote to {file}\n")

    dict_to_csv(passed_rxns, db.replace(".csv", "_out.csv"))
    dict_to_csv(failed_rxns, db.replace(".csv", "_err.csv"))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description=main.__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        # formatter_class=SortingHelpFormatter,
    )

    parser.add_argument(
        "file",
        help="csv file",
    )

    # parser.add_argument(
    #     # "-o",
    #     "--outfile",
    #     const="startidx",
    #     action="store_const",
    #     help="output filename",
    # )

    parser.add_argument(
        "-c",
        "--column",
        # const="col",
        # action="store_const",
        default="SMILES",
        # nargs="+",
        help="name of column containing SMILES to be processed",
    )
    parser.add_argument(
        "-C",
        "--extra-columns",
        const="cols",
        action="store_const",
        # nargs="+",
        help="specify additional column header(s) to be extracted from csv",
    )
    parser.add_argument(
        "-l",
        "--line",
        metavar="N",
        default=0,
        type=int,
        help="start from line (not counting header)",
    )
    parser.add_argument(
        "-s",
        "--strict",
        action="store_true",
        default=False,
        help="exit when a bad sanitization is encountered",
    )

    args = parser.parse_args()
    main(args.file, args.column, args.line, args.strict)
