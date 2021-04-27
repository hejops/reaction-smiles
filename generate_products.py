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


def main(db, startidx=0, strict=False):
    """
    Extract all reactants from a csv file and generate products

    By default, file.csv will be written to file_out.csv
    """

    if not Path(db).is_file():
        print(f"{db} is not a valid file!")
        sys.exit()

    # # TODO: move this to csvops.py
    # with open(db) as f:
    #     # TODO: detect if there is a header
    #     reactions = [row.split(",")[0] for row in f][1:]
    # reactions = reactions[startidx:]

    csv_data = []
    warn = 0
    start_time = time.time()

    reader = csv.DictReader(open(db))
    t = sum(1 for row in open(db))
    # print(row_count)
    # t = len(list(reader))
    # with open(db) as f:
    #     reader = csv.DictReader(f)

    for i, row in enumerate(reader):
        N = row["N"]
        reaction = row["SMILES"]

        # TODO: return name as list
        i += startidx
        print(f"\rReaction {i}: ", end="")
        full_smiles, name, sanit = get_full_reaction(reaction)

        # might not need to return sanit if lowercase check is used
        if strict and "c" in full_smiles:
            sys.exit()

        if sanit:
            print(f"{i} OK")

        else:
            print(f"{i} OK, but not sanitized")
            warn += 1

        print(dashes)

        csv_data.append(
            {
                "Line": i + 1,
                "Reaction (N)": N,
                "Reactant 1": name[0],
                "Reactant 2": name[1],
                "Sanitized": sanit,
                "Reactants": reaction,
                "Full": full_smiles,
            }
        )

        if i == 100:
            break

        # TODO: if smarts <8 colons, warn simple; this is not so easy because smarts is no longer returned

    # TODO: because there are duplicate reactions, results always < rows
    # fails = [x for x in reactions if x not in passed_rxns]

    print(
        f"""
    {len(csv_data)}/{t} rows processed successfully
    Warning: {warn} products were unsanitized
    Finished in {time.time() - start_time} seconds
"""
    )
    # Note: {warn} reactions were obtained with relatively simple rules

    # print(csv_data)

    outfile = db.replace(".csv", "_out.csv")
    headers = list(csv_data[0])
    with open(outfile, "w") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=headers)
        writer.writeheader()
        for data in csv_data:
            writer.writerow(data)

    print(f"Wrote to {outfile}\n")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=main.__doc__)

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
        "-C",
        "--columns",
        const="col",
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
    # print(args)
    main(args.file, args.line, args.strict)
