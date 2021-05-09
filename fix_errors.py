#!/usr/bin/env python3
# import csv
# import requests
# import urllib
from pyfzf.pyfzf import FzfPrompt
import logging
import os
import pandas as pd
import readline
import shutil
import sys

# logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)

from smartjoin import *
from drawops import *


def get_smarts(reactants: str, line: str):
    delim = " + "
    successful_rxns = guess_rxn(reactants)
    rxns = [delim.join(x["names"]) for x in list(successful_rxns.values())]

    fzf = FzfPrompt()
    chosen_rxn = fzf.prompt(rxns, f'--prompt="Reaction {line}: "')[0].split(delim)

    for val in successful_rxns.values():
        if chosen_rxn in val.values():
            # TODO: also return the rxns, so that the csv can be updated
            return val["smarts"]


def fix_rxn(line: str, x: list):

    row = df["Line"] == line
    bad_rxn = df.loc[row, "SMILES"].values.item().split(">>")[0]
    print(bad_rxn)
    N = df.loc[row, "N"].values.item()

    # if not smarts:
    smarts = get_smarts(bad_rxn, line)
    # TODO: currently, looping is not much use if same smarts is always used
    # use a try/except?

    fixed_rxn = use_known_smarts(bad_rxn, smarts)

    df.at[row, "SMILES"] = fixed_rxn
    imgpath = f"reactions/{line}_{N}.png"
    draw_mol(fixed_rxn, imgpath)

    print(f"Fixed reaction {line}")
    x.append(line)


def main():
    """
    Docstring
    """

    file = sys.argv[1]

    # TODO: argparse with int range https://stackoverflow.com/a/6512463
    # start = int(sys.argv[2])
    # if len(sys.argv) == 4:
    #     end = int(sys.argv[3])
    # else:
    #     end = start
    # lines = range(start, end + 1)
    # # sys.exit()

    # for line in lines:

    global df
    df = pd.read_csv(file)
    x = []

    while True:

        try:
            lines = input("Reaction number: ")

            if not lines:
                lines = x[-1] + 1

            elif "," in lines:
                lines = lines.split(",")

            elif "-" in lines:
                start = int(lines.split("-")[0])
                end = int(lines.split("-")[1]) + 1
                lines = range(start, end)
                print(lines)

            else:
                lines = [lines]

            lines = [int(line) for line in lines]

            for line in lines:
                fix_rxn(line, x)

        except KeyboardInterrupt:
            df.to_csv(file, index=False)
            print(f"Wrote: {file}")
            sys.exit()

    # df.to_csv(file, index=False)


if __name__ == "__main__":
    main()
