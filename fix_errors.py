#!/usr/bin/env python3
# import requests
# import urllib
from pyfzf.pyfzf import FzfPrompt
import csv
import logging
import os
import pandas as pd
import shutil
import sys

logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)

from smartjoin import *
from drawops import *


def get_smarts(reactants):
    delim = " + "
    successful_rxns = guess_rxn(reactants)
    rxns = [delim.join(x["names"]) for x in list(successful_rxns.values())]

    fzf = FzfPrompt()
    chosen_rxn = fzf.prompt(rxns, "--reverse")[0].split(delim)

    for val in successful_rxns.values():
        if chosen_rxn in val.values():
            # TODO: also return the rxns, so that the csv can be updated
            return val["smarts"]


def main(smarts=None):
    """
    Docstring
    """

    file = sys.argv[1]

    # TODO: argparse with int range https://stackoverflow.com/a/6512463
    start = int(sys.argv[2])

    if len(sys.argv) == 4:
        end = int(sys.argv[3])
    else:
        end = start

    lines = range(start, end + 1)
    # sys.exit()

    df = pd.read_csv(file)

    for line in lines:

        row = df["Line"] == line
        bad_rxn = df.loc[row, "SMILES"].values.item().split(">>")[0]
        N = df.loc[row, "N"].values.item()

        if not smarts:
            smarts = get_smarts(bad_rxn)

        fixed_rxn = use_known_smarts(bad_rxn, smarts)
        # TODO: exit if wrong smarts used (i.e. overshot)

        df.at[row, "SMILES"] = fixed_rxn
        # TODO: once this is robust, remember to reuse the filename
        df.to_csv("test.csv", index=False)
        imgpath = f"reactions/{line}_{N}.png"
        draw_mol(fixed_rxn, imgpath)

        print(f"Fixed reaction {line}")


if __name__ == "__main__":
    main()
