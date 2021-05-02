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


def main():
    """
    Docstring
    """

    file = sys.argv[1]

    # TODO: int range https://stackoverflow.com/a/6512463
    line = int(sys.argv[2])

    df = pd.read_csv(file)
    row = df["Line"] == line
    # row = df.loc[df["Line"] == line]
    bad_rxn = df.loc[row, "SMILES"].values.item().split(">>")[0]
    N = df.loc[row, "N"].values.item()

    smarts = get_smarts(bad_rxn)
    fixed_rxn = use_known_smarts(bad_rxn, smarts)

    df.at[row, "SMILES"] = fixed_rxn
    df.to_csv("test.csv", index=False)
    imgpath = f"reactions/{line}_{N}.png"
    draw_mol(fixed_rxn, imgpath)

    # sys.exit()


if __name__ == "__main__":
    main()
