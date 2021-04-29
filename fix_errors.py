#!/usr/bin/env python3
# import requests
# import urllib
from pyfzf.pyfzf import FzfPrompt
import csv
import logging
import os
import shutil
import sys

logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)

from smartjoin import *

def get_smarts(reactants):
    delim = " + "
    successful_rxns = guess_rxn(reactants)
    rxns = [delim.join(x["names"]) for x in list(successful_rxns.values())]

    fzf = FzfPrompt()
    chosen_rxn = fzf.prompt(rxns, '--reverse')[0].split(delim)

    for val in successful_rxns.values():
        if chosen_rxn in val.values():
            return val["smarts"]

def main():
    """
    Docstring
    """

    file = sys.argv[1]

    # TODO: int range https://stackoverflow.com/a/6512463
    line = sys.argv[2]
    # smarts = sys.argv[3]

    f = open(file)
    reader = csv.reader(f)
    row = [x for x in reader if x[0] == line][0]
    pdt_to_fix = [x.split(">>")[0] for x in row if ">>" in x][0]

    smarts = get_smarts(pdt_to_fix)
    use_known_smarts(pdt_to_fix, smarts)

if __name__ == "__main__":
    main()

