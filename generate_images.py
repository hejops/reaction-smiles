#!/usr/bin/env python3
# from rdkit import Chem
# import os
# import requests
# import shutil
# import urllib
from pathlib import Path
import csv
import logging
import sys

from drawops import *

logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)


def main(db):
    """
    Docstring
    """

    reader = csv.DictReader(open(db))
    # t = sum(1 for row in open(db))
    Path("reactions").mkdir(parents=True, exist_ok=True)

    for i, row in enumerate(reader):
        N = row["N"]
        reaction = row["SMILES"]
        reacs = ", ".join(sorted([row["Reactant 1"], row["Reactant 2"]]))
        Path(f"reactions/{reacs}").mkdir(parents=True, exist_ok=True)
        imgpath = f"reactions/{reacs}/{i+1}_{N}.png"
        draw_mol(reaction, imgpath)

        print(f"\rGenerated {i+1} images...", end="")
    print("Done")


if __name__ == "__main__":
    db = sys.argv[1]
    main(db)
