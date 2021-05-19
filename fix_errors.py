#!/usr/bin/env python3
# import csv
# import requests
# import urllib
from pyfzf.pyfzf import FzfPrompt
import argparse
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


def fix_rxn(column: str, line: str, x: list):

    # row = df["Line"] == line
    row = df[column] == line

    # if column != "Line":
    line = df.loc[row, "Line"].values.item()
    N = df.loc[row, "N"].values.item()
    # row.loc["N"]?

    bad_rxn = df.loc[row, "SMILES"].values.item()
    draw_mol(bad_rxn, show=True)
    print(bad_rxn)

    bad_rxn = bad_rxn.split(">>")[0]

    # if not smarts:
    smarts = get_smarts(bad_rxn, line)
    # TODO: currently, looping is not much use if same smarts is always used
    # use a try/except?

    fixed_rxn = use_known_smarts(bad_rxn, smarts)

    if "c" in fixed_rxn:
        force_fixed_rxn = input_with_prefill(
            "Unusual aromatic? Fix required: ", fixed_rxn
        )
        fixed_rxn = force_fixed_rxn

    df.at[row, "SMILES"] = fixed_rxn
    imgpath = f"reactions/{line}_{N}.png"
    draw_mol(fixed_rxn, imgpath)

    print(f"Fixed reaction {line}")
    # sys.exit()
    x.append(line)


def input_with_prefill(prompt: str, text: str) -> str:
    # https://stackoverflow.com/a/8505387
    def hook():
        readline.insert_text(text)
        readline.redisplay()

    readline.set_pre_input_hook(hook)
    result = input(prompt)
    readline.set_pre_input_hook()
    return result


def main():
    """
    Docstring
    """

    parser = argparse.ArgumentParser(
        description=main.__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        # formatter_class=SortingHelpFormatter,
    )

    parser.add_argument(
        "file",
        help="csv file",
    )

    parser.add_argument(
        "-c",
        "--column",
        # const="col",
        # action="store_const",
        default="Line",
        # nargs="+",
        help="name of column containing SMILES to be processed",
    )
    parser.add_argument(
        "-i",
        "--input-file",
        # const="list",
        # action="store_const",
        help="read numbers from a file",
    )

    args = parser.parse_args()

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
    df = pd.read_csv(args.file)
    print("Loaded", args.file)
    x = []

    while True:

        try:

            if args.input_file:
                with open(args.input_file) as f:
                    lines = [line.rstrip() for line in f if not line.startswith("#")]
                # print(lines)

            else:

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
                fix_rxn(args.column, line, x)

        # except KeyboardInterrupt:
        except:
            df.to_csv(args.file, index=False)
            print(f"Wrote: {args.file}")
            sys.exit()

    # df.to_csv(file, index=False)


if __name__ == "__main__":
    main()
