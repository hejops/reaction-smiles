#!/usr/bin/env python3
# from pathlib import Path
# from rdkit.Chem import AllChem
# from rdkit.Chem import Draw
# from rdkit.Chem.Draw import rdMolDraw2D
# import psutil
# import re
# import requests
# import urllib
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import rdChemReactions
import logging
import shutil
import sys

# from drawops import *
# from molops import *
from transforms import *

logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)

# suppress annoying sanitization warnings
# https://github.com/rdkit/rdkit/issues/2683#issuecomment-538872880
RDLogger.DisableLog('rdApp.*')

dashes = "-" * shutil.get_terminal_size().columns

# https://stackoverflow.com/q/55613225
# https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html
# https://www.rdkit.org/docs/GettingStartedInPython.html#chemical-reactions
# https://www.rdkit.org/docs/GettingStartedInPython.html#substructure-based-transformations
# https://www.rdkit.org/docs/RDKit_Book.html#reaction-smarts
# https://www.rdkit.org/docs/RDKit_Book.html#smarts-reference
# https://www.rdkit.org/docs/source/rdkit.Chem.rdChemReactions.html#rdkit.Chem.rdChemReactions.ChemicalReaction
# https://www.rdkit.org/docs/source/rdkit.Chem.rdChemReactions.html#rdkit.Chem.rdChemReactions.HasReactionSubstructMatch

# https://stackoverflow.com/a/59394291
# TODO: fix hinting? probably rdkit's fault

def guess_rxn(reactants: str) -> dict:
    """
    Arguments:
        reactants:	reactants (string: A.B)

    Attempt to apply every stored smarts to the reactants

    Returns all successful reactions as a nested dict:

        {
            product (rdkit tuple object): {
                "smarts": transformation smarts (str),
                "names": names (str),
            }
        }

    Note: An exit will be called if no suitable transformation is found
    """

    def get_product_tuple(reactant_tuple, names):
        # The reactant tuple must be reversed because smarts order is
        # important!
        product_tuple_1 = reaction.RunReactants(reactant_tuple)
        product_tuple_2 = reaction.RunReactants(reactant_tuple[::-1])

        if len(product_tuple_1) > 0:
            return product_tuple_1

        elif len(product_tuple_2) > 0:
            names = names[::-1]
            return product_tuple_2

    reac1 = reactants.split(".")[0]
    reac2 = reactants.split(".")[1]
    reactant_tuple = (Chem.MolFromSmiles(reac1), Chem.MolFromSmiles(reac2))

    successful_rxns = {}

    for smarts in transformations:

        # useSmiles=True leads to mangled products; don't use it!
        names = transformations[smarts]
        reaction = rdChemReactions.ReactionFromSmarts(smarts, useSmiles=False)
        product_tuple = get_product_tuple(reactant_tuple, names)
        # names = " + ".join(names)

        if product_tuple:
            successful_rxns[product_tuple] = { "smarts": smarts, "names": names }

        else:
            continue

    if not successful_rxns:
        print("No matching transformations found")
        print(reactants)
        sys.exit()

    return successful_rxns


def get_product_smiles(product) -> tuple[str, bool]:
    """
    Arguments:
        product (mol)

    Returns:
        product smiles (string)
        sanitization state (boolean)
    """

    # sanitize "partially"
    # http://www.qsar4u.com/files/rdkit_tutorial/rdkit.html#Structure-Sanitization
    # http://www.rdkit.org/docs/RDKit_Book.html#molecular-sanitization
    # https://www.rdkit.org/docs/Cookbook.html#explicit-valence-error-partial-sanitization
    # pdt.UpdatePropertyCache(strict=False)

    try:
        Chem.SanitizeMol(
            product,
            Chem.SanitizeFlags.SANITIZE_FINDRADICALS
            | Chem.SanitizeFlags.SANITIZE_KEKULIZE
            | Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
            | Chem.SanitizeFlags.SANITIZE_SETCONJUGATION
            | Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION
            | Chem.SanitizeFlags.SANITIZE_SYMMRINGS,
            # catchErrors=True,
        )
        product_smiles = Chem.MolToSmiles(product, kekuleSmiles=True)
        sanit = True

    except:
        product_smiles = Chem.MolToSmiles(product, kekuleSmiles=False)
        sanit = False

    return product_smiles, sanit


def get_full_reaction(reactants: str) -> tuple[str, str, bool]:
    """
    Arguments:
        reactants (string)

    Returns:
        full smiles (string)
        transformation name (string)
        sanitization state (boolean)
    """

    successful_rxns = guess_rxn(reactants)
    print(f"{len(successful_rxns)} matching transformations found")

    # Try to sanitize every product obtained. Note: smiles are NOT canonical by
    # default, which typically leads to an "inflated" dict.

    # { full_smiles : smarts }
    good_smiles = {}
    bad_smiles = {}

    # https://stackoverflow.com/a/53865188
    # ith reaction, jth product
    for i, (product_tuple, attempted_transformation) in enumerate(successful_rxns.items()):

        x = [p[0] for p in product_tuple]
        for j, product in enumerate(x):

            product_smiles, sanit = get_product_smiles(product)

            if not product_smiles:
                # If even an unsanitized product cannot be obtained, exit.
                # A new nucleophile or electrophile smarts must then be added.
                print(dashes)
                print("Products could not be obtained")
                print(f"Reactants: {reactants}")
                print("Last product:")
                print(get_product_smiles(product)[0].split(">")[-1])
                print(dashes)
                sys.exit()

            # TODO: the "smarts" key is not of much use, return it instead of sanit
            current_reaction = attempted_transformation["names"]

            if not sanit:
                full_smiles = f"{reactants}>>{product_smiles}"
                bad_smiles[product_smiles] = current_reaction
                continue

            full_smiles = f"{reactants}>>{product_smiles}"
            # print(get_product_smiles(product_tuple[0][1])[0])
        
            # TODO: with and (instead of or), products like Br- cannot be retrieved
            if attempted_transformation["smarts"] in two_products and len(product_tuple) > 1:
                ps2, _ = get_product_smiles(product_tuple[0][1])
                # print(ps2)
                full_smiles += f".{ps2}"

            good_smiles[full_smiles] = current_reaction

            if j + 1 == len(x):
                print(f"{i+1}: {current_reaction}, {j} products")

            # stop after 20 products obtained
            if len(good_smiles) == 20:
                break

    # logging.info(good_smiles)
    # sys.exit()

    if good_smiles:
        # TODO: remove duplicate keys in dict
        print(f"Up to {len(good_smiles)} good products found")

        """
        The most likely product is the first one that passed the smarts match
        and product sanitization.
    
        Return the first key (full smiles) of the dict, the transformation
        name, and sanitization state; transformation smarts is (probably) not
        needed from this point on.
    
        list(dict) is required because dicts cannot be hashed by index;
        set(dict) should not be used because it does not preserve order of
        addition.
        """

        result = list(good_smiles)[0]
        print(result)
        print(good_smiles[result])
        # print(len(list(good_smiles)))
        # print(len(set(good_smiles)))
        # sys.exit()
        return result, good_smiles[result], True

    elif bad_smiles:
        print("Warning: product not sanitized")

        result = list(bad_smiles)[0]
        print(result)
        print(bad_smiles[result])
        return result, bad_smiles[result], False

def main():
    """
    Arguments:
        SMILES [A.B]

    Generate a product from reactants based on a set of transformation rules.
    """

    if len(sys.argv) == 1:
        print(main.__doc__)
        sys.exit()

    if len(sys.argv) == 2:

        reactants = sys.argv[1]
        if ">" in reactants:
            print("Reaction detected (>>); product will be removed")
            reactants = reactants.split(">")[0]

        get_full_reaction(reactants)

if __name__ == "__main__":
    main()
