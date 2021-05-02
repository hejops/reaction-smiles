#!/bin/env python3
# from pathlib import Path
# import matplotlib.image as mpimg	    much slower than PIL
# import matplotlib.pyplot as plt
# import sys
from PIL import Image
from PIL import ImageDraw
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
import psutil


def close_img_viewer():
    # awfully roundabout way to close the window but ok
    # ironically, closing the image feels slower than opening it; this probably
    # depends on the image viewer used
    # https://stackoverflow.com/a/17489464
    for proc in psutil.process_iter():
        if proc.name() == "display":
            proc.kill()


def show_img(imgpath):
    close_img_viewer()
    # TODO: retain focus on parent?
    IMG = Image.open(imgpath)
    IMG.show()


def draw_mol(smiles, imgpath="structure.png", width=800, height=200, show=False):
    """
    Generate structure from SMILES
    """
    # addAtomIndices -- https://www.rdkit.org/docs/GettingStartedInPython.html#drawing-molecules
    # includeAtomNumbers? -- https://stackoverflow.com/a/53581867
    # example 5? https://www.programcreek.com/python/example/119298/rdkit.Chem.Draw.MolToImage
    # https://www.rdkit.org/docs/source/rdkit.Chem.Draw.rdMolDraw2D.html#rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions.addAtomIndices

    # mol2 = Chem.MolFromSmiles(reactant, sanitize=False)
    # mol = Chem.AddHs(mol)	    # probably unnecessary
    # TODO: extra arg idx -- include idx in img

    if ">>" not in smiles:

        mol = Chem.MolFromSmiles(smiles)
        d = rdMolDraw2D.MolDraw2DCairo(width,height)  # or MolDraw2DSVG to get SVGs
        # d.drawOptions().addAtomIndices = True
        d.DrawMolecule(mol)

    else:

        # https://www.rdkit.org/docs/GettingStartedInPython.html#drawing-chemical-reactions
        rxn = AllChem.ReactionFromSmarts(smiles, useSmiles=True)
        d = rdMolDraw2D.MolDraw2DCairo(width,height)
        d.DrawReaction(rxn)

    d.FinishDrawing()

    png = d.GetDrawingText()
    # TODO: annotate png
    # print(png)
    # draw = ImageDraw.Draw(png)
    # draw.text((0, 0),"Sample Text",(255,255,255))

    open(imgpath, "wb+").write(png)

    # d.WriteDrawingText(imgpath)

    if show:
        show_img(imgpath)

    # return mol


def show_all_products(products_tuple):
    close_img_viewer()
    # all must be sanitized
    img = Draw.MolsToGridImage([p[0] for p in products_tuple], molsPerRow=6)
    img.show()


def show_substruct_matches(reactants, sub):
    # https://www.rdkit.org/docs/GettingStartedInPython.html#drawing-molecules
    # p = Chem.MolFromSmiles(reactants)
    # subms = [x for x in ms if x.HasSubstructMatch(p)]

    close_img_viewer()
    mol = Chem.MolFromSmiles(reactants)
    patt = Chem.MolFromSmarts(sub)
    hit_ats = list(mol.GetSubstructMatch(patt))
    hit_bonds = []
    for bond in patt.GetBonds():
        aid1 = hit_ats[bond.GetBeginAtomIdx()]
        aid2 = hit_ats[bond.GetEndAtomIdx()]
        hit_bonds.append(mol.GetBondBetweenAtoms(aid1, aid2).GetIdx())
    d = rdMolDraw2D.MolDraw2DCairo(500, 500)  # or MolDraw2DCairo to get PNGs
    rdMolDraw2D.PrepareAndDrawMolecule(
        d, mol, highlightAtoms=hit_ats, highlightBonds=hit_bonds
    )
    # d.drawMolecule
    d.FinishDrawing()
    d.WriteDrawingText("atom_annotation_1.png")
    i = Image.open("atom_annotation_1.png")
    i.show()


# def main():

#     # sname = sys.argv[0]
#     if len(sys.argv) != 2:
#         print(draw_mol.__doc__)
#         sys.exit()

#     draw_mol()


# if __name__ == "__main__":
#     main()
