#!/usr/bin/env python3
# import os
# import requests
# import shutil
# import urllib
# import logging
# import sys
from rdkit import Chem

# logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)
# https://www.rdkit.org/docs/GettingStartedInPython.html#modifying-molecules
# https://www.rdkit.org/docs/GettingStartedInPython.html#writing-molecules


def map_atom(mol, idx):
    # print("B", mol, idx)
    # print("C", mapped_mol, idx)
    # return mol.GetAtomWithIdx(idx).SetAtomMapNum(1), 555
    return mol.GetAtomWithIdx(idx).SetAtomMapNum(1)


def auto_map_charged(mol):
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() == 1:
            atom.SetAtomMapNum(0)
    return mol


def unmap(mol):
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return mol


def remap(mol, idx):
    unmap(mol)
    map_atom(mol, idx)
    return mol


def neutralise(mol):
    # https://www.rdkit.org/docs/Cookbook.html?highlight=getformalcharge#neutralizing-molecules
    for atom in mol.GetAtoms():
        # atom.SetNoImplicit(False)
        atom.SetNumExplicitHs(0)
        atom.SetFormalCharge(0)
        # atom.SetNumRadicalElectrons(0)
        atom.UpdatePropertyCache()
    return mol
