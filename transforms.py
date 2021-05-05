#!/usr/bin/env python3
# import itertools
# import re
# import sys

# Construct a dictionary of reaction SMARTS to be iterated through

# https://github.com/rdkit/rdkit/discussions/4042#discussioncomment-621064
# https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html
# https://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html
# TODO: dump/load json?

# TODO:
# CC1=CC=C(C=C1)S(=O)(=O)[CH-]C1=CC=CC=C1.COC1=CC=C(C=C2C3=CC=CC=C3C(=O)C3=C2C=CC=C3)C=C1>>COC1=C(C(C2=CC=CC=C2)S(=O)(=O)C2=CC=C(C)C=C2)C=C(C=C2C3=CC=CC=C3C(=O)C3=C2C=CC=C3)C=C1

# C1CCC(=CC1)N1CCOCC1.ClC1=CC(Cl)(Cl)C(=O)C2=C1C=CC=N2

nucleophiles = {
    # these are not really needed since the alkene fallback works fine
    # "SnCC=C": ["[Sn:1][C:2][C:3]=[CX3:4]", "[Sn:1][C:2][C+:3][CX3:4]"],
    "SiCC=C": ["[Si:1][C:2][C:3]=[C:4]", "[Si:1][C:2][C+:3][C:4]"],
    # "silyl, from O": ["[Si:1][O:2][C:3]=[C:4]", "[Si:1][O+:2]=[C:3][C:4]"],
    "silyl, SiOC=C": ["[Si:1][O:2][C:3]=[CX3:4]", "[Si:1][O:2][C+:3][CX3:4]"],
    # "silyl, sanit": ["[Si:1][O:2][C:3]=[C:4]", "[Si:1][O:2][C:3][C:4]"],
    "enolate": ["[C:1]=[C:2][O-:3]", "[C:1]=[C:2][O:3]"],
    "malonic acid: O=CCC=O": [
        "[O:1]=[C:2][CH:3]([C:4]=[O:5])",
        "[O:1]=[C:2][C+:3]([C:4]=[O:5])",
    ],
    "carboxylate": ["[O:1]=[C:2][O-:3]", "[O:1]=[C:2][O+0:3]"],
    "azide, CNN-/NNN-": ["[N+:2]=[N-:3]", "[N+:2]=[N+0:3]"],
    # "azide, NNC-": ["[N:1]#[N:2]=[CH:3]", "[N:1]#[N+:2][CH:3]"],
    "amide": ["[nH:1][c:2]=[O:3]", "[nH+:1]=[c:2][O:3]"],
    # "imidazole, H": ["[nH:1][c:2][n:3]", "[nH+:1][c:2][n:3]"],
    "imidazole": ["[n:1][c:2][n:3]", "[n+1:1][c:2][n:3]"],
    # "imidazole, de-arom": ["[n:1][c:2][n:3]", "[N+1:1]=[C:2][N:3]"],
    "pyridine": ["[n:1]([c:2])([c:3])", "[n+:1]([c:2])([c:3])"],
    "enamine": ["[N:1][C:2]=[C:3]", "[N+:1]=[C:2][C:3]"],
    "enamine, aromatic": ["[nX3:1][c:2][cH:3]", "[NX3+:1]=[C:2][CH:3]"],
    "NC=N": ["[N:1][C:2]=[N:3]", "[N+:1]=[C:2][N:3]"],
    "nitrile": ["[*+0:1][C:2]#[N:3]", "[*+0:1][C+:2]=[N:3]"],
    "S ylide": ["[S:1]=[CX3:2]", "[S+:1][CX3:2]"],
    "S ylide, N": ["[S:1][NH-:2]", "[S:1][NH+0:2]"],
    "S ylide, disconnect?": ["[S+:1][CX3:2]", "[S+0:1].[CX3:2]"],
    "alkene, C=CH2 or C=CH": ["[C:1]=[CX3;H2,H1:2]", "[C+:1][CX4;H2,H1:2]"],
    "alkene, terminal": ["[C:1]=[CX3H2:2]", "[C+:1][CX4H2:2]"],
    "alkene": ["[C:1]=[CX3:2]", "[C+:1][CX3:2]"],
    "alkene (dearomatise)": ["[cH;c:1][cH1:2]", "[CH+;c+:1][CH:2]"],
    "N=": ["[C:1]=[NX2:2]", "[C+:1][N:2]"],
    "C-": ["[C-:1]", "[C+0:1]"],
    # "C-": ["[CX3-:1]", "[CX3+0:1]"],
    "P": ["[P:1]", "[P+1:1]"],
    # "aromatic NH": ["[nH:1]", "[n+:1]"],
    # "aromatic N, sanit": ["[n:1]", "[n:1]"],
    # https://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html#N
    "amine (not amide!)": ["[NX3;H2,H1;!$(NC=O);!$(NC=S):1]", "[NX3+1:1]"],
    # "amine (pri/sec)": ["[NX3;H2,H1:1]", "[NX3;H2,H1+:1]"],
    "amine (pri)": ["[NH2:1]", "[NH2+:1]"],
    "amine (sec)": ["[NH1:1]", "[NH1+:1]"],
    "amine (any)": ["[NX3:1]", "[NX3+:1]"],
    "alkoxide": ["[O-:1]", "[O+0:1]"],
    "alcohol": ["[OH:1]", "[OH+:1]"],
    "water": ["[OH2:1]", "[OH2+:1]"],
    "aryl": ["[cH:1]", "[c+0:1]"],
    "HF": ["[FH:1]", "[FH+:1]"],
    "NH-": ["[NH-:1]", "[NH+0:1]"],
    "furan": ["[o:1][c:2][cH:3]", "[O+:1]=[C:2][CH:3]"],
    "thiol": ["[SH:1]", "[SH+:1]"],
    "P=CH (HWE, neutral)": ["[P:1]=[CX3:2]", "[P+:1][CX3:2]"],
    "thiofuran": ["[s:1]", "[S+:1]"],
    "NHC": ["[CX2H0:1]([N:2])([N:3])", "[C:1]([N:2])([N:3])"],
    # "<++>": ["<++>", "<++>"],
    # "<++>": ["<++>", "<++>"],
    # "<++>": ["<++>", "<++>"],
}

# atom indices are incremented to allow simple concatenation
electrophiles = {
    "aniline": [
        "[C:1]=[C:2][C:3]=[C:4][C:5]=[N+1:6]",
        "[C:1][C:2]=[C:3][C:4]=[C:5][N+0:6]",
    ],
    # C[Si](C)(C)OC1=CCCC1.CN(C)C1=CC2=C(C=C1)C=C1C=CC(C=C1O2)=[N+](C)C
    "aniline, aromatic": [
        "[cH:1][c:2][c:3][c:4][c:5]=[N+1:6]",
        "[cH:1][c:2][c:3][c:4][c:5][N+0:6]",
    ],
    "phenol": [
        "[C:1]=[C:2][C:3]=[C:4][C:5]=[O:6]",
        "[C:1][C:2]=[C:3][C:4]=[C:5][O-:6]",
    ],
    "phenol, semi-aromatic": [
        "[C:1]=[C:2][c:3][c:4][C:5]=[O:6]",
        "[C:1][C:2][c:3][c:4][C:5][O-:6]",
    ],
    "nitro vinyl": [
        "[C:1]=[C:2][N+:3]([O-:4])=[O:5]",
        "[C:1][C:2]=[N+:3]([O-:4])[O-:5]",
    ],
    "nitrile vinyl": ["[C:1]=[C:2][C:3]#[N:4]", "[C:1][C:2]=[C:3]=[N-:4]"],
    # see https://github.com/rdkit/rdkit/discussions/4083#discussioncomment-655125
    # COC1=NC2=NO[N+]([O-])=C2C=N1.C1CCC(=CC1)N1CCOCC1
    "c=cc=o+": ["[cH:1][c:2][c:3][o+1:4]", "[cH:1][c:2][c:3][o+0:4]"],
    "c=cc=n+": ["[cH:1][c:2][c:3][n+1:4]", "[cH:1][c:2][c:3][n+0:4]"],
    "c=cc=n (?)": ["[cH:1][c:2][c:3][n:4]", "[cH:1][c:2][c:3][n-:4]"],
    "n=cc=n+": ["[n:1][c:2][c:3][n+1:4]", "[n:1][c:2][c:3][n+0:4]"],
    "enone (aromatic)": ["[cH:1][c:2][c:3]=[O:4]", "[cH:1][c:2][c:3][O-:4]"],
    "azodicarb": ["[N:1]=[N:2][C:3]=[O:4]", "[N:1][N:2]=[C:3][O-:4]"],
    "enone, C=CC=O": ["[C:1]=[C:2][C:3]=[#8:4]", "[C:1][C:2]=[C:3][#8-:4]"],
    # TODO
    # C[Si](C)(C)OC1=CCCO1.CN1C(=O)[C@H](CC2=CNC3=C2C=CC=C3)\[N+](=C\C=C\C2=CC=CC=C2)C1(C)C
    # C[Si](C)(C)OC1=CCCO1.CN1C(=O)[C@H](CC2=CNC3=C2C=CC=C3)\N(C=CC(C)\C2=CC=CC=C2)C1(C)C
    "enone, C=CC=N+": [
        "[CH:1]=[C:2][C:3]=[NX3+:4]",
        "[#6H:1][#6:2]=[#6:3][#7X3+0:4]",
    ],
    # "enone (n)": ["[c:1][c:2][c:3][N+:4]", "[c:1][c:2][c:3][n:4]"],
    "C+ vinyl": ["[C+:1][C:2]=[C:3]", "[C+0:1][C:2]=[C:3]"],
    "C+ (2 aromatic neighbours)": ["[C+:1]([c:2])[c:3]", "[C+0:1]([c:2])[c:3]"],
    "C+ (1 aromatic neighbour)": ["[c:1][C+:2][C:3]", "[C+0:2]([c:1])[C:3]"],
    "NH=NH+": ["[NH:1]=[NH+:2]", "[NH:1][NH+0:2]"],
    "aryl+": ["[cH+:1]", "[cH+0:1]"],
    # this is more sensible, but also gives false positives
    "aryl (dearomatise)": ["[cH:1][cX3:2]", "[CH:1][C-:2]"],
    "C+": ["[CX3+:1]", "[C+0:1]"],
    "C+, cyclopropene": ["[c+:1]1[c:2][c:3]1", "[c+0:1]1[c:2][c:3]1"],
    # why does this override aniline?
    "C=N+": ["[CX3:1]=[NX3+:2]", "[CX4:1][NX3+0:2]"],
    "c=n+": ["[cH:1][n+:2]", "[CH:1][N+0:2]"],
    # does this match anything?
    "carbonyl": ["[CX3:1]=[O:2]", "[C:1][O-:2]"],
    # C\C=C/C.C1=CC=C(C=C1)[C+]=C=C(C1=CC=CC=C1)C1=CC=CC=C1
    "allene": ["[C+:1]=[C:2]=[C:3]", "[C+0:1]=[C:2]=[C:3]"],
    # C[CH-][N+]([O-])=O.CC1=CC=C(C=C2C=CC=C2)C=C1
    # TODO: this should not be triggered
    # CC(C)N.C1CN2CCCC3=C2C(C1)=CC([CH+]\C=C\C1=CC2=C4N(CCCC4=C1)CCC2)=C3>>CC(C)[NH2+]C([CH-][CH+]C1=CC2=C3C(=C1)CCCN3CCC2)C1=CC2=C3C(=C1)CCCN3CCC2
    "alkene, cC=C": ["[CH:1]([c:2])=[C:3]", "[CH:1]([c:2])[C-:3]"],
    "sulfate, C=CS=O": ["[C:1]=[C:2][S:3]=[O:4]", "[C:1][C:2]=[S:3][O-:4]"],
    "Br": ["[CX4:1][Br:2]", "[C:1].[Br-:2]"],
    "diazonium": ["[N:1]#[N+:2]", "[N:1]=[N+0:2]"],
    "C(Cl)2": ["[Cl:1][C:2]([Cl:3])[C:4]=[O:5]", "[Cl:1].[c:2]([Cl:3])[c:4][O-:5]"],
    "SMe+": ["[CH3:1][s+1:2]", "[CH3:1].[s+0:2]"],
    # "<++>": ["<++>", "<++>"],
    # "<++>": ["<++>", "<++>"],
}

# edge cases
two_products = {
    # C[Si](C)(C)OC1=CCCO1.ClC1=CC(Cl)(Cl)C(=O)C2=C1C=CC=N2
    # C1CCC(=CC1)N1CCOCC1.ClC1=CC(Cl)(Cl)C(=O)C2=C1C=CC=N2
    "[Si:1][O:2][C:3]=[C:4].[C:5]([Cl:6])([Cl:7])[C:8]=[O:9]>>[Si:1][O:2][C+:3][C:4][Cl:6].[C:5]([Cl:7])=[C:8][O-:9]": [
        "silyl",
        "C(Cl)2",
    ],
    "[N:1][C:2]=[CX3:3].[C:5]([Cl:6])([Cl:7])[C:8]=[O:9]>>[N+:1]=[C:2][C:3][Cl:6].[C:5]([Cl:7])=[C:8][O-:9]": [
        "enamine",
        "C(Cl)2",
    ],
    # TODO: main product is not getting returned in full smiles
    # O.C[S+]1C2=C(C=CC=C2)C2=C1C=CC=C2
    # CO.C[S+]1C2=C(C=CC=C2)C2=C1C=CC=C2
    "[OX2:1].[CH3:2][s+:3]>>[OX2+:1][CH3:2].[s+0:3]": [
        "water",
        "CH3S+",
    ],
    # "[C-:1].[C:2][Br:3]>>[C+0:1][C:2].[Br-:3]": [
    #     "C-",
    #     "Br",
    # ],
    # # CCCN.BrCC1=CC=CC=C1
    # "[NX3:1].[C:4][Br:5]>>[NX3+:1][C:4].[Br-:5]": [
    #     "amine",
    #     "Br",
    # ],
    # enamine ccn
    # BrCC1=CC=CC=C1.N1C=CC2=C1C=CC=C2
}

important_reactants = ["C-:", "P", "S", "s", "H0"]
transformations = {}
important_transformations = {}
ugly_transformations = {}

for e_name, e in electrophiles.items():

    for n_name, n in nucleophiles.items():
        n_reac = n[0]
        n_pdt = n[1]

        # As long as reactants have <= 10 atoms, an actual increment is
        # probably unnecessary...
        e_reac = e[0].replace(":", ":5")
        e_pdt = e[1].replace(":", ":5")
        reacs = f"{n_reac}.{e_reac}"
        smarts = f"{reacs}>>{n_pdt}{e_pdt}"

        # if full < 8 colons -> simple

        # https://stackoverflow.com/a/3389611
        # https://stackoverflow.com/a/6531704
        if any(reac in reacs for reac in important_reactants):
            print(reacs)
            important_transformations[smarts] = [n_name, e_name]

        # elif n_name in important_reactants or e_name in important_reactants:
        #     important_transformations[smarts] = [n_name, e_name]

        # might not be needed
        # if "." in n_pdt or "." in e_pdt:
        #     two_products[smarts] = [n_name, e_name]

        # this is mostly just to give alkene electrophiles low priority
        elif "-" in e_pdt:
            ugly_transformations[smarts] = [n_name, e_name]

        else:
            transformations[smarts] = [n_name, e_name]

"""
Sort the dictionary so that longest smarts are tried first. While this
generally produces the most sensible transformation on the first attempt, some
reactants may end up getting undue precedence over others, which can result in
products that can never be sanitized (an irritating example being pyridine >
imidazole). This is why it is necessary to attempt all possible transformations
later.

Reactants with relatively obvious reactive centres (and short smarts) are
forcibly placed higher in the dictionary.
"""

# TODO: sort the important ones also!

# https://stackoverflow.com/a/11753945
# https://stackoverflow.com/a/53962042
def sort_dict(Dict: dict) -> dict:
    sorted_dict = {}
    for k in sorted(Dict, key=len, reverse=True):
        sorted_dict[k] = Dict[k]
    return sorted_dict


transformations = (
    sort_dict(important_transformations)
    # probably not necessary anymore
    # | two_products
    | sort_dict(transformations)
    | ugly_transformations
)

# TODO: find out which nucleophiles/electrophiles were never used, and remove
# them to reduce iterations

if __name__ == "__main__":
    print(
        f"""
{len(transformations)} transformations generated
    ({len(nucleophiles)} nucleophiles, {len(electrophiles)} electrophiles)

{len(important_transformations)} reactions are marked important
Important reactants: {important_reactants}
"""
    )

# enone (N)
# bad arom??
# C[Si](C)(C)OC1=CCCO1.CN1C(=O)[C@H](CC2=CNC3=C2C=CC=C3)\[N+](=C\C=C\C2=CC=CC=C2)C1(C)C

# aryl nuc must be sanitised, no choice
# CC1=CC(C)=CC=C1.CC1=CC=C([CH+]C2=CC=CC=C2)C=C1

# https://github.com/rdkit/rdkit/discussions/4083
# arom O???
# O.C1=CC2=CC3=CC=CC=C3[O+]=C2C=C1
# product should be
# C1=CC2C([OH2+])C3=CC=CC=C3OC=2C=C1

# malonic
# CCOC(=O)C(C)C(=O)OCC.CC(C)(C)C1=CC(=CC2=CC3=C4N(CCCC4=C2)CCC3)C=C(C1=O)C(C)(C)C

# pyridine
# CN1CCC[C@H]1C1=CN=CC=C1.CN(CC(F)(F)F)C1=CC=C([CH+]C2=CC=C(C=C2)N(C)CC(F)(F)F)C=C1>>CN(CC(F)(F)F)C1=CC=C(C(C2=CC=C(N(C)CC(F)(F)F)C=C2)[N+]2(C)CCC[C@H]2C2=CN=CC=C2)C=C1

# enolate
# CN1CCCC2=CC([CH+]C3=CC=C4N(C)CCCC4=C3)=CC=C12.C\C(=C(\[O-])C1=[N+](CCN1C1=C(C)C=C(C)C=C1C)C1=C(C)C=C(C)C=C1C)C1=CC=CC=C1

# CCCCOP(OCCCC)OCCCC.CN(C)C1=CC=C(C=C1)C1=CC=CC=C[CH+]1
# matching the imidazole reactant is basically impossible if uppercase is used
# CN1C=CN=C1.CN(C)C1=CC=C(C=C2C=CC(C=C2)=[N+](C)C)C=C1
# CC1=NC=CN1.CN(C)C1=CC=C(C=C2C=CC(C=C2)=[N+](C)C)C=C1
# CC1=NC=CN1.C1CC2=C3N1CCCC3=CC([CH+]C1=CC3=C4N(CC3)CCCC4=C1)=C2
# CC1=CNC=N1.CN1CCCC2=CC([CH+]C3=CC=C4N(C)CCCC4=C3)=CC=C12>>CC1=CN(C(C2=CC3=C(C=C2)N(C)CCC3)C2=CC=C3C(=C2)CCCN3C)C=N1
# because nucleophile is somewhat ambiguous (aliph/arom), only generic atom # number is specified
# O=C1NC=CC=C1.COC1=CC=C(C=C2C=C(C(=O)C(=C2)C2=CC=CC=C2)C2=CC=CC=C2)C=C1
# same for this
# O=C1NC=CC=C1.CN1C(=O)N(C)C(=O)C(=CC2=CC3=C4N(CCCC4=C2)CCC3)C1=O
# [O-]C(=O)[C@@H]1CCCN1C1=CCCCC1.CC(C)(C)C1=CC(=CC2=CC=C(C=C2)[N+]([O-])=O)C=C(C1=O)C(C)(C)C
# CC(C)N.COC1=CC=C(C=C2C=C(C(=O)C(=C2)C2=CC=CC=C2)C2=CC=CC=C2)C=C1
# OC1=CC=C(Cl)C=C1.COC1=CC=C(C=C2C=C(C(=O)C(=C2)C2=CC=CC=C2)C2=CC=CC=C2)C=C1
# CC(C)N.CN(C)C1=CC=C(C=C2C=CC(C=C2)=[N+](C)C)C=C1

# [O-][N+](=O)\C=C\C1=CC=CC=C1.[O-]C(=O)[C@@H]1CCCN1C1=CCCCC1
# C1CN2CCCN=C2C1.O=C1C(=CC2=CC=CC=C2)C(=O)C2=C1C=CC=C2
# [O-]C(=O)[C@@H]1CCCN1C1=CCCCC1.CN1CCCC2=CC([CH+]C3=CC=C4N(C)CCCC4=C3)=CC=C12
# CC(C)N.C1CN2CCCC3=C2C(C1)=CC([CH+]\C=C\C1=CC2=C4N(CCCC4=C1)CCC2)=C3
# C1COCCN1.CCOC(=O)\N=N\C(=O)OCC
# NCC1=CC=CC=C1.CCOC(=O)C(=CC1=CC=C(C)C=C1)C(=O)OCC
# simplest transformations are the fallback; this should avoid false positives
# CC1(C)OC(=O)[C-](CC2=CC=CC=C2)C(=O)O1.C1CC2=C3N1CCCC3=CC([CH+]C1=CC3=C4N(CC3)CCCC4=C1)=C2
# CC(C)N.C1CC2=C3N1CCCC3=CC([CH+]C1=CC3=C4N(CC3)CCCC4=C1)=C2
# reactions from here on are "worked up" because attempting to write a
# smarts that retains charges was too difficult
# C1CCC(=CC1)N1CCOCC1.[O-][N+]1=C2C=C(C=C(C2=NO1)[N+]([O-])=O)[N+]([O-])=O
# CC1=CC(C)=CN1.C1CC2=C3N1CCCC3=CC([CH+]C1=CC3=C4N(CC3)CCCC4=C1)=C2
