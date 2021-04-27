from smartjoin import get_full_rxn

def test_difficult_enamine():
    rxn = "CC1=CC(C)=CN1.C1CC2=C3N1CCCC3=CC([CH+]C1=CC3=C4N(CC3)CCCC4=C1)=C2"
    pdt = "<++>"
    try:
        p = get_full_rxn(rxn, 0)
    except:
        assert False
    # assert p == pdt


# def test_<++>():
#     rxn = "<++>"
#     pdt = "<++>"
#     try:
#         p = guess_rxn(rxn, 0)
#     except:
#         assert False
#     # assert p == pdt


# def test_<++>():
#     rxn = "<++>"
#     pdt = "<++>"
#     try:
#         p = guess_rxn(rxn, 0)
#     except:
#         assert False
#     # assert p == pdt


# def test_<++>():
#     rxn = "<++>"
#     pdt = "<++>"
#     try:
#         p = guess_rxn(rxn, 0)
#     except:
#         assert False
#     # assert p == pdt


def test_difficult_aromatic_electrophile():
    rxn = "C1CCC(=CC1)N1CCOCC1.[O-][N+]1=C2C=C(C=C(C2=NO1)[N+]([O-])=O)[N+]([O-])=O"
    pdt = "C1CCC(=CC1)N1CCOCC1.[O-][N+]1=C2C=C(C=C(C2=NO1)[N+]([O-])=O)[N+]([O-])=O>>O=[N+]([O-])C1=C(C2=C(N3CCOCC3)CCCC2)C2=[N+]([O-])ON=C2C([N+](=O)[O-])=C1"
    try:
        p = get_full_rxn(rxn, 0)
    except:
        assert False
    assert p == pdt


def test_amine_aniline():
    rxn = "CC(C)N.CN(C)C1=CC=C(C=C2C=CC(C=C2)=[N+](C)C)C=C1"
    pdt = "CC(C)N.CN(C)C1=CC=C(C=C2C=CC(C=C2)=[N+](C)C)C=C1>>CC(C)[NH2+]C(C1=CC=C(N(C)C)C=C1)C1=CC=C(N(C)C)C=C1"
    try:
        p = get_full_rxn(rxn, 0)
    except:
        assert False
    assert p == pdt


def test_amide_phenol():
    rxn = "O=C1NC=CC=C1.COC1=CC=C(C=C2C=C(C(=O)C(=C2)C2=CC=CC=C2)C2=CC=CC=C2)C=C1"
    pdt = "O=C1NC=CC=C1.COC1=CC=C(C=C2C=C(C(=O)C(=C2)C2=CC=CC=C2)C2=CC=CC=C2)C=C1>>COC1=CC=C(C(OC2=[NH+]C=CC=C2)C2=CC(C3=CC=CC=C3)=C([O-])C(C3=CC=CC=C3)=C2)C=C1"
    try:
        p = get_full_rxn(rxn, 0)
    except:
        assert False
    assert p == pdt
