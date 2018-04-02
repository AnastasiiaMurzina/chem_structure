from biopandas.mol2 import PandasMol2

if __name__ == '__main__':
    pmol = PandasMol2().read_mol2('benzene.mol2')
    print(pmol)

    '''Notes about library:
    pmol = PandasMol2().read_mol2('benzene.mol2')
    pmol.mol2_text[:500] - text of file
    pmol.df.head(n) - first n (or tail(n) - lasts) atoms of file: xyz, and info (3 x 9 table + charge)
    pmol.df['atom_type'] != 'H'] - filter in []
    pmol.df[pmol.df['atom_type'] == 'O.2'][['x', 'y', 'z']] - filers example
    
    
    Typos:
    Code       Definition
C.3        carbon sp3
C.2        carbon sp2
C.1        carbon sp
C.ar       carbon aromatic
C.cat      cabocation (C+) used only in a guadinium group
N.3        nitrogen sp3
N.2        nitrogen sp2
N.1        nitrogen sp
N.ar       nitrogen aromatic
N.am       nitrogen amide
N.pl3      nitrogen trigonal planar
N.4        nitrogen sp3 positively charged
O.3        oxygen sp3
O.2        oxygen sp2
O.co2      oxygen in carboxylate and phosphate groups
O.spc      oxygen in Single Point Charge (SPC) water model
O.t3p      oxygen in Transferable Intermolecular Potential (TIP3P) water model
S.3        sulfur sp3
S.2        sulfur sp2
S.O        sulfoxide sulfur
S.O2/S.o2  sulfone sulfur
P.3        phosphorous sp3
F          fluorine
H          hydrogen
H.spc      hydrogen in Single Point Charge (SPC) water model
H.t3p      hydrogen in Transferable Intermolecular Potential (TIP3P) water model
LP         lone pair
Du         dummy atom
Du.C       dummy carbon
Any        any atom
Hal        halogen
Het        heteroatom = N, O, S, P
Hev        heavy atom (non hydrogen)
Li         lithium
Na         sodium
Mg         magnesium
Al         aluminum
Si         silicon
K          potassium
Ca         calcium
Cr.thm     chromium (tetrahedral)
Cr.oh      chromium (octahedral)
Mn         manganese
Fe         iron
Co.oh      cobalt (octahedral)
Cu         copper
    '''