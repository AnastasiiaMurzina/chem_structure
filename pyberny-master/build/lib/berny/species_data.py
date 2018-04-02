# Any copyright is dedicated to the Public Domain.
# http://creativecommons.org/publicdomain/zero/1.0/


def get_property(idx, name):
    if isinstance(idx, int):
        select, value = 'number', str(idx)
    else:
        select, value = 'symbol', idx
    return next(row for row in species_data if row[select] == value)[name]


null = None  # make the lines below a valid JSON
species_data = [
    {"number": 1, "symbol": "H", "name": "hydrogen", "vdw_radius": 1.2, "covalent_radius": 0.38, "mass": 1.0079},
    {"number": 2, "symbol": "He", "name": "helium", "vdw_radius": 1.4, "covalent_radius": 0.32, "mass": 4.0026},
    {"number": 3, "symbol": "Li", "name": "lithium", "vdw_radius": 1.82, "covalent_radius": 1.34, "mass": 6.941},
    {"number": 4, "symbol": "Be", "name": "beryllium", "vdw_radius": 1.53, "covalent_radius": 0.9, "mass": 9.0122},
    {"number": 5, "symbol": "B", "name": "boron", "vdw_radius": 1.92, "covalent_radius": 0.82, "mass": 10.811},
    {"number": 6, "symbol": "C", "name": "carbon", "vdw_radius": 1.7, "covalent_radius": 0.77, "mass": 12.0107},
    {"number": 7, "symbol": "N", "name": "nitrogen", "vdw_radius": 1.55, "covalent_radius": 0.75, "mass": 14.0067},
    {"number": 8, "symbol": "O", "name": "oxygen", "vdw_radius": 1.52, "covalent_radius": 0.73, "mass": 15.9994},
    {"number": 9, "symbol": "F", "name": "fluorine", "vdw_radius": 1.47, "covalent_radius": 0.71, "mass": 18.9984},
    {"number": 10, "symbol": "Ne", "name": "neon", "vdw_radius": 1.54, "covalent_radius": 0.69, "mass": 20.1797},
    {"number": 11, "symbol": "Na", "name": "sodium", "vdw_radius": 2.27, "covalent_radius": 1.54, "mass": 22.9897},
    {"number": 12, "symbol": "Mg", "name": "magnesium", "vdw_radius": 1.73, "covalent_radius": 1.3, "mass": 24.305},
    {"number": 13, "symbol": "Al", "name": "aluminium", "vdw_radius": 1.84, "covalent_radius": 1.18, "mass": 26.9815},
    {"number": 14, "symbol": "Si", "name": "silicon", "vdw_radius": 2.1, "covalent_radius": 1.11, "mass": 28.0855},
    {"number": 15, "symbol": "P", "name": "phosphorus", "vdw_radius": 1.8, "covalent_radius": 1.06, "mass": 30.9738},
    {"number": 16, "symbol": "S", "name": "sulfur", "vdw_radius": 1.8, "covalent_radius": 1.02, "mass": 32.065},
    {"number": 17, "symbol": "Cl", "name": "chlorine", "vdw_radius": 1.75, "covalent_radius": 0.99, "mass": 35.453},
    {"number": 18, "symbol": "Ar", "name": "argon", "vdw_radius": 1.88, "covalent_radius": 0.97, "mass": 39.948},
    {"number": 19, "symbol": "K", "name": "potassium", "vdw_radius": 2.75, "covalent_radius": 1.96, "mass": 39.0983},
    {"number": 20, "symbol": "Ca", "name": "calcium", "vdw_radius": 2.31, "covalent_radius": 1.74, "mass": 40.078},
    {"number": 21, "symbol": "Sc", "name": "scandium", "vdw_radius": 2.11, "covalent_radius": 1.44, "mass": 44.9559},
    {"number": 22, "symbol": "Ti", "name": "titanium", "vdw_radius": null, "covalent_radius": 1.36, "mass": 47.867},
    {"number": 23, "symbol": "V", "name": "vanadium", "vdw_radius": null, "covalent_radius": 1.25, "mass": 50.9415},
    {"number": 24, "symbol": "Cr", "name": "chromium", "vdw_radius": null, "covalent_radius": 1.27, "mass": 51.9961},
    {"number": 25, "symbol": "Mn", "name": "manganese", "vdw_radius": null, "covalent_radius": 1.39, "mass": 54.938},
    {"number": 26, "symbol": "Fe", "name": "iron", "vdw_radius": null, "covalent_radius": 1.25, "mass": 55.845},
    {"number": 27, "symbol": "Co", "name": "cobalt", "vdw_radius": null, "covalent_radius": 1.26, "mass": 58.9332},
    {"number": 28, "symbol": "Ni", "name": "nickel", "vdw_radius": 1.63, "covalent_radius": 1.21, "mass": 58.6934},
    {"number": 29, "symbol": "Cu", "name": "copper", "vdw_radius": 1.4, "covalent_radius": 1.38, "mass": 63.546},
    {"number": 30, "symbol": "Zn", "name": "zinc", "vdw_radius": 1.39, "covalent_radius": 1.31, "mass": 65.39},
    {"number": 31, "symbol": "Ga", "name": "gallium", "vdw_radius": 1.87, "covalent_radius": 1.26, "mass": 69.723},
    {"number": 32, "symbol": "Ge", "name": "germanium", "vdw_radius": 2.11, "covalent_radius": 1.22, "mass": 72.64},
    {"number": 33, "symbol": "As", "name": "arsenic", "vdw_radius": 1.85, "covalent_radius": 1.19, "mass": 74.9216},
    {"number": 34, "symbol": "Se", "name": "selenium", "vdw_radius": 1.9, "covalent_radius": 1.16, "mass": 78.96},
    {"number": 35, "symbol": "Br", "name": "bromine", "vdw_radius": 1.85, "covalent_radius": 1.14, "mass": 79.904},
    {"number": 36, "symbol": "Kr", "name": "krypton", "vdw_radius": 2.02, "covalent_radius": 1.1, "mass": 83.8},
    {"number": 37, "symbol": "Rb", "name": "rubidium", "vdw_radius": 3.03, "covalent_radius": 2.11, "mass": 85.4678},
    {"number": 38, "symbol": "Sr", "name": "strontium", "vdw_radius": 2.49, "covalent_radius": 1.92, "mass": 87.62},
    {"number": 39, "symbol": "Y", "name": "yttrium", "vdw_radius": null, "covalent_radius": 1.62, "mass": 88.9059},
    {"number": 40, "symbol": "Zr", "name": "zirconium", "vdw_radius": null, "covalent_radius": 1.48, "mass": 91.224},
    {"number": 41, "symbol": "Nb", "name": "niobium", "vdw_radius": null, "covalent_radius": 1.37, "mass": 92.9064},
    {"number": 42, "symbol": "Mo", "name": "molybdenum", "vdw_radius": null, "covalent_radius": 1.45, "mass": 95.94},
    {"number": 43, "symbol": "Tc", "name": "technetium", "vdw_radius": null, "covalent_radius": 1.56, "mass": 98},
    {"number": 44, "symbol": "Ru", "name": "ruthenium", "vdw_radius": null, "covalent_radius": 1.26, "mass": 101.07},
    {"number": 45, "symbol": "Rh", "name": "rhodium", "vdw_radius": null, "covalent_radius": 1.35, "mass": 102.9055},
    {"number": 46, "symbol": "Pd", "name": "palladium", "vdw_radius": 1.63, "covalent_radius": 1.31, "mass": 106.42},
    {"number": 47, "symbol": "Ag", "name": "silver", "vdw_radius": 1.72, "covalent_radius": 1.53, "mass": 107.8682},
    {"number": 48, "symbol": "Cd", "name": "cadmium", "vdw_radius": 1.58, "covalent_radius": 1.48, "mass": 112.411},
    {"number": 49, "symbol": "In", "name": "indium", "vdw_radius": 1.93, "covalent_radius": 1.44, "mass": 114.818},
    {"number": 50, "symbol": "Sn", "name": "tin", "vdw_radius": 2.17, "covalent_radius": 1.41, "mass": 118.71},
    {"number": 51, "symbol": "Sb", "name": "antimony", "vdw_radius": 2.06, "covalent_radius": 1.38, "mass": 121.76},
    {"number": 52, "symbol": "Te", "name": "tellurium", "vdw_radius": 2.06, "covalent_radius": 1.35, "mass": 127.6},
    {"number": 53, "symbol": "I", "name": "iodine", "vdw_radius": 1.98, "covalent_radius": 1.33, "mass": 126.9045},
    {"number": 54, "symbol": "Xe", "name": "xenon", "vdw_radius": 2.16, "covalent_radius": 1.3, "mass": 131.293},
    {"number": 55, "symbol": "Cs", "name": "caesium", "vdw_radius": 3.43, "covalent_radius": 2.25, "mass": 132.9055},
    {"number": 56, "symbol": "Ba", "name": "barium", "vdw_radius": 2.68, "covalent_radius": 1.98, "mass": 137.327},
    {"number": 57, "symbol": "La", "name": "lanthanum", "vdw_radius": null, "covalent_radius": 1.69, "mass": 138.9055},
    {"number": 58, "symbol": "Ce", "name": "cerium", "vdw_radius": null, "covalent_radius": null, "mass": 140.116},
    {"number": 59, "symbol": "Pr", "name": "praseodymium", "vdw_radius": null, "covalent_radius": null, "mass": 140.9077},
    {"number": 60, "symbol": "Nd", "name": "neodymium", "vdw_radius": null, "covalent_radius": null, "mass": 144.24},
    {"number": 61, "symbol": "Pm", "name": "promethium", "vdw_radius": null, "covalent_radius": null, "mass": 145},
    {"number": 62, "symbol": "Sm", "name": "samarium", "vdw_radius": null, "covalent_radius": null, "mass": 150.36},
    {"number": 63, "symbol": "Eu", "name": "europium", "vdw_radius": null, "covalent_radius": null, "mass": 151.964},
    {"number": 64, "symbol": "Gd", "name": "gadolinium", "vdw_radius": null, "covalent_radius": null, "mass": 157.25},
    {"number": 65, "symbol": "Tb", "name": "terbium", "vdw_radius": null, "covalent_radius": null, "mass": 158.9253},
    {"number": 66, "symbol": "Dy", "name": "dysprosium", "vdw_radius": null, "covalent_radius": null, "mass": 162.5},
    {"number": 67, "symbol": "Ho", "name": "holmium", "vdw_radius": null, "covalent_radius": null, "mass": 164.9303},
    {"number": 68, "symbol": "Er", "name": "erbium", "vdw_radius": null, "covalent_radius": null, "mass": 167.259},
    {"number": 69, "symbol": "Tm", "name": "thulium", "vdw_radius": null, "covalent_radius": null, "mass": 168.9342},
    {"number": 70, "symbol": "Yb", "name": "ytterbium", "vdw_radius": null, "covalent_radius": null, "mass": 173.04},
    {"number": 71, "symbol": "Lu", "name": "lutetium", "vdw_radius": null, "covalent_radius": 1.6, "mass": 174.967},
    {"number": 72, "symbol": "Hf", "name": "hafnium", "vdw_radius": null, "covalent_radius": 1.5, "mass": 178.49},
    {"number": 73, "symbol": "Ta", "name": "tantalum", "vdw_radius": null, "covalent_radius": 1.38, "mass": 180.9479},
    {"number": 74, "symbol": "W", "name": "tungsten", "vdw_radius": null, "covalent_radius": 1.46, "mass": 183.84},
    {"number": 75, "symbol": "Re", "name": "rhenium", "vdw_radius": null, "covalent_radius": 1.59, "mass": 186.207},
    {"number": 76, "symbol": "Os", "name": "osmium", "vdw_radius": null, "covalent_radius": 1.28, "mass": 190.23},
    {"number": 77, "symbol": "Ir", "name": "iridium", "vdw_radius": null, "covalent_radius": 1.37, "mass": 192.217},
    {"number": 78, "symbol": "Pt", "name": "platinum", "vdw_radius": 1.75, "covalent_radius": 1.28, "mass": 195.078},
    {"number": 79, "symbol": "Au", "name": "gold", "vdw_radius": 1.66, "covalent_radius": 1.44, "mass": 196.9665},
    {"number": 80, "symbol": "Hg", "name": "mercury", "vdw_radius": 1.55, "covalent_radius": 1.49, "mass": 200.59},
    {"number": 81, "symbol": "Tl", "name": "thallium", "vdw_radius": 1.96, "covalent_radius": 1.48, "mass": 204.3833},
    {"number": 82, "symbol": "Pb", "name": "lead", "vdw_radius": 2.02, "covalent_radius": 1.47, "mass": 207.2},
    {"number": 83, "symbol": "Bi", "name": "bismuth", "vdw_radius": 2.07, "covalent_radius": 1.46, "mass": 208.9804},
    {"number": 84, "symbol": "Po", "name": "polonium", "vdw_radius": 1.97, "covalent_radius": null, "mass": 209},
    {"number": 85, "symbol": "At", "name": "astatine", "vdw_radius": 2.02, "covalent_radius": null, "mass": 210},
    {"number": 86, "symbol": "Rn", "name": "radon", "vdw_radius": 2.2, "covalent_radius": 1.45, "mass": 222},
    {"number": 87, "symbol": "Fr", "name": "francium", "vdw_radius": 3.48, "covalent_radius": null, "mass": 223},
    {"number": 88, "symbol": "Ra", "name": "radium", "vdw_radius": 2.83, "covalent_radius": null, "mass": 226},
    {"number": 89, "symbol": "Ac", "name": "actinium", "vdw_radius": null, "covalent_radius": null, "mass": 227},
    {"number": 90, "symbol": "Th", "name": "thorium", "vdw_radius": null, "covalent_radius": null, "mass": 232.0381},
    {"number": 91, "symbol": "Pa", "name": "protactinium", "vdw_radius": null, "covalent_radius": null, "mass": 231.0359},
    {"number": 92, "symbol": "U", "name": "uranium", "vdw_radius": 1.86, "covalent_radius": null, "mass": 238.0289}
]
