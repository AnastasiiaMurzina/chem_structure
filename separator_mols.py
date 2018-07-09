import os, subprocess, shutil, tempfile
from mol2_worker import xyz_names_bonds

FNULL = open(os.devnull, 'w')

# def get_adjency(file_name):


if __name__ == '__main__':
    initial_xyz = 'mols_dir/Ethene.xyz'
    tmp_dir = tempfile.mkdtemp()
    initial_mol2 = os.path.join(tmp_dir, 'Ethene.mol2')
    subprocess.call(['babel', '-ixyz', initial_xyz, '-omol2', initial_mol2], stdout=FNULL)
    bonds, _  = xyz_names_bonds(initial_mol2)
    connected_sets = [{bonds[0][0], bonds[0][1]}]
    for i in sorted(bonds, key=lambda x: min(x[0], x[1]))[1:]:
        flag_new = True
        for j in connected_sets:
            if i[0] in j:
                j.add(i[1])
                flag_new = False
                break
            elif i[1] in j:
                flag_new = False
                j.add(i[0])
                break
        if flag_new:
            connected_sets.append({})
        pass
    # print(bonds)
    # for i in bonds:
    #     flag_new = True
    #     for j in set_of_dicts:
    #         if i[0] in j:
    #             flag_new = False
    #             break

    # mols = 'mols_dir/Ethene.xyz'
    shutil.rmtree(tmp_dir)