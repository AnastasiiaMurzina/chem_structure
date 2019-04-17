import os
import sys
import subprocess

alias = '/opt/mopac/run_script.sh'
shut_up = os.devnull

def write_mop_html(fname_xyz, fmop = ''):
    fmop = fname_xyz.split('.')[0] + '.mop' if fmop == '' else fmop
    lines = []
    with open(fname_xyz, 'r') as xyz:
        n = int(xyz.readline())
        xyz.readline()
        for _ in range(n): lines.append(xyz.readline())
    with open(fmop, 'w') as wf:
        wf.write('HTML\n{0}\n\n'.format(fname_xyz))
        for i in lines:
            wf.write("{}\t{}\t1\t{}\t1\t{}\t1\n".format(*i.split()))
    return fmop


def write_locate_mop(arc_product, arc_reactant, flocate):
    open(flocate, 'w').write("""geo_dat="{}" geo_ref="{}" opt geo-ok locate-ts(C:10)
 Locate and refine the transition state for the reaction""".format(arc_reactant, arc_product))

def write_saddle_mop(mop_saddle, dat, ref):
    open(mop_saddle, 'w').write("""geo_dat="{} 10p0 first.mop" geo_ref="{} 10p0 second.mop" saddle geo-ok
 Locate and refine the transition state for the reaction""".format(dat, ref))


if __name__ == '__main__':
    freactant, fproduct = sys.argv[1], sys.argv[2] # reactant.xyz, product.xyz
    cur_dir = os.path.dirname(freactant) #probably choose tmp???

    reactant = str(os.path.basename(freactant).split('.')[0])
    fm1 = write_mop_html(freactant, fmop=str(reactant)+'.mop')
    # subprocess.call([alias, fm1, shut_up]) # html of reactant

    product = str(os.path.basename(fproduct).split('.')[0])
    fm2 = write_mop_html(fproduct, fmop=str(product)+'.mop')
    # subprocess.call([alias, fm2, shut_up]) # html of product

    arc_prod = os.path.join(cur_dir, product+'.arc')
    arc_reac = os.path.join(cur_dir, reactant+'.arc')
    loc_mop = os.path.join(cur_dir, 'locate_{}_{}.mop'.format(reactant, product))
    if os.path.isfile(arc_prod) and os.path.isfile(arc_reac):
        write_locate_mop(arc_prod, arc_reac, loc_mop)
        # subprocess.call([alias, loc_mop]) # locate ts

        data_ref = '.'.join(loc_mop.split('.')[:-1:])
        saddle_mop = data_ref+'step2.mop'
        write_saddle_mop(saddle_mop, data_ref, data_ref)
        subprocess.call([alias, saddle_mop]) # saddle point

        # TODO ts
        # ts  geo_dat="C2H4+C4H6-cyclohexene Step 2.arc" ts html
        # Given a rough transition state, refine it using TS

        # oldgeo force

        #TODO irc

        #irc=1* html x-priority=0.2 geo_dat="C2H4+C4H6-cyclohexene_ts_then_force.arc"
        # Diels-Alder addition: C2H4+C4H6 => Cyclohexene

    else:
        print("Something wrong with mop reactant or product precalculations")
