import os, sys

tracy_dir = os.getenv('TRACY_LIB')
sys.path.append(tracy_dir+'/tracy/lib')

import libtracy as scsi


def printf(format, *args): sys.stdout.write(format % args); sys.stdout.flush()

def prt_mat(mat):
    for x in mat:
        for i, y in enumerate(x):
            printf('%14.6e', y)
        print('')

def prt_elems(scsi, lat):
    printf('\n')
    for k in range(0, lat.conf.Cell_nLoc+1):
        printf('%4d ', k)
        if not True:
            lat.elems[k].prt_elem()
            printf('\n')
        elif not True:
            lat.elems[k].print()
        else:
            printf('%2d\n', lat.elems[k].Pkind)
            if lat.elems[k].Pkind == scsi.PartsKind.drift:
                printf('%-8s %5.3f\n', lat.elems[k].Name, lat.elems[k].PL)
            elif lat.elems[k].Pkind == scsi.PartsKind.Mpole:
                printf('%-8s %5.3f %5.3f\n',
                       lat.elems[k].Name, lat.elems[k].PL,
                       lat.elems[k].PBpar[scsi.HOMmax+2])

def set_config():
    lat = scsi.LatticeType()

    lat.conf.trace        = False
    lat.conf.reverse_elem = True
    lat.conf.mat_meth     = not False

    lat.conf.H_exact    = False; lat.conf.quad_fringe    = False
    lat.conf.Cavity_on  = False; lat.conf.radiation      = False
    lat.conf.emittance  = False; lat.conf.IBS            = False
    lat.conf.pathlength = False; lat.conf.bpm            = 0
    lat.conf.Cart_Bend  = False; lat.conf.dip_edge_fudge = True

    return lat

def get_lat(file_name):
    lat = set_config()
    lat.Lat_Read_py(file_name);
    lat.Lat_Init();
    lat.ChamberOff();
    lat.conf.CODimax = 10;

    return lat
   
def main(file_name):
    lat = get_lat(file_name)

    lat.Ring_GetTwiss(True, 0e-3); lat.print();
    print('\nnu: ', lat.conf.TotalTune)

    print('\nLinear Poincare map:')
    prt_mat(lat.conf.OneTurnMat)

    print('\nA:')
    prt_mat(lat.conf.Ascr)

    lastpos = 0;
    lat.getcod(1e-3, lastpos);
    print('\nCod:', lat.conf.CODvect)

    print('\nPartsKind.Mpole:',
          scsi.PartsKind.Mpole.name, scsi.PartsKind.Mpole.value)

    if False:
        lat.prt_fams()
        lat.prt_elems()

    printf('\nHOMmax = %d\n', scsi.HOMmax)

    if False:
        prt_elems(scsi, lat)
    
    lat.prt_lat2('"linlat1.out', True);
    lat.prt_lat4('linlat.out', True, 10);
#    lat.prtmfile('flat_file.dat');


lat_dir = os.getenv('LAT')
if not True:
    file_name = lat_dir+'/BESSY-III/NoTG-TGRB-B60-6bend-6sx_JB_tracy'
else:
    file_name = lat_dir+'/DIAMOND-II/vmx_symopt_Qx205_Qy364'

main(file_name)
