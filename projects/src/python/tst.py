import os, sys

tracy_dir = os.getenv('TRACY_LIB')
sys.path.append(tracy_dir+'/tracy/lib')
import libtracy as scsi

HOMmax = 21

def printf(format, *args): sys.stdout.write(format % args); sys.stdout.flush()

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
                       lat.elems[k].PBpar[HOMmax+2])

def set_config():
    lat = scsi.LatticeType()

    lat.conf.trace        = False;
    lat.conf.reverse_elem = not False;
    lat.conf.mat_meth     = not False;

    lat.conf.H_exact    = False; lat.conf.quad_fringe    = False;
    lat.conf.Cavity_on  = False; lat.conf.radiation      = False;
    lat.conf.emittance  = False; lat.conf.IBS            = False;
    lat.conf.pathlength = False; lat.conf.bpm            = 0;
    lat.conf.Cart_Bend  = False; lat.conf.dip_edge_fudge = True;

    return lat

def get_lat(file_name):
    lat = set_config()
    lat.Lat_Read(file_name);
    lat.Lat_Init();
    lat.ChamberOff();
    lat.conf.CODimax = 10;

    return lat
   
def main(file_name):
    lat = get_lat(file_name)
    lat.Ring_GetTwiss(True, 0e-3); lat.print();
    print('\nnu:', lat.conf.TotalTune)
    print('\nPartsKind.Mpole:',
          scsi.PartsKind.Mpole.name, scsi.PartsKind.Mpole.value)
    if False:
        lat.prt_fams()
        lat.prt_elems()
    prt_elems(scsi, lat)

    return lat


lat_dir = os.getenv('LAT')
file_name = lat_dir+'/BESSY-III/NoTG-TGRB-B60-6bend-6sx_JB_tracy'
lat = main(file_name)
