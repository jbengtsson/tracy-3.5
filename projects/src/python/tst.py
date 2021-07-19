import os, sys

tracy_dir = os.getenv('TRACY_LIB')
sys.path.append(tracy_dir+'/tracy/lib')
import libtracy as scsi

lat = scsi.LatticeType()

lat.conf.trace        = False;
lat.conf.reverse_elem = not False;
lat.conf.mat_meth     = not False;

lat.conf.H_exact    = False; lat.conf.quad_fringe    = False;
lat.conf.Cavity_on  = False; lat.conf.radiation      = False;
lat.conf.emittance  = False; lat.conf.IBS            = False;
lat.conf.pathlength = False; lat.conf.bpm            = 0;
lat.conf.Cart_Bend  = False; lat.conf.dip_edge_fudge = True;

lat_dir = os.getenv('LAT')
lat.Lat_Read(lat_dir+'/BESSY-III/NoTG-TGRB-B60-6bend-6sx_JB_tracy');
lat.Lat_Init();

lat.ChamberOff();

lat.conf.CODimax = 10;

lat.Ring_GetTwiss(True, 0e-3); lat.print();
