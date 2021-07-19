import os, sys

tracy_dir = os.getenv('TRACY_LIB')
sys.path.append(tracy_dir+'/tracy/lib')
import libtracy as scsi

lat = scsi.LatticeType()

lat_dir = os.getenv('LAT')
lat.Lat_Read(lat_dir+'/BESSY-III/NoTG-TGRB-B60-6bend-6sx_JB_tracy');

lat.Lat_Init();

lat.ChamberOff();

lat.conf.CODimax = 10;

lat.Ring_GetTwiss(True, 0e-3); lat.print();
