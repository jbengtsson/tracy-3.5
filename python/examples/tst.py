import numpy as np
import tracy.lib as scsi


# Constants.

x_ = 0
px_ = 1
y_ = 2
py_ = 3
delta_ = 4
ct_ = 5


def prt_mat(mat):
    """

    Todo:
        Check if matrix should not be transposed before printing
    """
    for row in mat:
        row_txt = ["{:14.6e}".format(val) for val in row]
        print(" ".join(row_txt))


def prt_elems(scsi, lat):
    print("\n")
    for k in range(0, lat.conf.Cell_nLoc + 1):
        print("{:4d} ".format(k))
        if not True:
            lat.elems[k].prt_elem()
            print("\n")
        elif not True:
            # lat.elems[k].print()
            pass
        else:
            print("{:2d}\n".format(lat.elems[k].Pkind))
            if lat.elems[k].Pkind == scsi.PartsKind.drift:
                print("{:-8s} {:5.3f}\n".format(lat.elems[k].Name, lat.elems[k].PL))
            elif lat.elems[k].Pkind == scsi.PartsKind.Mpole:
                print(
                    "{:-8s} {:5.3f} {:5.3f}\n".format(
                        lat.elems[k].Name,
                        lat.elems[k].PL,
                        lat.elems[k].PBpar[scsi.HOMmax + 2],
                    )
                )


def set_config():
    lat = scsi.LatticeType()

    lat.conf.trace = False
    lat.conf.reverse_elem = True
    lat.conf.mat_meth = not False

    lat.conf.H_exact = False
    lat.conf.quad_fringe = False
    lat.conf.Cavity_on = False
    lat.conf.radiation = False
    lat.conf.emittance = False
    lat.conf.IBS = False
    lat.conf.pathlength = False
    lat.conf.bpm = 0
    lat.conf.Cart_Bend = False
    lat.conf.dip_edge_fudge = True

    return lat


def get_lat(file_name):
    lat = set_config()
    lat.Lat_Read(file_name)
    lat.Lat_Init()
    lat.ChamberOff()
    lat.conf.CODimax = 10

    return lat


def main(file_name):
    lat = get_lat(file_name)

    lat.Ring_GetTwiss(True, 0e-3)
    lat.print("")
    print("\nnu: ", lat.conf.TotalTune)

    # Compute radiation effects.
    if lat.conf.mat_meth:
        eps_x, sigma_delta, U_0, J, tau, I = lat.get_eps_x(True)
    else:
        lat.GetEmittance(lat.ElemIndex("cav"), True)

    print("\nE [GeV]: {:3.1f}\n", lat.conf.Energy)

    print("\nLinear Poincare map:")
    prt_mat(lat.conf.OneTurnMat)

    print("\nA:")
    prt_mat(lat.conf.Ascr)

    # Compute off-momentum closed orbit.
    lastpos = 0
    lat.getcod(1e-3, lastpos)

    fmt = """
    Cod :            {cod}
    HOM max :        {hom_max:d}
    PartsKind.Mpole: {mname} = {mvalue}
    """
    print(
        fmt.format(
            cod=lat.conf.CODvect,
            hom_max=scsi.HOMmax,
            mname=scsi.PartsKind.Mpole.name,
            mvalue=scsi.PartsKind.Mpole.value,
        )
    )

    if False:
        lat.prt_fams()
        lat.prt_elems()

    if False:
        prt_elems(scsi, lat)

    # Locate 1st dipole & print transport matrix.
    for k in range(0, lat.conf.Cell_nLoc + 1):
        if (lat.elems[k].Pkind == scsi.PartsKind.Mpole) and (
            lat.elems[k].n_design == scsi.MpoleKind.Dip
        ):
            print("\n{} M_elem:\n".format(lat.elems[k].Name))
            prt_mat(lat.elems[k].M_elem)
            break

    lat.prt_lat2("linlat1.out", True)
    lat.prt_lat4("linlat.out", True, 10)
    lat.prtmfile("flat_file.dat")

    # Polymorphic number class.

    ps = scsi.tps()
    ps.print("\n")
    ps = scsi.tps(1e0)
    ps.print("\n")
    ps = scsi.tps(0e0, 2)
    ps.print("\n")

    ps = scsi.ss_vect_double()
    ps.zero()
    ps.print("\nps:-----\n")

    if False:
        Id = scsi.ss_vect_tps()
        print(Id)
        Id.identity()
        Id.print("\nId:\n")
        # print('\nId:\n', Id)

        Id[3].print("\n")
        Id[0] = scsi.tps(0e0, 3)
        Id[0].print("\n")
        Id[0] = Id[3]
        Id[0].print("\n")

    # Compute on-momentum closed orbit.
    lastpos = 0
    lat.getcod(0e0, lastpos)
    print("\nCod:", lat.conf.CODvect)

    if False:
        # Track one turn.
        ps = scsi.ss_vect_double()
        ps.zero()
        ps[x_] = 1e-3
        ps[y_] = -1e-3
        lat.Cell_Pass1(0, lat.conf.Cell_nLoc, ps, lastpos)
        ps.print("\n")

    if False:
        # Compute linear Poincare map.
        M = scsi.ss_vect_tps()
        M.identity()
        lat.Cell_Pass2(0, lat.conf.Cell_nLoc, M, lastpos)
        M.print("\n")


if __name__ == "__main__":
    import os
    import sys

    # tracy_dir = os.getenv("TRACY_LIB")
    # sys.path.append(tracy_dir + "/tracy/lib")
    # Lattices.
    lat_dir = os.getenv("LAT")
    lattices = {
        # "bessy-iii": lat_dir + "/BESSY-III/NoTG-TGRB-B60-6bend-6sx_JB_tracy",
        # "diamond": lat_dir + "/DIAMOND-II/vmx_symopt_Qx205_Qy364",
        "bessy-ii": "tests/lattices/b2_stduser_beamports_blm_tracy_corr",
    }
    main(lattices["bessy-ii"])
