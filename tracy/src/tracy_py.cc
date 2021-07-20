#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace py = pybind11;

int no_tps = 1;

PYBIND11_MODULE(libtracy, scsi) {
    scsi.doc() = "Self-Consistent Symplectic Integrator (SCSI)";

    py::enum_<PartsKind>(scsi, "PartsKind")
      .value("Wigl",       PartsKind::Wigl)
      .value("Mpole",      PartsKind::Mpole)
      .value("Cavity",     PartsKind::Cavity)
      .value("marker",     PartsKind::marker)
      .value("undef",      PartsKind::undef)
      .value("Insertion",  PartsKind::Insertion)
      .value("FieldMap",   PartsKind::FieldMap)
      .value("Spreader",   PartsKind::Spreader)
      .value("Recombiner", PartsKind::Recombiner)
      .value("Solenoid",   PartsKind::Solenoid)
      .value("Map",        PartsKind::Map)
      .export_values();

    py::class_<ConfigType> (scsi, "ConfigType")
      .def_readwrite("trace",          &ConfigType::trace)
      .def_readwrite("reverse_elem",   &ConfigType::reverse_elem)
      .def_readwrite("mat_meth",       &ConfigType::mat_meth)
      .def_readwrite("H_exact",        &ConfigType::H_exact)
      .def_readwrite("quad_fringe",    &ConfigType::quad_fringe)
      .def_readwrite("Cavity_on",      &ConfigType::Cavity_on)
      .def_readwrite("radiation",      &ConfigType::radiation)
      .def_readwrite("emittance",      &ConfigType::emittance)
      .def_readwrite("IBS",            &ConfigType::IBS)
      .def_readwrite("pathlength",     &ConfigType::pathlength)
      .def_readwrite("bpm",            &ConfigType::bpm)
      .def_readwrite("Cart_Bend",      &ConfigType::Cart_Bend)
      .def_readwrite("dip_edge_fudge", &ConfigType::dip_edge_fudge)
      .def_readwrite("Cell_nLoc",      &ConfigType::Cell_nLoc)
      .def_readwrite("Elem_nFam",      &ConfigType::Elem_nFam)
      .def_readwrite("CODimax",        &ConfigType::CODimax)
      .def_readwrite("TotalTune",      &ConfigType::TotalTune)
      .def_readwrite("alpha_rad",      &ConfigType::alpha_rad)
      .def_readwrite("D_rad",          &ConfigType::D_rad)
      .def_readwrite("J",              &ConfigType::J)
      .def_readwrite("tau",            &ConfigType::tau)
      .def_readwrite("D_IBS",          &ConfigType::D_IBS)
      .def_readwrite("eps",            &ConfigType::eps)
      .def_readwrite("epsp",           &ConfigType::epsp)
      .def_readwrite("Chrom",          &ConfigType::Chrom)
      .def(py::init<>());

    py::class_<CellType>   (scsi, "CellType")
      .def_readwrite("Fnum",    &CellType::Fnum)
      .def_readwrite("Knum",    &CellType::Knum)
      .def_readwrite("S",       &CellType::S)
      .def_readwrite("Alpha",   &CellType::Alpha)
      .def_readwrite("Beta",    &CellType::Beta)
      .def_readwrite("Eta",     &CellType::Eta)
      .def_readwrite("Etap",    &CellType::Etap)
      .def_readwrite("BeamPos", &CellType::BeamPos)
      .def(py::init<>());

    py::class_<ElemType>   (scsi, "ElemType")
      .def(py::init<>());

    py::class_<ElemFamType>   (scsi, "ElemFamType")
      .def(py::init<>());

    py::class_<LatticeType>(scsi, "LatticeType")
      .def_readwrite("elemf", &LatticeType::elemf)
      .def_readwrite("elems", &LatticeType::elems)
      .def_readwrite("conf",  &LatticeType::conf)
      .def(py::init<>())
      .def("Lat_Read",      &LatticeType::Lat_Read)
      .def("Lat_Init",      &LatticeType::Lat_Init)
      .def("ChamberOff",    &LatticeType::ChamberOff)
      .def("Ring_GetTwiss", &LatticeType::Ring_GetTwiss)
      .def("print",         &LatticeType::print);
}
