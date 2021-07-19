
namespace py = pybind11;

int no_tps = 1;

PYBIND11_MODULE(libtracy, scsi) {
    scsi.doc() = "Self-Consistent Symplectic Integrator (SCSI)";

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
      .def_readwrite("CODimax",        &ConfigType::CODimax)
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
