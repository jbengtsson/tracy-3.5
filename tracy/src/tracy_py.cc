
namespace py = pybind11;

int no_tps = 1;

PYBIND11_MODULE(libtracy, scsi) {
    scsi.doc() = "Self-Consistent Symplectic Integrator (SCSI)";

    py::class_<ConfigType> (scsi, "ConfigType")
      .def_readwrite("CODimax", &ConfigType::CODimax)
      .def(py::init<>());

    py::class_<CellType>   (scsi, "CellType")
      .def(py::init<>());

    py::class_<ElemType>   (scsi, "ElemType")
      .def(py::init<>());
    
    scsi.def("printglob", &printglob);

    py::class_<LatticeType>(scsi, "LatticeType")
      .def_readwrite("conf", &LatticeType::conf)
      .def(py::init<>())
      .def("Lat_Read", &LatticeType::Lat_Read)
      .def("Lat_Init", &LatticeType::Lat_Init)
      .def("ChamberOff", &LatticeType::ChamberOff)
      .def("Ring_GetTwiss", &LatticeType::Ring_GetTwiss);
}
