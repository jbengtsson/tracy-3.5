#include <pybind11/stl.h>
#include <pybind11/pybind11.h>
// #include <pybind11/operators.h>
#include <vector>
#include <tuple>

#include "tracy_lib.h"

namespace py = pybind11;


// Polymorphic number template class wrapper.

template<typename T>
void declare_field(py::module &scsi, const std::string &typestr) {
  using Class = ss_vect<T>;
  std::string pyclass_name = std::string("ss_vect") + typestr;
  py::class_<Class>(scsi, pyclass_name.c_str(), py::buffer_protocol(),
		    py::dynamic_attr())
    .def(py::init<>())
    // Python print redirect.
    .def("__repr__", [](const ss_vect<T> &a) {
		       std::stringstream stream;
		       stream << a;
		       return stream.str();
		     })
    .def("print",      &Class::print)
    .def("__getitem__", [](const ss_vect<T> &a, const int k) {
			  return a[k];
			})
    .def("__setitem__", [](ss_vect<T> &a, const int k, const T &x) {
			  a[k] = x;
			})
    .def("zero",       &Class::zero)
    .def("identity",   &Class::identity);
}


/*
 * get_eps_x: provide the memory and return it to python
 *            TODO: check for memory leaks!
 */
class LatticeTypePy : public LatticeType {

public:
    std::tuple<double, double, double, std::vector<double>,std::vector<double>,std::vector<double>>
    get_eps_x_py(const bool prt){
		double eps_x, sigma_delta, U_0;
		// TODO: allocate appropriate space
		std::vector<double> J(3, 0.0), tau(3, 0.0), I(3, 0.0);

		get_eps_x(eps_x, sigma_delta, U_0, J, tau, I, prt);

		std::tuple<double, double, double, std::vector<double>,std::vector<double>,std::vector<double>>
		  tup( make_tuple(eps_x, sigma_delta, U_0, J, tau, I));
		return tup;
	}
};


PYBIND11_MODULE(libtracy, scsi) {
    scsi.doc() = "Self-Consistent Symplectic Integrator (SCSI)";

    // Polymorphic number class.

    // Enums.

    py::enum_<spatial_ind>(scsi, "spatial_ind")
      .value("X_", spatial_ind::X_)
      .value("Y_", spatial_ind::Y_)
      .value("Z_", spatial_ind::Z_);

    py::enum_<phase_space_ind>(scsi, "phase_space_ind")
      .value("x_",     phase_space_ind::x_)
      .value("px_",    phase_space_ind::px_)
      .value("y_",     phase_space_ind::y_)
      .value("py_",    phase_space_ind::py_)
      .value("delta_", phase_space_ind::delta_)
      .value("ct_",    phase_space_ind::ct_);

    // Classes.

    py::class_<tps>(scsi, "tps")
      .def(py::init<>())
      .def(py::init<const double>())
      .def(py::init<const double, const int>())
      // Python print redirect.
      .def("print", &tps::print);

    declare_field<double>(scsi, "_double");
    declare_field<tps>(scsi, "_tps");

    // Beamline class.

    // Constants.

    scsi.attr("HOMmax") = HOMmax;
    scsi.attr("c0")     = c0;
    scsi.attr("q_e")    = q_e;
    scsi.attr("m_e")    = m_e;
    scsi.attr("mu_0")   = mu_0;
    scsi.attr("eps_0")  = eps_0;
    scsi.attr("r_e")    = r_e;
    scsi.attr("h_bar")  = h_bar;

    // Enums.

    py::enum_<PartsKind>(scsi, "PartsKind")
      .value("drift",      PartsKind::drift)
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

    py::enum_<MpoleKind>(scsi, "MpoleKind")
      .value("All",   MpoleKind::All)
      .value("Dip",   MpoleKind::Dip)
      .value("Quad",  MpoleKind::Quad)
      .value("Sext",  MpoleKind::Sext)
      .value("Oct",   MpoleKind::Oct)
      .value("Dec",   MpoleKind::Dec)
      .value("Dodec", MpoleKind::Dodec);

    py::enum_<PlaneKind>(scsi, "Horizontal")
      .value("Horizontal", PlaneKind::Horizontal)
      .value("Vertical",   PlaneKind::Vertical);

    // Classes.

    py::class_<ConfigType>(scsi, "ConfigType")
      // bool
      .def_readwrite("trace",          &ConfigType::trace)
      .def_readwrite("reverse_elem",   &ConfigType::reverse_elem)
      .def_readwrite("stable",         &ConfigType::stable)
      .def_readwrite("ErrFlag",        &ConfigType::ErrFlag)
      .def_readwrite("Cavity_on",      &ConfigType::Cavity_on)
      .def_readwrite("radiation",      &ConfigType::radiation)
      .def_readwrite("emittance",      &ConfigType::emittance)
      .def_readwrite("quad_fringe",    &ConfigType::quad_fringe)
      .def_readwrite("H_exact",        &ConfigType::H_exact)
      .def_readwrite("Cart_Bend",      &ConfigType::Cart_Bend)
      .def_readwrite("dip_edge_fudge", &ConfigType::dip_edge_fudge)
      .def_readwrite("pathlength",     &ConfigType::pathlength)
      .def_readwrite("Aperture_on",    &ConfigType::Aperture_on)
      .def_readwrite("EPU",            &ConfigType::EPU)
      .def_readwrite("mat_meth",       &ConfigType::mat_meth)
      .def_readwrite("IBS",            &ConfigType::IBS)
      .def_readwrite("tuneflag",       &ConfigType::tuneflag)
      .def_readwrite("chromflag",      &ConfigType::chromflag)
      .def_readwrite("codflag",        &ConfigType::codflag)
      .def_readwrite("mapflag",        &ConfigType::mapflag)
      .def_readwrite("passflag",       &ConfigType::passflag)
      .def_readwrite("overflag",       &ConfigType::overflag)
      .def_readwrite("chambre",        &ConfigType::chambre)
      // long int
      .def_readwrite("Cell_nLoc",      &ConfigType::Cell_nLoc)
      .def_readwrite("Elem_nFam",      &ConfigType::Elem_nFam)
      .def_readwrite("CODimax",        &ConfigType::CODimax)
      // int
      .def_readwrite("bpm",            &ConfigType::bpm)
      .def_readwrite("hcorr",          &ConfigType::hcorr)
      .def_readwrite("vcorr",          &ConfigType::vcorr)
      .def_readwrite("qt",             &ConfigType::qt)
      .def_readwrite("gs",             &ConfigType::gs)
      .def_readwrite("ge",             &ConfigType::ge)
      .def_readwrite("RingType",       &ConfigType::RingType)
      .def_readwrite("lossplane",      &ConfigType::lossplane)
      // double
      .def_readwrite("dPcommon",       &ConfigType::dPcommon)
      .def_readwrite("dPparticle",     &ConfigType::dPparticle)
      .def_readwrite("delta_RF",       &ConfigType::delta_RF)
      .def_readwrite("Omega",          &ConfigType::Omega)
      .def_readwrite("U0",             &ConfigType::U0)
      .def_readwrite("Alphac",         &ConfigType::Alphac)
      .def_readwrite("Energy",         &ConfigType::Energy)
      .def_readwrite("dE",             &ConfigType::dE)
      .def_readwrite("CODeps",         &ConfigType::CODeps)
      .def_readwrite("Qb",             &ConfigType::Qb)
      .def_readwrite("alpha_z",        &ConfigType::alpha_z)
      .def_readwrite("beta_z",         &ConfigType::beta_z)
      .def_readwrite("beta0",          &ConfigType::beta0)
      .def_readwrite("gamma0",         &ConfigType::gamma0)
      // std::vector<double>:
      .def_readwrite("TotalTune",      &ConfigType::TotalTune)
      .def_readwrite("Chrom",          &ConfigType::Chrom)
      .def_readwrite("alpha_rad",      &ConfigType::alpha_rad)
      .def_readwrite("D_rad",          &ConfigType::D_rad)
      .def_readwrite("J",              &ConfigType::J)
      .def_readwrite("tau",            &ConfigType::tau)
      .def_readwrite("D_IBS",          &ConfigType::D_IBS)
      .def_readwrite("eps",            &ConfigType::eps)
      .def_readwrite("epsp",           &ConfigType::epsp)
      .def_readwrite("CODvect",        &ConfigType::CODvect)
      .def_readwrite("wr",             &ConfigType::wr)
      .def_readwrite("wi",             &ConfigType::wi)
      // std::vector< std::vector<double> >
      .def_readwrite("OneTurnMat",     &ConfigType::OneTurnMat)
      .def_readwrite("Ascr",           &ConfigType::Ascr)
      .def_readwrite("Ascrinv",        &ConfigType::Ascrinv)
      .def_readwrite("Vr",             &ConfigType::Vr)
      .def_readwrite("Vi",             &ConfigType::Vi)
      .def(py::init<>());

    py::class_<CellType>(scsi, "CellType")
      // int
      .def_readwrite("Fnum",       &CellType::Fnum)
      .def_readwrite("Knum",       &CellType::Knum)
      // double
      .def_readwrite("S",          &CellType::S)
      .def_readwrite("curly_dH_x", &CellType::curly_dH_x)
      // std::vector<double>:
      .def_readwrite("dI",         &CellType::dI)
      .def_readwrite("dS",         &CellType::dS)
      .def_readwrite("dT",         &CellType::dT)
      .def_readwrite("Eta",        &CellType::Eta)
      .def_readwrite("Etap",       &CellType::Etap)
      .def_readwrite("Alpha",      &CellType::Alpha)
      .def_readwrite("Beta",       &CellType::Beta)
      .def_readwrite("Nu",         &CellType::Nu)
      .def_readwrite("BeamPos",    &CellType::BeamPos)
      // std::vector< std::vector<double> >
      .def_readwrite("maxampl",    &CellType::maxampl)
      .def_readwrite("A",          &CellType::A)
      .def_readwrite("sigma",      &CellType::sigma)
      // CellType
      .def(py::init<>());

    py::class_<ElemType>(scsi, "ElemType")
      // std::string
      .def_readwrite("Name",    &ElemType::Name)
      // bool
      .def_readwrite("Reverse", &ElemType::Reverse)
      // double
      .def_readwrite("PL",      &ElemType::PL)
      // PartsKind
      .def_readwrite("Pkind",   &ElemType::Pkind)
      // Virtual class: no constructor.
      // .def(py::init<>())
      .def("prt_elem", &ElemType::prt_elem)
      .def("print",    &ElemType::print);

    py::class_<ElemFamType>(scsi, "ElemFamType")
      // ElemType
      // int
      .def_readwrite("nKid",    &ElemFamType::nKid)
      .def_readwrite("NoDBN",   &ElemFamType::NoDBN)
      // std::vector<int>
      .def_readwrite("KidList", &ElemFamType::KidList)
      // std::vector<string>
      .def_readwrite("DBNlist", &ElemFamType::DBNlist)
      .def(py::init<>());

    py::class_<LatticeTypePy>(scsi, "LatticeType")
      // std::vector<ElemFamType>
      .def_readwrite("elemf", &LatticeTypePy::elemf)
      // std::vector<ElemType*>
      .def_readwrite("elems", &LatticeTypePy::elems)
      // ConfigType
      .def_readwrite("conf",  &LatticeTypePy::conf)
      .def(py::init<>())
      .def("Lat_Init",      &LatticeTypePy::Lat_Init)
      .def("prt_fams",      &LatticeTypePy::prt_fams)
      .def("prt_elems",     &LatticeTypePy::prt_elems)
      .def("GetnKid",       &LatticeTypePy::GetnKid)
      .def("ElemIndex",     &LatticeTypePy::ElemIndex)
      .def("Elem_GetPos",   &LatticeTypePy::Elem_GetPos)
      .def("Lat_Read",      &LatticeTypePy::Lat_Read)
      .def("prtmfile",      &LatticeTypePy::prtmfile)
      .def("rdmfile",       &LatticeTypePy::rdmfile)
      .def("prt_lat1",      &LatticeTypePy::prt_lat1)
      .def("prt_lat2",      &LatticeTypePy::prt_lat2)
      .def("prt_lat3",      &LatticeTypePy::prt_lat3)
      .def("prt_lat4",      &LatticeTypePy::prt_lat4)
      .def("getcod",        &LatticeTypePy::getcod)
      .def("Ring_GetTwiss", &LatticeTypePy::Ring_GetTwiss)
      .def("ChamberOff",    &LatticeTypePy::ChamberOff)
      .def("print",         &LatticeTypePy::print)
      .def("get_eps_x",     &LatticeTypePy::get_eps_x_py)
      .def("GetEmittance",  &LatticeTypePy::GetEmittance)
      .def("Cell_Pass1",    &LatticeTypePy::Cell_Pass1)
      .def("Cell_Pass2",    &LatticeTypePy::Cell_Pass2);

    py::class_<DriftType, ElemType>(scsi, "DriftType")
      .def(py::init<>());

    py::class_<MpoleType, ElemType>(scsi, "MpoleType")
      // int
      .def_readwrite("Pmethod",  &MpoleType::Pmethod)
      .def_readwrite("PN",       &MpoleType::PN)
      .def_readwrite("Porder",   &MpoleType::Porder)
      .def_readwrite("n_design", &MpoleType::n_design)
      // double
      .def_readwrite("PdTpar",   &MpoleType::PdTpar)
      .def_readwrite("PdTsys",   &MpoleType::PdTsys)
      .def_readwrite("PdTrms",   &MpoleType::PdTrms)
      .def_readwrite("PdTrnd",   &MpoleType::PdTrnd)
      .def_readwrite("PTx1",     &MpoleType::PTx1)
      .def_readwrite("PTx2",     &MpoleType::PTx2)
      .def_readwrite("Pgap",     &MpoleType::Pgap)
      .def_readwrite("Pirho",    &MpoleType::Pirho)
      .def_readwrite("Pc0",      &MpoleType::Pc0)
      .def_readwrite("Pc1",      &MpoleType::Pc1)
      .def_readwrite("Ps1",      &MpoleType::Ps1)
      // std::vector<double>
      .def_readwrite("PdSsys",   &MpoleType::PdSsys)
      .def_readwrite("PdSrms",   &MpoleType::PdSrms)
      .def_readwrite("PdSrnd",   &MpoleType::PdSrnd)
      // MpoleArray
      .def_readwrite("PBpar",    &MpoleType::PBpar)
      .def_readwrite("PBsys",    &MpoleType::PBsys)
      .def_readwrite("PBrms",    &MpoleType::PBrms)
      .def_readwrite("PBrnd",    &MpoleType::PBrnd)
      .def_readwrite("PB",       &MpoleType::PB)
      // pthicktype
      .def_readwrite("Pthick",   &MpoleType::Pthick)
      // std::vector< std::vector<double> >
      .def_readwrite("M_elem",   &MpoleType::M_elem)
      .def(py::init<>());

    py::class_<CavityType,     ElemType>(scsi, "CavityType")
      .def(py::init<>());

    py::class_<MarkerType,     ElemType>(scsi, "MarkerType")
      .def(py::init<>());

    py::class_<WigglerType,    ElemType>(scsi, "WigglerType")
      .def(py::init<>());

    py::class_<InsertionType,  ElemType>(scsi, "InsertionType")
      .def(py::init<>());

    py::class_<FieldMapType,   ElemType>(scsi, "FieldMapType")
      .def(py::init<>());

    py::class_<SpreaderType,   ElemType>(scsi, "SpreaderType")
      .def(py::init<>());

    py::class_<RecombinerType, ElemType>(scsi, "RecombinerType")
      .def(py::init<>());

    py::class_<SolenoidType,   ElemType>(scsi, "SolenoidType")
      .def(py::init<>());

    py::class_<MapType,        ElemType>(scsi, "MapType")
      .def(py::init<>());
}
