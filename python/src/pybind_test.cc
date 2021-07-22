// -*- c-set-style: linux -*-
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/iostream.h>
#include <iostream>
#include <sstream>


class TestElement
{
	double length;
	std::string name;

public:
	TestElement(std::string name, double length): name(name){
		this->length = length;
	}
	const std::string &getName() const { return this->name; }

	void setName(std::string &name){
		this->name = name;
	}
	const double &getLength() const{
		return this->length;
	}
	void setLength(double length){
		this->length = length;
	}


	friend std::ostream& operator<<(std::ostream& os, const TestElement& elem);

};

std::ostream& operator<<(std::ostream& os, const TestElement& elem)
{
	os << "TestElement(" << elem.name << ", " << elem.length << ")";
	return os;
}

namespace py = pybind11;


PYBIND11_MODULE(pybind_test, m) {
	m.doc() = "pybind11 example ";

	py::class_<TestElement>(m, "TestElement")
		.def(py::init<const std::string &, double>())
		.def_property("name", &TestElement::getName, &TestElement::setName)
		.def_property("length", &TestElement::getLength, &TestElement::setLength)
		.def("__repr__",
		     [](const TestElement &elem){
			     std::stringstream repr;
			     repr << "TestElement(" << elem.getName()
				  <<  ", " << elem.getLength()
				  << ")";
			     return repr.str();
		     }
		 );

	// from pybind11 documentation
	py::add_ostream_redirect(m, "ostream_redirect");
}
