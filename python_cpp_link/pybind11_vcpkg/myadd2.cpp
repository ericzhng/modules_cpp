
#include <pybind11/pybind11.h>
namespace py = pybind11;

int add(int a, int b) {
    return a + b;
}

PYBIND11_MODULE(myadd2, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("add", &add, "A function that adds two numbers",
        py::arg("i"), py::arg("j"));
}
