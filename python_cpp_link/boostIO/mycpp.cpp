#include <boost/python.hpp>
#include <iostream>

class Greeter {
public:
    void say_hello(const std::string& name) {
        std::cout << "Hello, " << name << " from C++!" << std::endl;
    }
};

// Expose the class to Python
BOOST_PYTHON_MODULE(mycpp) {
    using namespace boost::python;
    class_<Greeter>("Greeter")
        .def("say_hello", &Greeter::say_hello);
}
