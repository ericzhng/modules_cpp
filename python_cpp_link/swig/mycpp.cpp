#include <iostream>

class Greeter {
public:
    void say_hello(const std::string& name) {
        std::cout << "Hello, " << name << " from C++!" << std::endl;
    }
};
