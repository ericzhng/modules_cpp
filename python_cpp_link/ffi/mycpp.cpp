#include <iostream>

extern "C" {
    __declspec(dllexport) void hello(int x) {
        std::cout << "Hello from C++! Received: " << x << std::endl;
    }
}
