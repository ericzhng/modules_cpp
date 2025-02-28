

# # Set the C source code (can be in a separate file)
# ffi.set_source("mycpp_", """
# #include <iostream>

# extern "C" {
    # __declspec(dllexport) void hello(int x) {
        # std::cout << "Hello from C++! Received: " << x << std::endl;
    # }
# }
    # """)
    
# # Compile the C code
# ffi.compile()

# import mycpp_
# mycpp_.hello(42)

