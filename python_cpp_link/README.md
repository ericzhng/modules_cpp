
## windows compiler
"C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.43.34808\bin\Hostx64\x64\dumpbin.exe" /dependents mycpp.dll

## SWIG
"C:\Users\ZhangHui\AppData\Roaming\Python\Python312\Scripts\swig.exe" -python -c++ mycpp.i


## pybind11, swig, boost-python
all above still have some issues, python code cannot find the python package


The package python3 is compatible with built-in CMake targets:

    find_package(Python3 COMPONENTS Development REQUIRED)
    target_link_libraries(main PRIVATE Python3::Python)

The package python3 provides a python interpreter that supports virtual environments:

    >tools\python3\python.exe -m venv c:\path\to\venv
    >set VIRTUAL_ENV=c:\path\to\venv
    >set PATH=c:\path\to\venv\Scripts;%PATH%
    >set PYTHONHOME=

    See https://docs.python.org/3/library/venv.html for more details.



C:\dev\vcpkg\installed\x64-windows\tools\python3\python.exe -m venv D:\2-code\modules_cpp\extern\venv-python
set VIRTUAL_ENV=D:\2-code\modules_cpp\extern\venv-python
set PATH=D:\2-code\modules_cpp\extern\venv-python\Scripts;%PATH%
set PYTHONHOME=


    find_package(Python COMPONENTS Interpreter Development)
    find_package(pybind11 CONFIG)

    # pybind11 method:
    pybind11_add_module(MyModule1 src1.cpp)

    # Python method:
    Python_add_library(MyModule2 src2.cpp)
    target_link_libraries(MyModule2 PRIVATE pybind11::headers)
    set_target_properties(MyModule2 PROPERTIES
        INTERPROCEDURAL_OPTIMIZATION ON
        CXX_VISIBILITY_PRESET ON
        VISIBILITY_INLINES_HIDDEN ON
    )



    find_package(Boost REQUIRED COMPONENTS python)
    target_link_libraries(main PRIVATE Boost::python)

or the generated cmake configs via:

    find_package(boost_python REQUIRED CONFIG)
    target_link_libraries(main PRIVATE Boost::python)

