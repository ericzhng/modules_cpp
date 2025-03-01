
## windows compiler
"C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.43.34808\bin\Hostx64\x64\dumpbin.exe" /exports mycpp.dll


## SWIG
"C:\Users\ZhangHui\AppData\Roaming\Python\Python312\Scripts\swig.exe" -python -c++ mycpp.i


The package python3 is compatible with built-in CMake targets:

    find_package(Python3 COMPONENTS Development REQUIRED)
    target_link_libraries(main PRIVATE Python3::Python)

The package python3 provides a python interpreter that supports virtual environments:

    >tools\python3\python.exe -m venv c:\path\to\venv
    >set VIRTUAL_ENV=c:\path\to\venv
    >set PATH=c:\path\to\venv\Scripts;%PATH%
    >set PYTHONHOME=

    See https://docs.python.org/3/library/venv.html for more details.

## venv

C:\dev\vcpkg\installed\x64-windows\tools\python3\python.exe -m venv D:\2-code\modules_cpp\extern\venv-python
set VIRTUAL_ENV=D:\2-code\modules_cpp\extern\venv-python
set PATH=D:\2-code\modules_cpp\extern\venv-python\Scripts;%PATH%
set PYTHONHOME=


D:\2-code\1-cpp-modules\extern\venv-python312\Scripts\activate.bat




