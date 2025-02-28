
## windows compiler
"C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.43.34808\bin\Hostx64\x64\dumpbin.exe" /dependents mycpp.dll

## SWIG
"C:\Users\ZhangHui\AppData\Roaming\Python\Python312\Scripts\swig.exe" -python -c++ mycpp.i


## pybind11, swig, boost-python
all above still have some issues, python code cannot find the python package
