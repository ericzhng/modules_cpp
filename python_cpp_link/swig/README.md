
## SWIG
D:\2-code\1-cpp-modules\extern\venv-python312\Scripts\activate.bat
pip install swig
swig -python -c++ mycpp.i

You can only build in Release mode; Debug mode generate errors of cannot find python312_d.lib.

After compilation, you need to copy both _mycpp.pyd and mycpp.py files to the python script directory. Then you can execute test.py

### Note
For the 2nd method of add_custom_command, you need to modify name of the final pyd file to match with the name of the python file, adding a _ ahead.



## use following to check library
"C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.43.34808\bin\Hostx64\x64\dumpbin.exe" /exports build/Release/_mycpp.pyd
