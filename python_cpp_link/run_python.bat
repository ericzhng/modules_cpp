@echo off
REM Set Anaconda base path (modify as needed)
SET ANACONDA_PATH=C:\Software\anaconda3

REM Initialize Conda (only needed if 'conda' is not in PATH)
CALL "%ANACONDA_PATH%\Scripts\activate.bat" "%ANACONDA_PATH%"

REM Activate the virtual environment (replace 'myenv' with your environment name)
CALL conda activate env_pipe

REM Run the Python script
python ctypes/test.py
python ffi/test.py
python pybind11/test.py
python boost-python/test.py
python swig/test.py

REM (Optional) Deactivate the environment
CALL conda deactivate
