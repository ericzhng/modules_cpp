@echo off
REM Set Anaconda base path (modify as needed)
SET ANACONDA_PATH=C:\Software\anaconda3

REM Initialize Conda (only needed if 'conda' is not in PATH)
CALL "%ANACONDA_PATH%\Scripts\activate.bat" "%ANACONDA_PATH%"

REM Activate the virtual environment (replace 'myenv' with your environment name)
CALL conda activate torch

REM Run the Python script
python python_reader.py

REM (Optional) Deactivate the environment
CALL conda deactivate
