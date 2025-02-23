@echo off
REM Set Anaconda base path (modify as needed)
SET ANACONDA_PATH=C:\Software\anaconda3

REM Initialize Conda (only needed if 'conda' is not in PATH)
CALL "%ANACONDA_PATH%\Scripts\activate.bat" "%ANACONDA_PATH%"

REM Activate the virtual environment (replace 'myenv' with your environment name)
REM CALL conda activate env_test
CALL conda activate env_pipe

REM Run the Python script
python python/client.py

REM (Optional) Deactivate the environment
CALL conda deactivate
