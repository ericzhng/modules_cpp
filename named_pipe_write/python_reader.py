
import os
import time

pipe_name = r'\\.\pipe\my_named_pipe' if os.name == 'nt' else '/tmp/my_named_pipe'

print("Waiting for the pipe connection...")

while True:
    try:
        with open(pipe_name, 'r') as pipe:
            print("Connected to the pipe. Reading data...\n")
            while True:
                line = pipe.readline().strip()
                if not line:
                    print("Pipe closed. Exiting...")
                    break
                print(f"Received: {line}")
    except FileNotFoundError:
        time.sleep(1)
