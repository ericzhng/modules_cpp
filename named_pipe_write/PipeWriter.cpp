
#include "PipeWriter.h"

bool NamedPipeWriter::start(int max_retries, int retry_delay_ms)
{
    if (max_retries <= 0 || retry_delay_ms <= 0) {
        std::cerr << "Invalid arguments: max_retries and retry_delay_ms must be positive integers.\n";
        return false;
    }

#ifdef _WIN32
    std::wstring wPipeName = convertToWideString(pipeName);

    int attempts = 0;
    while (attempts < max_retries)
    {
        // Create named pipe
        hPipe = CreateNamedPipe(
            wPipeName.c_str(),
            PIPE_ACCESS_OUTBOUND,
            PIPE_TYPE_BYTE | PIPE_READMODE_BYTE | PIPE_WAIT,
            1, 0, 0, 0, NULL
        );

        if (hPipe == INVALID_HANDLE_VALUE) {

            std::cerr << "Failed to create named pipe. Error: " << GetLastError()
                << " (Attempt " << (attempts + 1) << "/" << max_retries << ")\n";

            std::this_thread::sleep_for(std::chrono::milliseconds(retry_delay_ms));
            attempts++;
        }
        else {
            break;
        }
    }

    if (hPipe == NULL || hPipe == INVALID_HANDLE_VALUE) {
        std::cerr << "Failed to connect named pipe" << std::endl;
        return false;
    }

    std::cout << "Waiting for client to connect...\n";

    while (true) {
        BOOL connected = ConnectNamedPipe(hPipe, NULL);
        if (connected) {
            break;  // Successfully connected
        }

        DWORD err = GetLastError();
        if (err == ERROR_PIPE_CONNECTED) {
            break;  // Client connected before `ConnectNamedPipe`
        }
        else if (err == ERROR_NO_DATA || err == ERROR_BROKEN_PIPE) {
            std::cerr << "Client disconnected. Waiting for reconnection...\n";
            DisconnectNamedPipe(hPipe);
            continue;
        }
        else {
            std::cerr << "Failed to connect named pipe. Error: " << err << std::endl;
            CloseHandle(hPipe);
            return false;
        }
    }


#else
    if (mkfifo(pipeName.c_str(), 0666) == -1 && errno != EEXIST) {
        perror("Failed to create FIFO");
        return false;
    }

    int attempts = 0;
    while (attempts < max_retries)
    {
        fd = open(pipeName.c_str(), O_WRONLY);
        if (fd == -1) {
            std::cerr << "Failed to open FIFO for writing. Error: " << strerror(errno)
                << " (Attempt " << (attempts + 1) << "/" << max_retries << ")\n";
            std::this_thread::sleep_for(std::chrono::milliseconds(retry_delay_ms));
            attempts++;
        }
        else {
            break;
        }
    }

    if (fd == -1) {
        return false;
    }
#endif

    std::cout << "Client connected. Ready to write.\n";
    return true;
}


bool NamedPipeWriter::write(const std::string& message)
{
    if (message.empty()) {
        return true;
    }

#ifdef _WIN32
    if (hPipe == INVALID_HANDLE_VALUE) {
        std::cerr << "Write failed: Pipe is not initialized!\n";
        return false;
    }

    DWORD bytesWritten;
    size_t msgSize = message.size();

    if (msgSize > MAXDWORD) {  // MAXDWORD = 0xFFFFFFFF (largest possible DWORD)
        std::cerr << "Write failed: Message too large for WriteFile()!\n";
        return false;
    }

    BOOL success = WriteFile(hPipe, message.c_str(), static_cast<DWORD>(msgSize), &bytesWritten, NULL);
    if (!success)
    {
        DWORD err = GetLastError();
        if (err == ERROR_BROKEN_PIPE) {
            std::cerr << "Client disconnected. Reconnecting...\n";
            start();  // Attempt to reconnect
            return false;
        }
        std::cerr << "Write failed. Error: " << err << std::endl;
        return false;
    }
#else
    if (fd == -1) {
        std::cerr << "Write failed: FIFO is not initialized!\n";
        return false;
    }

    ssize_t written = ::write(fd, message.c_str(), message.size());
    if (written == -1) {
        std::cerr << "Write failed. Error: " << strerror(errno) << std::endl;
        return false;
    }
#endif

    return true;
}

void NamedPipeWriter::close()
{
    std::cout << "Closing pipe...\n";

#ifdef _WIN32
    if (hPipe != NULL && hPipe != INVALID_HANDLE_VALUE) {
        FlushFileBuffers(hPipe);
        DisconnectNamedPipe(hPipe);
        CloseHandle(hPipe);
        hPipe = INVALID_HANDLE_VALUE;
    }
#else
    if (fd != -1) {
        ::close(fd);
        fd = -1;
    }
#endif
}
