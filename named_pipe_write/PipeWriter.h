
#pragma once

#include <iostream>
#include <string>
#include <thread>
#include <stdexcept>

#ifdef _WIN32
    #include <windows.h>
    #include <stringapiset.h>  // For MultiByteToWideChar
#else
    #include <fcntl.h>
    #include <unistd.h>
    #include <sys/stat.h>
    #include <errno.h>
#endif


class NamedPipeWriter {

private:
    std::string pipeName;

#ifdef _WIN32
    HANDLE hPipe;
#else
    int fd;
#endif

#ifdef _WIN32

    std::wstring convertToWideString(const std::string& str)
    {
        int size_needed = MultiByteToWideChar(CP_UTF8, 0, str.c_str(), -1, NULL, 0);
        std::wstring wstr(size_needed, 0);
        MultiByteToWideChar(CP_UTF8, 0, str.c_str(), -1, &wstr[0], size_needed);
        return wstr;
    }

#endif


public:
    explicit NamedPipeWriter(const std::string& name)
        : pipeName(name)
#ifdef _WIN32
        , hPipe(INVALID_HANDLE_VALUE)
#else
        , fd(-1)
#endif
    {
    };

    ~NamedPipeWriter() {
        close();
    };

    bool start(int max_retries = 5, int retry_delay_ms = 1000);

    bool write(const std::string& message);

    void close();
};
