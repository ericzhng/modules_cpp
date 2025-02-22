
#include "PipeWriter.h"


int main()
{
    const char* pipePath =
#ifdef _WIN32
        R"(\\.\pipe\my_named_pipe)";
#else
        "/tmp/my_named_pipe";
#endif

    NamedPipeWriter pipeWriter(pipePath);

    // Failed to create or connect the pipe
    if (!pipeWriter.start()) {
        std::cerr << "Failed to initialize named pipe. Exiting.\n";
        return -1;
    }

    for (int i = 1; i <= 100; ++i) {

        std::string msg = std::to_string(i) + "\n";
        if (!pipeWriter.write(msg)) {
            std::cerr << "Write failed. Stopping.\n";
            break;
        }
        std::cout << "Sent: " << i << std::endl;
        std::this_thread::sleep_for(std::chrono::seconds(1));
    }

    return 0;
}
