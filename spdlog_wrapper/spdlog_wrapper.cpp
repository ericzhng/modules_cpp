
#include "Logger.h"

int main() {
    // Initialize console logger
    Logger::InitConsoleLogger();

    // Initialize file loggers
    Logger::InitFileLoggers("logs/info.log", "logs/error.log");

    // Logging messages
    Logger::Info("Application started.");
    Logger::Warn("This is a warning message.");
    Logger::Error("An error occurred: {}", 404);
    Logger::Debug("This is a debug message.");
    Logger::Trace("This is a trace message.");

    return 0;
}
