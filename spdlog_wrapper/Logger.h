#pragma once

#define FMT_UNICODE 0

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <memory>
#include <vector>


class Logger {
public:
    static void InitConsoleLogger() {
        auto consoleSink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
        consoleSink->set_pattern("[%^%l%$] %Y-%m-%d %H:%M:%S %v");
        consoleLogger = std::make_shared<spdlog::logger>("console_logger", consoleSink);
        consoleLogger->set_level(spdlog::level::info);
        spdlog::register_logger(consoleLogger);
    }

    static void InitFileLoggers(const std::string& infoLogFile = "logs/info.log", const std::string& errorLogFile = "logs/error.log") {
        auto infoFileSink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(infoLogFile, true);
        infoFileSink->set_pattern("[%l] %Y-%m-%d %H:%M:%S %v");

        auto errorFileSink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(errorLogFile, true);
        errorFileSink->set_pattern("[%l] %Y-%m-%d %H:%M:%S %v");
        errorFileSink->set_level(spdlog::level::err);

        infoLogger = std::make_shared<spdlog::logger>("info_logger", infoFileSink);
        errorLogger = std::make_shared<spdlog::logger>("error_logger", errorFileSink);

        infoLogger->set_level(spdlog::level::info);
        errorLogger->set_level(spdlog::level::err);

        spdlog::register_logger(infoLogger);
        spdlog::register_logger(errorLogger);
    }

    template<typename... Args>
    static void Info(const std::string& msg, Args&&... args) {
        consoleLogger->info(msg, std::forward<Args>(args)...);
        infoLogger->info(msg, std::forward<Args>(args)...);
    }

    template<typename... Args>
    static void Warn(const std::string& msg, Args&&... args) {
        consoleLogger->warn(msg, std::forward<Args>(args)...);
        infoLogger->warn(msg, std::forward<Args>(args)...);
    }

    template<typename... Args>
    static void Error(const std::string& msg, Args&&... args) {
        consoleLogger->error(msg, std::forward<Args>(args)...);
        errorLogger->error(msg, std::forward<Args>(args)...);
    }

    template<typename... Args>
    static void Debug(const std::string& msg, Args&&... args) {
        consoleLogger->debug(msg, std::forward<Args>(args)...);
    }

    template<typename... Args>
    static void Trace(const std::string& msg, Args&&... args) {
        consoleLogger->trace(msg, std::forward<Args>(args)...);
    }

    static void SetLevel(spdlog::level::level_enum level) {
        consoleLogger->set_level(level);
        infoLogger->set_level(level);
        errorLogger->set_level(level);
    }

private:
    static std::shared_ptr<spdlog::logger> consoleLogger;
    static std::shared_ptr<spdlog::logger> infoLogger;
    static std::shared_ptr<spdlog::logger> errorLogger;
};

std::shared_ptr<spdlog::logger> Logger::consoleLogger = nullptr;
std::shared_ptr<spdlog::logger> Logger::infoLogger = nullptr;
std::shared_ptr<spdlog::logger> Logger::errorLogger = nullptr;
