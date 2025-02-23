
#pragma once

#include <spdlog/spdlog.h>

#include <sstream>

#ifdef _WIN32
#include <windows.h>
#include <dbghelp.h>
#pragma comment(lib, "Dbghelp.lib")
#else
#include <execinfo.h>
#include <dlfcn.h>
#include <cxxabi.h>
#include <unistd.h>
#endif

#include <cpptrace/cpptrace.hpp>
#include <cpptrace/formatting.hpp>

class tLog
{
private:
    static std::shared_ptr<spdlog::logger> logger_console;
    static std::shared_ptr<spdlog::logger> logger_file;

    std::string logPath;
    std::string traceFile;
    std::string errorFile;

	bool isConsole = true;
	bool isInfo = false;
	bool isError = false;
    static bool iWarnTrace;

private:
    static std::ostringstream log_stack_trace();

public:
    tLog(const std::string path = "logs", const std::string traceFile = "trace.log", const std::string errorFile = "error.log") :
        isConsole(true), isInfo(false), isError(false)
    {
        logPath = path;
        this->traceFile = traceFile;
        this->errorFile = errorFile;
    }

    ~tLog() {
        // shutdown all registered loggers
        spdlog::shutdown();
    #ifdef _WIN32
        SymCleanup(GetCurrentProcess());
    #endif
    }

    void set_path(const std::string& path) {
        logPath = path;
    }

    void set_info(bool flag = true, const std::string& file = "trace.log") {
        isInfo = flag;
        if (isInfo)
            traceFile = file;
    }

    void set_warn_trace(bool flag = true) {
        iWarnTrace = flag;
    }

    void set_error(bool flag = true, const std::string& file = "error.log") {
        isError = flag;
        if (isError)
            errorFile = file;
    }

	// initialize the logger
    void init();

	// log functions
    template<typename... Args>
    static void Trace(const std::string& msg, Args&&... args) {
        logger_console->trace(msg, std::forward<Args>(args)...);
        logger_file->trace(msg, std::forward<Args>(args)...);
    }

    template<typename... Args>
    static void Debug(const std::string& msg, Args&&... args) {
        logger_console->debug(msg, std::forward<Args>(args)...);
        logger_file->debug(msg, std::forward<Args>(args)...);
    }

    template<typename... Args>
    static void Info(const std::string& msg, Args&&... args) {
        logger_console->info(msg, std::forward<Args>(args)...);
        logger_file->info(msg, std::forward<Args>(args)...);
    }

    template<typename... Args>
    static void Warn(const std::string& msg, Args&&... args) {
        logger_console->warn(msg, std::forward<Args>(args)...);
        logger_file->warn(msg, std::forward<Args>(args)...);

        //std::ostringstream oss = log_stack_trace();
        //logger_console->warn(oss.str());

        if (iWarnTrace) {
            auto formatter = cpptrace::formatter{}
                .header("Stack trace (most recent call first):")
                .addresses(cpptrace::formatter::address_mode::object)
                .paths(cpptrace::formatter::path_mode::full)
                .snippets(true);

            logger_console->warn(formatter.format(cpptrace::generate_trace(1), true));
            logger_file->warn(formatter.format(cpptrace::generate_trace(1), false));
        }
    }

    template<typename... Args>
    static void Error(const std::string& msg, Args&&... args) {
        logger_console->error(msg, std::forward<Args>(args)...);
        logger_file->error(msg, std::forward<Args>(args)...);

        //std::ostringstream oss = log_stack_trace();
        //logger_console->error(oss.str());

        auto formatter = cpptrace::formatter{}
            .header("Stack trace (most recent call first):")
            .addresses(cpptrace::formatter::address_mode::object)
			.paths(cpptrace::formatter::path_mode::full)
            .snippets(true);

        logger_console->error(formatter.format(cpptrace::generate_trace(1), true));
        logger_file->error(formatter.format(cpptrace::generate_trace(1), false));
    }

};
