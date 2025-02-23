
#include "spdlog_wrap.h"

#include <iomanip>

#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>


// initialize the logger
std::shared_ptr<spdlog::logger> tLog::logger_console = nullptr;
std::shared_ptr<spdlog::logger> tLog::logger_file = nullptr;
bool tLog::iWarnTrace = false;

void tLog::init()
{
	auto sink_console = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
	sink_console->set_pattern("[%^%l%$] %H:%M %v");
	sink_console->set_level(spdlog::level::debug);

	std::string sinkName;
	sinkName = logPath + "/" + traceFile;
	auto sink_trace = std::make_shared<spdlog::sinks::basic_file_sink_mt>(sinkName.c_str(), true);
	sink_trace->set_pattern("[%l] %Y-%m-%d %H:%M:%S %v");
	sink_trace->set_level(spdlog::level::trace);

	sinkName = logPath + "/" + errorFile;
	auto sink_error = std::make_shared<spdlog::sinks::basic_file_sink_mt>(sinkName.c_str(), true);
	sink_error->set_pattern("[%l] %Y-%m-%d %H:%M:%S %v");
	sink_error->set_level(spdlog::level::err);

	// add sinks to a vector and create the logger
	std::vector<spdlog::sink_ptr> sinks{ sink_console };
	logger_console = std::make_shared<spdlog::logger>("console_sink", begin(sinks), end(sinks));
	logger_console->set_level(spdlog::level::trace);   // set to be lower than the lowest level of the sinks
	spdlog::register_logger(logger_console);

	// add sinks to a vector and create the logger
	std::vector<spdlog::sink_ptr> sinks_file{ sink_trace, sink_error };
	logger_file = std::make_shared<spdlog::logger>("multi_sink", begin(sinks_file), end(sinks_file));
	logger_file->set_level(spdlog::level::trace);
	spdlog::register_logger(logger_file);

	// for native windows stack trace
#ifdef _WIN32
	// Initialize symbol handler for stack tracing
	SymSetOptions(SYMOPT_LOAD_LINES | SYMOPT_UNDNAME | SYMOPT_DEFERRED_LOADS);
	SymInitialize(GetCurrentProcess(), nullptr, TRUE);
#endif

}


std::ostringstream tLog::log_stack_trace()
{
	constexpr int maxFrames = 20;

	std::ostringstream oss;
	oss << "Traceback (from parent to child calls):\n";

#ifdef _WIN32
	void* stack[maxFrames];
	USHORT nframes = CaptureStackBackTrace(0, maxFrames, stack, nullptr);
	if (nframes == 0) {
		oss << "  Empty (failed to capture stack trace)\n";
		return oss;
	}

	HANDLE process = GetCurrentProcess();

	SYMBOL_INFO* symbol = (SYMBOL_INFO*)calloc(sizeof(SYMBOL_INFO) + 256 * sizeof(char), 1);
	if (symbol == NULL) {
		oss << "  Empty (failed to allocate memory for symbol info)\n";
		return oss;
	}
	symbol->MaxNameLen = 255;
	symbol->SizeOfStruct = sizeof(SYMBOL_INFO);

	DWORD displacement;

	IMAGEHLP_LINE64 line;
	line.SizeOfStruct = sizeof(IMAGEHLP_LINE64);
	std::string funcName = "UNK";
	std::string fileName = "UNK";
	int lineNumber = 0;

	// Start from 1 to skip log_stack_trace itself
	for (USHORT i = 2; i < nframes - 6; i++) {

		DWORD64 address = (DWORD64)(stack[i]);

		if (SymFromAddr(process, address, 0, symbol)) {
			funcName = symbol->Name;
			if (funcName.empty()) {
				funcName = "[UNK]";
			}
		}

		if (SymGetLineFromAddr64(process, address, &displacement, &line)) {
			fileName = line.FileName;
			lineNumber = line.LineNumber;
		}

		oss << "  #" << std::left << std::setw(2) << i - 1 << " 0x" << std::hex << address << " in ";
		oss << funcName << " at " << fileName << ":" << "" << std::dec << lineNumber << "\n";
	}
	free(symbol);

#else

	void* stack[maxFrames];
	int nframes = backtrace(stack, maxFrames);

	if (nframes == 0) {
		oss << "  Empty (failed to capture stack trace)\n";
		return oss;
	}

	char** symbols = backtrace_symbols(stack, nframes);
	if (!symbols) {
		oss << "  Empty (failed to retrieve symbols)\n";
		return oss;
	}

	// Start from 1 to skip log_stack_trace itself
	for (int i = 1; i < nframes; i++)
	{
		Dl_info info;
		if (dladdr(stack[i], &info) && info.dli_sname) {
			int status = 0;
			char* demangled = abi::__cxa_demangle(info.dli_sname, nullptr, nullptr, &status);
			oss << i << ": " << (status == 0 ? demangled : info.dli_sname)
				<< " - " << (info.dli_fname ? info.dli_fname : "[Unknown file]")
				<< " (0x" << stack[i] << ")\n";

			oss << "  #" << std::left << std::setw(2) << i - 1 << " 0x" << stack[i] << " in ";
			oss << symbols[i] << "\n";

			free(demangled);
		}
	}

	free(symbols);
#endif

	return oss;
}
