
#include "spdlog_wrap.h"

void test1() {
    tLog::Error("An error occurred: {}", 404);
}

void test2() {
    tLog::Error("An warning occurred: {}", 404);
}

void test3() {
    tLog::Warn("An warning occurred: {}", 404);
}

void test_wrap() {
    test1();
    //test2();
    //test3();
}

int main() {

	tLog tLog("logs", "trace.log", "error.log");

    // Initialize file loggers
    tLog.init();
	tLog.set_warn_trace(true);

    // Logging messages
    tLog::Debug("This is a debug message.");
    tLog::Trace("This is a trace message.");
    tLog::Info("Application started.");

    test_wrap();

    return 0;
}
