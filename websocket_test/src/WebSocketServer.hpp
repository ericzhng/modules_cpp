
#ifndef WEBSOCKET_SERVER_HPP
#define WEBSOCKET_SERVER_HPP

#include <websocketpp/config/asio_no_tls.hpp>
#include <websocketpp/server.hpp>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <array>
#include <memory>
#include <thread>

typedef websocketpp::server<websocketpp::config::asio> server;

class WebSocketServer {
public:
    WebSocketServer();
    ~WebSocketServer();

    void start();
    void stop();

private:
    void execute_command(const std::string& cmd, std::string& result);
    void on_message(websocketpp::connection_hdl hdl, server::message_ptr msg);

    server m_server;
    std::thread m_thread;
    bool m_running;

    static void run_server(WebSocketServer* instance);
};

#endif // WEBSOCKET_SERVER_HPP