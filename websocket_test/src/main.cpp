
#include "WebSocketServer.hpp"
#include <iostream>

int main() {
    WebSocketServer ws_server;

    std::cout << "Starting WebSocket server..." << std::endl;
    ws_server.start();
    std::cout << "WebSocket server started. Press Enter to stop..." << std::endl;

    // Wait for the user to press Enter
    std::cin.get();

    std::cout << "Stopping WebSocket server..." << std::endl;
    ws_server.stop();
    std::cout << "WebSocket server stopped." << std::endl;

    return 0;
}
