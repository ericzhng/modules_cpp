import asyncio
import websockets
import logging
import sys

# If no event loop exists, create and set one
loop = asyncio.new_event_loop()
asyncio.set_event_loop(loop)

async def send_echo(websocket, i: int):
    try:
        message = f"echo {i}"
        print(f"Sending message: {message}")
        await websocket.send(message)  # Send message to WebSocket
        response = await websocket.recv()  # Receive the echo response
        print(f"Received echo: {response}")
    except websockets.exceptions.ConnectionClosed as e:
        logging.error(f"WebSocket connection closed unexpectedly: {e}")
        sys.exit(1)  # Stop the program when connection is closed
    except Exception as e:
        logging.error(f"An error occurred: {e}")
        sys.exit(1)  # Stop the program when connection is closed

async def main():
    uri = "ws://localhost:9002"
    try:
        async with websockets.connect(uri) as websocket:
            for i in range(1, 101):
                await send_echo(websocket, i)
                await asyncio.sleep(0.1)  # Simulate small delay between commands
    except websockets.exceptions.InvalidURI as e:
        logging.error(f"Invalid WebSocket URI: {e}")
    except websockets.exceptions.WebSocketException as e:
        logging.error(f"WebSocket connection failed: {e}")
    except Exception as e:
        logging.error(f"Unexpected error in connection: {e}")

loop.run_until_complete(main())
