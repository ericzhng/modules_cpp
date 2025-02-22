import asyncio
import websockets

async def send_command():
    uri = "ws://localhost:9002"
    async with websockets.connect(uri) as websocket:
        command = "echo Hello, World!"
        await websocket.send(command)

        result = await websocket.recv()
        print(f"Command output: {result}")

asyncio.get_event_loop().run_until_complete(send_command())
