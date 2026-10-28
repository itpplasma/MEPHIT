"""Exercise delayed vector framing independently and with real BinaryIO."""
import os
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import threading
import time

client, interpreter, script = sys.argv[1:]
with tempfile.TemporaryDirectory() as directory:
    fifo = Path(directory) / "messages"
    os.mkfifo(fifo)
    errors = []
    long_size = struct.calcsize("@l")
    def read_flag():
        with fifo.open("rb", buffering=0) as stream:
            return struct.unpack("@l", stream.read(long_size))[0]
    def write_message(data, delay=False):
        with fifo.open("wb", buffering=0) as stream:
            stream.write(data)
            # Keep the completed writer open while the client consumes the
            # vector. Its acknowledgment must precede opening another reply.
            if delay:
                time.sleep(0.03)
        if delay:
            time.sleep(0.03)
    def delayed_peer():
        try:
            for _ in range(4):
                write_message(struct.pack("@l", read_flag()))
            for cycle in range(100):
                for state in (-1, -2):
                    assert read_flag() == state
                    write_message(struct.pack("@l", state))
                    with fifo.open("rb", buffering=0) as stream:
                        assert stream.read(long_size) == struct.pack("@l", 4)
                        data = stream.read(32)
                        assert len(data) == 32
                    for reply_index in range(3 if state == -1 else 1):
                        reply = struct.pack("@l", 4) + data if state == -1 else struct.pack("@ld", 1, 7.25)
                        write_message(reply, delay=(cycle == 0))
                        if state == -1 and reply_index < 2:
                            assert read_flag() == 0
            assert read_flag() == -3
            write_message(struct.pack("@l", -3))
        except Exception as error:
            errors.append(error)
    peer = threading.Thread(target=delayed_peer, daemon=True)
    peer.start()
    try:
        result = subprocess.run([client, "protocol", str(fifo)], capture_output=True,
                                text=True, timeout=10)
    except subprocess.TimeoutExpired as error:
        raise AssertionError(f"delayed peer errors: {errors}") from error
    peer.join(timeout=1)
    assert result.returncode == 0, result.stderr + repr(errors)
    assert not peer.is_alive() and not errors, errors
    daemon = subprocess.Popen([interpreter, "-nw", "-v", "0", script, "-P", str(fifo)],
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    try:
        result = subprocess.run([client, "protocol", str(fifo)], capture_output=True,
                                text=True, timeout=10)
        output, _ = daemon.communicate(timeout=5)
        assert result.returncode == 0, result.stderr + output
        assert daemon.returncode == 0, output
    finally:
        if daemon.poll() is None:
            daemon.kill()
            daemon.wait()
print("delayed vector framing + 100 actual BinaryIO scalar/vector cycles: PASS")
