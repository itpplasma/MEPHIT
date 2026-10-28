"""Exercise supervisor exit attribution and actual truncated pipe messages."""
import os
from pathlib import Path
import struct
import subprocess
import sys
import tempfile


def run(command):
    return subprocess.run(command, capture_output=True, text=True, timeout=5)


supervisor, client = sys.argv[1:]
with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    daemon = root / "daemon.sh"
    for config, body, diagnostic in [
        ("wait", "exit 8", "FreeFem exited with error code 8"),
        ("wait", "kill -KILL $$", "FreeFem terminated abnormally"),
        ("fail", "exec sleep 30", "MEPHIT exited with error code 23"),
        ("ok", "exit 0", ""),
    ]:
        daemon.write_text("#!/bin/sh\n" + body + "\n")
        daemon.chmod(0o700)
        result = run([supervisor, "0", config, "", directory, str(daemon)])
        assert result.returncode == (1 if diagnostic else 0), result
        assert diagnostic in result.stderr, result.stderr
        if config == "wait":
            assert "FreeFem failed; stopping MEPHIT" in result.stderr, result.stderr
            assert "freefem.log" in result.stderr, result.stderr
        assert not list(root.glob("MEPHIT_*.dat")), "supervisor left its FIFO"

    data = root / "data"
    for kind, content, succeeds in [
        ("scalar", struct.pack("@l", 42), True),
        ("scalar", b"", False),
        ("scalar", b"abc", False),
        ("vector", struct.pack("@l2d", 2, 1.25, -2.5), True),
        ("vector", struct.pack("@l", 2) + b"abc", False),
    ]:
        data.write_bytes(content)
        result = run([client, kind, str(data)])
        assert (result.returncode == 0) == succeeds, result
        if not succeeds:
            assert "FreeFem communication failed" in result.stderr, result.stderr

    fifo = root / "broken_fifo"
    os.mkfifo(fifo)
    writer = subprocess.Popen([client, "write", str(fifo)], stderr=subprocess.PIPE,
                              text=True)
    with fifo.open("rb", buffering=0) as reader:
        assert len(reader.read(struct.calcsize("@l"))) == struct.calcsize("@l")
    _, error = writer.communicate(timeout=5)
    assert writer.returncode == 1, (writer.returncode, error)
    assert "FreeFem communication failed while writing" in error, error
print("supervisor exits/signals, valid/truncated transfers, broken pipe: PASS")
