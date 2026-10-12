"""Signal one rank or all ranks; require one message and a nonzero exit."""
import os
import signal
import subprocess
import sys

launcher, executable, count, mode, signal_name = sys.argv[1:]
stop_signal = getattr(signal, signal_name)
expected_message = {
    'SIGINT': 'Simulation interrupted by Ctrl+C.',
    'SIGTERM': 'Simulation stopped after SIGTERM (scheduler termination or cancellation).',
    'SIGUSR1': 'Simulation stopped after SIGUSR1 (scheduler warning).',
}[signal_name]
count = int(count)
process = subprocess.Popen(
    [launcher, '-n', str(count), executable],
    stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
)
try:
    ranks = {}
    while len(ranks) < count:
        line = process.stdout.readline()
        if not line:
            raise AssertionError('MPI ranks exited before becoming ready: ' + process.stderr.read())
        if line.startswith('READY '):
            _, rank, pid = line.split()
            ranks[int(rank)] = int(pid)
    targets = ranks.values() if mode == 'all' else [ranks[count - 1]]
    for pid in targets:
        os.kill(pid, stop_signal)
    stdout, stderr = process.communicate(timeout=15)
    output = stdout + stderr
    assert output.count(expected_message) == 1, output
    assert process.returncode == 128 + stop_signal, output
    assert 'forrtl:' not in output and 'Backtrace' not in output, output
    print(f'{signal_name}, {count} ranks, {mode}: one stop message; exit {process.returncode}')
finally:
    if process.poll() is None:
        process.kill()
        process.communicate()
