import re
import subprocess
import syngular
import random
import warnings

from pathlib import Path
from packaging.version import Version, InvalidVersion


def execute_singular_command(singular_command, timeout='default', verbose=False):
    if timeout == 'default':
        timeout = syngular.TIMEOUT
    if isinstance(singular_command, list):
        singular_command = "\n".join(singular_command)
    # print(singular_command)
    if len(singular_command) >= 131072:  # 2 ** 17
        Path("/tmp/.singular_commands/").mkdir(parents=True, exist_ok=True)
        random_integer = random.randint(0, 2**64 - 1)
        file_path = f"/tmp/.singular_commands/singular_command_{random_integer}"
        with open(file_path, "w+") as file:
            file.write(singular_command)
        test = subprocess.Popen(["timeout", "--verbose", str(timeout), "Singular", "--quiet", file_path],
                                stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    else:
        test = subprocess.Popen(["timeout", "--verbose", str(timeout), "Singular", "--quiet", "--execute", singular_command],
                                stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    try:
        output, stderr = test.communicate()
    except KeyboardInterrupt:
        print("Keyboard interrupt received. Terminating the Singular process.")
        if test.poll() is None:  # Check if the process is still running
            import os
            import signal
            os.kill(test.pid, signal.SIGTERM)
        raise KeyboardInterrupt
    output = output.decode("utf-8")
    if stderr is not None:
        stderr = stderr.decode("utf-8")
    if len(output) == 0 and stderr is None:
        raise SingularException(f"Empty output while executing: {singular_command}")
    if len(output) > 0 and output[-1] == "\n":
        output = output[:-1]
    if 'halt' in output or (stderr is not None and 'timeout' in stderr):
        raise TimeoutError(f"{timeout} s")
    if '   ? `Q` is undefined' in output:
        output = output.replace("   ? `Q` is undefined", "")
        print("Singular Error: Q is undefined")
    if (('error' in output and 'groebner base computations with inexact coefficients can not be trusted due to rounding errors' not in output) or
       'is not defined' in output or 'INVALID' in output or 'can not convert' in output):
        raise SingularException(f"{output}\n\n\nError occured while executing:\n{singular_command}")
    if '//options: redSB degBound redefine usage prompt\n' in output:
        output = output.replace("//options: redSB degBound redefine usage prompt\n", "")
        print("Singular Warning: redSB degBound redefine usage prompt")
    if '//options: degBound redTail redThrough redefine usage prompt\n' in output:
        output = output.replace("//options: degBound redTail redThrough redefine usage prompt\n", "")
        print("Singular Warning: degBound redTail redThrough redefine usage prompt")
    output = output.replace('// ** minbase applies only to the local or homogeneous case over coefficient fields\n', '')
    output = output.replace('//options: redSB degBound redTail redThrough intStrategy redefine usage prompt\n', '')
    output = output.replace('//options: degBound redTail redThrough intStrategy redefine usage prompt\n', '')
    output = re.sub(r'// \*\* .* is no standard basis\n', '', output)
    if 'groebner base computations with inexact coefficients can not be trusted due to rounding errors' in output:
        if verbose:
            print('Singular Warning: groebner base computations with inexact coefficients can not be trusted due to rounding errors')
        output = output.replace('// ** groebner base computations with inexact coefficients can not be trusted due to rounding errors\n', '')
    if len(singular_command) >= 131072:
        Path(file_path).unlink()
    if syngular.DEBUG:
        print("DEBUG - Command received:\n", singular_command)
        print("DEBUG - Output obtained:\n", output)
    return output


class SingularException(Exception):
    pass


class RootNotInFieldError(Exception):
    pass


class RootPrecisionError(Exception):
    pass


try:
    test = subprocess.Popen(["timeout", "5", "Singular", "--dump-versiontuple"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    output = test.communicate()[0]
    output = output.decode("utf-8").replace("\n", "")
    try:
        Singular_version = Version(output)
    except InvalidVersion:
        warnings.warn("Could not determine Singular version. Are you sure Singular is installed?", stacklevel=2)
        Singular_version = None
except FileNotFoundError as e:
    warnings.warn(f"\nAre you sure timeout is installed for use in the command line?\n{e}\nCould not determine Singular version.", stacklevel=2)
    Singular_version = None
