import re
import subprocess
import syngular
import random
import warnings
import shutil
import functools

from pathlib import Path
from packaging.version import Version, InvalidVersion


def execute_singular_command(singular_command, timeout='default', verbose=False):
    _check_singular_available()
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
    if syngular.DEBUG:
        if len(singular_command) >= 131072:
            print(f"DEBUG - Command of length {len(singular_command)} was written to file: {file_path}")
        else:
            print(f"DEBUG - Command of length {len(singular_command)} was executed directly.")
        print("DEBUG - Command received:\n", singular_command)
        print("DEBUG - Output obtained:\n", output)
        print("DEBUG - Stderr obtained:\n", stderr)
    if len(output) == 0 and stderr is None:
        raise SingularException(f"Empty output while executing: {singular_command}")
    if len(output) > 0 and output[-1] == "\n":
        output = output[:-1]
    if 'halt' in output or (stderr is not None and 'timeout' in stderr):
        raise TimeoutError(f"{timeout} s")
    if (('error' in output and 'groebner base computations with inexact coefficients can not be trusted due to rounding errors' not in output) or
       'is not defined' in output or 'is undefined' in output or 'INVALID' in output or 'can not convert' in output):
        raise SingularException(f"{output}\n\n\nError occured while executing:\n{singular_command}")
    option_lines = re.findall(r"^//options:[^\n]*(?:\n|$)", output, flags=re.MULTILINE)
    if syngular.DEBUG:
        for line in option_lines:
            print(f"Singular options: {line.removeprefix('//options:').strip()}")
    output = re.sub(r"^//options:[^\n]*(?:\n|$)", "", output, flags=re.MULTILINE)
    output = output.replace('// ** minbase applies only to the local or homogeneous case over coefficient fields\n', '')
    output = re.sub(r'// \*\* .* is no standard basis\n', '', output)
    if 'groebner base computations with inexact coefficients can not be trusted due to rounding errors' in output:
        if verbose:
            print('Singular Warning: groebner base computations with inexact coefficients can not be trusted due to rounding errors')
        output = output.replace('// ** groebner base computations with inexact coefficients can not be trusted due to rounding errors\n', '')
    if len(singular_command) >= 131072:
        Path(file_path).unlink()
    return output


class SingularException(Exception):
    pass


class RootNotInFieldError(Exception):
    pass


class RootPrecisionError(Exception):
    pass


def _check_singular_available(raise_error=True):
    if shutil.which("Singular") is not None:
        return True
    message = "Singular is required for this operation but was not found in PATH."
    if raise_error:
        raise FileNotFoundError(message)
    warnings.warn(message, stacklevel=2)
    return False


@functools.cache
def get_singular_version(raise_error=True):
    if not _check_singular_available(raise_error=raise_error):
        return None
    try:
        result = subprocess.run(["Singular", "--dump-versiontuple"], capture_output=True, text=True, timeout=5, check=True, )
        return Version(result.stdout.strip())
    except (subprocess.SubprocessError, InvalidVersion) as exception:
        if raise_error:
            raise
        warnings.warn(f"Could not determine Singular version, despite Singular being available: {exception}", stacklevel=2)
        return None


def __getattr__(name):
    if name == "Singular_version":
        warnings.warn(
            "Singular_version is deprecated; use get_singular_version() instead. "
            "Pass raise_error=False to retain the previous non-raising behaviour.",
            DeprecationWarning, stacklevel=2,
        )
        value = get_singular_version(raise_error=False)
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
