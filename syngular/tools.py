import re
import subprocess
import syngular
import random

from pathlib import Path


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
        test = subprocess.Popen(["timeout", str(timeout), "Singular", "--quiet", file_path], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    else:
        test = subprocess.Popen(["timeout", str(timeout), "Singular", "--quiet", "--execute", singular_command], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    output = test.communicate()[0]
    output = output.decode("utf-8")
    if output[-1] == "\n":
        output = output[:-1]
    if 'halt' in output:
        raise TimeoutError(f"{timeout} s")
    if 'error' in output and 'groebner base computations with inexact coefficients can not be trusted due to rounding errors' not in output:
        raise SingularException(f"{output}\n\n\nError occured while executing:\n{singular_command}")
    if '   ? `Q` is undefined' in output:
        output = output.replace("   ? `Q` is undefined", "")
        print("Singular Error: Q is undefined")
    if '//options: redSB degBound redefine usage prompt\n' in output:
        output = output.replace("//options: redSB degBound redefine usage prompt\n", "")
        print("Singular Warning: redSB degBound redefine usage prompt")
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
    return output


class SingularException(Exception):
    pass


class RootNotInFieldError(Exception):
    pass


class RootPrecisionError(Exception):
    pass
