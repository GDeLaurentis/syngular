import re
import subprocess
import syngular


def execute_singular_command(singular_command, timeout=syngular.TIMEOUT, verbose=False):
    if isinstance(singular_command, list):
        singular_command = "\n".join(singular_command)
    # print(singular_command)
    if len(singular_command) >= 131072:  # 2 ** 17
        with open("/tmp/.singular_command", "w") as file:
            file.write(singular_command)
        test = subprocess.Popen(["timeout", str(timeout), "Singular", "--quiet", "/tmp/.singular_command"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    else:
        test = subprocess.Popen(["timeout", str(timeout), "Singular", "--quiet", "--execute", singular_command], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    output = test.communicate()[0]
    output = output.decode("utf-8")
    if output[-1] == "\n":
        output = output[:-1]
    if 'halt' in output:
        raise TimeoutError(f"{timeout} s")
    if 'error' in output and 'groebner base computations with inexact coefficients can not be trusted due to rounding errors' not in output:
        raise SingularException(output)
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
    return output


class SingularException(Exception):
    pass
