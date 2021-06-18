import subprocess


def execute_singular_command(singular_command, timeout=60):
    if isinstance(singular_command, list):
        singular_command = "\n".join(singular_command)
    if len(singular_command) >= 131072:  # 2 ** 17
        with open("/tmp/.singular_command", "w") as file:
            file.write(singular_command)
        test = subprocess.Popen(["timeout", str(timeout), "Singular", "--quiet", "/tmp/.singular_command"], stdout=subprocess.PIPE)
    else:
        test = subprocess.Popen(["timeout", str(timeout), "Singular", "--quiet", "--execute", singular_command], stdout=subprocess.PIPE)
    output = test.communicate()[0]
    output = output.decode("utf-8")
    if output[-1] == "\n":
        output = output[:-1]
    if 'halt' in output:
        raise TimeoutError
    if 'error' in output:
        raise SingularException(output)
    return output


class SingularException(Exception):
    pass
