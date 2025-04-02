import importlib
import functools


class TemporarySetting:
    def __init__(self, module_name, setting_name, new_value):
        self.module_name = module_name
        self.setting_name = setting_name
        self.new_value = new_value
        self.old_value = None

    def __enter__(self):
        module = importlib.import_module(self.module_name)
        if hasattr(module, self.setting_name):
            self.old_value = getattr(module, self.setting_name)
            setattr(module, self.setting_name, self.new_value)
        else:
            raise AttributeError(f"Setting '{self.setting_name}' does not exist in module '{self.module_name}'")

    def __exit__(self, exc_type, exc_val, exc_tb):
        module = importlib.import_module(self.module_name)
        if hasattr(module, self.setting_name):
            setattr(module, self.setting_name, self.old_value)


def with_cm(cm):
    """
    A decorator that applies a context manager to the function.
    This allows a `with` statement to be used as a decorator.
    """
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            with cm:
                return func(*args, **kwargs)
        return wrapper
    return decorator


def with_other_cas_compatible_str(func):
    @with_cm(TemporarySetting("syngular", "FORCECDOTS", True))
    @with_cm(TemporarySetting("syngular", "CDOTCHAR", "*"))
    @with_cm(TemporarySetting("syngular", "UNICODEPOWERS", False))
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)
    return wrapper
