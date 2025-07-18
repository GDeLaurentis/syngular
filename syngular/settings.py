import importlib
import functools


class TemporarySetting:
    def __init__(self, module_or_module_name, setting_name, new_value):
        if isinstance(module_or_module_name, str):
            self.module = importlib.import_module(module_or_module_name)
        else:
            self.module = module_or_module_name
        self.setting_name = setting_name
        self.new_value = new_value
        self.old_value = None

    def __enter__(self):
        if hasattr(self.module, self.setting_name):
            self.old_value = getattr(self.module, self.setting_name)
            setattr(self.module, self.setting_name, self.new_value)
        else:
            raise AttributeError(f"Setting '{self.setting_name}' does not exist in module '{self.module.__name__}'")

    def __exit__(self, exc_type, exc_val, exc_tb):
        if hasattr(self.module, self.setting_name):
            setattr(self.module, self.setting_name, self.old_value)


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
