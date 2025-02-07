import importlib


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
