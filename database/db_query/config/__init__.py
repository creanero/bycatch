"""Expose base pipeline settings from the workspace-level config.py file."""

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

_CONFIG_PATH = Path(__file__).resolve().parent.parent / "config.py"
_SPEC = spec_from_file_location("_base_config", _CONFIG_PATH)

if _SPEC is None or _SPEC.loader is None:
	raise ImportError(f"Could not load configuration module at {_CONFIG_PATH}")

_MODULE = module_from_spec(_SPEC)
_SPEC.loader.exec_module(_MODULE)

# Re-export conventional config symbols (ALL_CAPS) for `from config import ...`.
for _name in dir(_MODULE):
	if _name.isupper():
		globals()[_name] = getattr(_MODULE, _name)

__all__ = [_name for _name in globals() if _name.isupper()]
