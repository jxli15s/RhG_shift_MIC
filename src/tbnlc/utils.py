from __future__ import annotations

from dataclasses import asdict, is_dataclass
from typing import Any, Mapping


def get_opt(opts: Any, name: str, default: Any) -> Any:
    """MATLAB-style option getter.

    Supports dict-like object or dataclass.
    """
    if opts is None:
        return default
    if is_dataclass(opts):
        d = asdict(opts)
        v = d.get(name, default)
        return default if v is None else v
    if isinstance(opts, Mapping):
        v = opts.get(name, default)
        return default if v is None else v
    if hasattr(opts, name):
        v = getattr(opts, name, default)
        return default if v is None else v
    return default
