"""Input adapters mapping external formats into the internal model."""

from __future__ import annotations

from .hfx import HFXError, HFXOptions, load_schema, read_hfx

__all__ = ["HFXError", "HFXOptions", "load_schema", "read_hfx"]
