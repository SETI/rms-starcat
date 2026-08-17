################################################################################
# tests/test_package.py
################################################################################

from __future__ import annotations

import builtins
import importlib
import sys
from typing import Any

import pytest

import starcat


def test_public_names() -> None:
    for name in starcat.__all__:
        assert hasattr(starcat, name), name

    assert 'Star' in starcat.__all__
    assert 'StarCatalog' in starcat.__all__
    assert 'SpiceStar' in starcat.__all__
    assert 'SpiceStarCatalog' in starcat.__all__
    assert 'UCAC4Star' in starcat.__all__
    assert 'UCAC4StarCatalog' in starcat.__all__
    assert 'YBSCStar' in starcat.__all__
    assert 'YBSCStarCatalog' in starcat.__all__


def test_version() -> None:
    assert isinstance(starcat.__version__, str)
    assert starcat.__version__ != ''


def test_version_without_version_module(monkeypatch: pytest.MonkeyPatch) -> None:
    # _version.py is generated at build time, so it is missing when the package
    # is used directly from a source checkout
    real_import = builtins.__import__

    def blocked_import(name: str, *args: Any, **kwargs: Any) -> Any:
        if name in ('_version', 'starcat._version'):
            raise ImportError(name)
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, '__import__', blocked_import)
    monkeypatch.delitem(sys.modules, 'starcat._version', raising=False)

    try:
        importlib.reload(starcat)
        assert starcat.__version__ == 'Version unspecified'
    finally:
        monkeypatch.undo()
        importlib.reload(starcat)

    assert isinstance(starcat.__version__, str)
