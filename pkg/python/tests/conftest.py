"""Shared test support.

The important thing here is :func:`external_service` — the Python counterpart of mzLib's
``ExternalServiceTestHelper``. Tests that touch a live service must be able to tell two failures
apart:

* the service is unavailable (down, rate-limited, 5xx, timed out) — **not** our bug, so the test
  should be **skipped** with a message saying so; versus
* the service answered but the contract is broken (wrong URL, response no longer parses, an
  expected value missing) — a real regression that must **fail**.

Without that distinction a red build is ambiguous, and an ambiguous red build gets ignored, which
is how a genuine contract break goes unnoticed for a month.

The classification itself lives in the bridge, not here: it raises ``ServiceUnavailableError``
for availability failures, so every consumer of the wire format benefits, not just this suite.
"""

from __future__ import annotations

import contextlib

import pytest

import pymzlib


@contextlib.contextmanager
def external_service(service_name: str = "PRIDE Archive"):
    """Skip rather than fail when an external service is unavailable.

    Args:
        service_name: Named in the skip message, so the CI log says which service was down.

    Example:
        >>> with external_service():                       # doctest: +SKIP
        ...     files = pymzlib.pride.list_files("PXD000001")
        ...     assert files
    """
    try:
        yield
    except pymzlib.ServiceUnavailableError as exc:
        pytest.skip(
            f"Skipping external-service test: {service_name} unavailable ({exc}). "
            "This is a third-party availability problem, not a code failure."
        )
