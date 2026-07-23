# API reference

Generated from the docstrings in the source, so it cannot fall out of step with the code.

## `pymzlib.pride`

::: pymzlib.pride
    options:
      show_root_full_path: false
      heading_level: 3

## `pymzlib.peptidoform`

::: pymzlib.peptidoform
    options:
      show_root_full_path: false
      heading_level: 3

## `pymzlib.flashlfq`

::: pymzlib.flashlfq
    options:
      show_root_full_path: false
      heading_level: 3

## Errors and diagnostics

Everything pyMzLib raises inherits from `PyMzLibError`, so a single `except` catches all of it.

::: pymzlib._bridge
    options:
      show_root_full_path: false
      heading_level: 3
      members:
        - PyMzLibError
        - UsageError
        - BridgeError
        - ServiceUnavailableError
        - BridgeTimeoutError
        - BridgeNotFoundError
        - bridge_path
        - bridge_version
