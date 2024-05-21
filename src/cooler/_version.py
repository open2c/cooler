try:
    from importlib.metadata import PackageNotFoundError, version
except ImportError:
    from importlib_metadata import PackageNotFoundError, version

try:
    __version__ = version("cooler")
except PackageNotFoundError:
    __version__ = "unknown"

__format_version__ = 3
__format_version_mcool__ = 2
__format_version_scool__ = 1
