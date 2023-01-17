import sys
from warnings import warn

from . import parallel

deprecated_names = ["partition", "split", "lock", "MultiplexDataPipe"]


if sys.version_info[0] == 3 and sys.version_info[1] >= 7:

    def __getattr__(name):
        if name in deprecated_names:
            warn(
                "The `cooler.tools` module is deprecated in v0.9 and will be "
                "removed in v0.10. Use `cooler.parallel` instead.",
                category=FutureWarning,
                stacklevel=2,
            )
            return getattr(parallel, name)
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

else:
    from .parallel import *  # noqa
