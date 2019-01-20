import warnings

def ls(*args, **kwargs):
    warnings.warn(
        "`cooler.io` is deprecated in 0.8, will be removed in 0.9. "
        "Use `cooler.fileops.list_coolers()` instead.",
        category=FutureWarning, stacklevel=2)
    from .fileops import list_coolers
    return list_coolers(*args, **kwargs)


def is_cooler(*args, **kwargs):
    warnings.warn(
        "`cooler.io` is deprecated in 0.8, will be removed in 0.9. "
        "Use `cooler.fileops.is_cooler()` instead.",
        category=FutureWarning, stacklevel=2)
    from .fileops import is_cooler
    return is_cooler(*args, **kwargs)


def create(*args, **kwargs):
    warnings.warn(
        "`cooler.io.create()` is deprecated in 0.8, will be removed in 0.9. "
        "Use `cooler.create_cooler()` with ordered=True instead.",
        category=FutureWarning, stacklevel=2)
    from .create import create
    return create(*args, **kwargs)


def create_from_unordered(*args, **kwargs):
    warnings.warn(
        "`cooler.io.create_from_unordered()` is deprecated in 0.8, "
        "will be removed in 0.9. Use `cooler.create_cooler()` instead.",
        category=FutureWarning, stacklevel=2)
    from .create import create_from_unordered
    return create_from_unordered(*args, **kwargs)
