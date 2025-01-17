from __future__ import annotations

# from textwrap import dedent
import json
import os
import warnings
from collections.abc import Callable
from datetime import datetime
from numbers import Number
from typing import Any, Literal

import h5py
import numpy as np
from asciitree import BoxStyle, LeftAligned
from asciitree.traversal import Traversal

from .create import MAGIC, MAGIC_SCOOL
from .util import natsorted, parse_cooler_uri

__all__ = ["cp", "is_cooler", "is_multires_file", "list_coolers", "ln", "mv"]


def json_dumps(o: object) -> str:
    """Write JSON in a consistent, human-readable way."""
    return json.dumps(
        o, indent=4, sort_keys=True, ensure_ascii=True, separators=(",", ": ")
    )


def json_loads(s: str) -> object:
    """Read JSON in a consistent way."""
    return json.loads(s)


def decode_attr_value(
    obj: str | np.generic | np.ndarray | h5py.Dataset
) -> str | Number | datetime | list | dict:
    """
    Decode a HDF5 attribute value (numpy object or string) into something
    JSON serializable.
    """
    if hasattr(obj, "item"):
        o = obj.item()
    elif hasattr(obj, "tolist"):
        o = obj.tolist()
    elif isinstance(obj, str):
        try:
            o = datetime.strptime(obj, "%Y-%m-%dT%H:%M:%S.%f")
        except ValueError:
            try:
                o = json.loads(obj)
            except ValueError:
                o = obj
    else:
        o = obj
    return o


class TreeNode:
    def __init__(self, obj: Any, depth: int = 0, level: int | None = None):
        self.obj = obj
        self.depth = depth
        self.level = level

    def get_type(self) -> str:
        return type(self.obj).__name__

    def get_children(self) -> list[TreeNode]:
        if hasattr(self.obj, "values"):
            if self.level is None or self.depth < self.level:
                depth = self.depth + 1
                children = self.obj.values()
                return [
                    self.__class__(o, depth=depth, level=self.level) for o in children
                ]
        return []

    def get_text(self) -> str:
        name = self.obj.name.split("/")[-1] or "/"
        if hasattr(self.obj, "shape"):
            name += f" {self.obj.shape} {self.obj.dtype}"
        return name


class AttrNode(TreeNode):
    def get_text(self) -> str:
        return self.obj.name.split("/")[-1] or "/"


def visititems(
    group: h5py.Group,
    func: Callable[[str, h5py.Group | h5py.Dataset], Any],
    level: int | None = None
) -> dict[str, Any]:
    """Like :py:method:`h5py.Group.visititems`, but much faster somehow."""

    def _visititems(node, func, result=None):
        children = node.get_children()
        if children:
            for child in children:
                result[child.obj.name] = func(child.obj.name, child.obj)
                _visititems(child, func, result)
        return result

    root = TreeNode(group, level=level)
    return _visititems(root, func, {})


def _is_cooler(grp: h5py.Group) -> bool:
    fmt = grp.attrs.get("format", None)
    if fmt == MAGIC:
        keys = ("chroms", "bins", "pixels", "indexes")
        if not all(name in grp.keys() for name in keys):
            warnings.warn(
                f"Cooler path {grp.name} appears to be corrupt",
                stacklevel=2,
            )
        return True
    return False


def is_cooler(uri: str) -> bool:
    """
    Determine if a URI string references a cooler data collection.
    Returns False if the file or group path doesn't exist.

    """
    filepath, grouppath = parse_cooler_uri(uri)
    if not h5py.is_hdf5(filepath):
        return False
    with h5py.File(filepath) as f:
        return _is_cooler(f[grouppath])


def is_multires_file(filepath: str, min_version: int = 1) -> bool:
    """
    Determine if a file is a multi-res cooler file.
    Returns False if the file doesn't exist.

    """
    if not h5py.is_hdf5(filepath):
        return False

    with h5py.File(filepath) as f:
        fmt = f.attrs.get("format", None)
        if "resolutions" in f.keys() and len(f["resolutions"].keys()) > 0:
            name = next(iter(f["resolutions"].keys()))
            if fmt == "HDF5::MCOOL" and _is_cooler(f["resolutions"][name]):
                return True
        elif "0" in f.keys() and _is_cooler(f["0"]) and min_version < 2:
            return True

    return False


def is_scool_file(filepath: str) -> bool:
    """
    Determine if a file is a single-cell cooler file.
    Returns False if the file doesn't exist.

    """
    if not h5py.is_hdf5(filepath):
        raise OSError(f"'{filepath}' is not an HDF5 file.")
        return False

    with h5py.File(filepath) as f:
        fmt = f.attrs.get("format", None)
        if fmt == MAGIC_SCOOL:
            keys = ("chroms", "bins", "cells")
            if not all(name in f.keys() for name in keys):
                warnings.warn(
                    "Scool file appears to be corrupt",
                    stacklevel=2,
                )
                return False
            if "cells" in f.keys() and len(f["cells"].keys()) > 0:
                for cells in f["cells"].keys():
                    if not _is_cooler(f["cells"][cells]):
                        return False
                return True
    return False


def list_coolers(filepath: str) -> list[str]:
    """
    List group paths to all cooler data collections in a file.

    Parameters
    ----------
    filepath : str

    Returns
    -------
    list
        Cooler group paths in the file.

    """
    if not h5py.is_hdf5(filepath):
        raise OSError(f"'{filepath}' is not an HDF5 file.")

    listing = []

    def _check_cooler(pth, grp):
        if _is_cooler(grp):
            listing.append("/" + pth if not pth.startswith("/") else pth)

    with h5py.File(filepath, "r") as f:
        _check_cooler("/", f)
        visititems(f, _check_cooler)

    return natsorted(listing)


def list_scool_cells(filepath: str) -> list[str]:
    """
    List the paths to all single-cell cool matrices in a file scool file.

    Parameters
    ----------
    filepath : str

    Returns
    -------
    list
        Cooler group paths of all cells in the file.

    """

    def _check_cooler(pth, grp):
        if _is_cooler(grp):
            listing.append("/" + pth if not pth.startswith("/") else pth)

    if is_scool_file(filepath):
        listing = []
        with h5py.File(filepath, "r") as f:
            _check_cooler("/", f)
            visititems(f, _check_cooler)
        if "/" in listing:
            listing.remove("/")
        return natsorted(listing)
    else:
        raise OSError(f"'{filepath}' is not a scool file.")


def ls(uri: str) -> list[str]:
    """
    Get all groups and datasets in an HDF5 file.

    Parameters
    ----------
    uri : str

    Returns
    -------
    list
        Group and dataset paths.

    """
    filepath, grouppath = parse_cooler_uri(uri)
    if not h5py.is_hdf5(filepath):
        raise OSError(f"'{filepath}' is not an HDF5 file.")

    listing = []

    def _check_all(pth, grp):
        listing.append("/" + pth if not pth.startswith("/") else pth)

    with h5py.File(filepath, "r") as f:
        _check_all(grouppath, f)
        visititems(f[grouppath], _check_all)

    return listing


def _copy(
    src_uri: str,
    dst_uri: str,
    overwrite: bool,
    link: bool,
    rename:bool,
    soft_link: bool
) -> None:
    src_path, src_group = parse_cooler_uri(src_uri)
    dst_path, dst_group = parse_cooler_uri(dst_uri)

    if sum([link, rename, soft_link]) > 1:
        raise ValueError('Must provide at most one of: "link", "rename", "soft_link"')

    if not os.path.isfile(dst_path) or overwrite:
        dst_write_mode = "w"
    else:
        dst_write_mode = "r+"

    if src_path == dst_path:
        src_write_mode = "r+"
    else:
        src_write_mode = "r"

    with h5py.File(src_path, src_write_mode) as src, \
         h5py.File(dst_path, dst_write_mode) as dst:  # fmt: skip
        if src_path == dst_path:
            if link or rename:
                src[dst_group] = src[src_group]
                if rename:
                    del src[src_group]
            elif soft_link:
                src[dst_group] = h5py.SoftLink(src_group)
            else:
                src.copy(src_group, dst_group)
        else:
            if link:
                raise OSError("Can't hard link between two different files.")
            elif soft_link:
                dst[dst_group] = h5py.ExternalLink(src_path, src_group)
            else:
                if dst_group == "/":
                    for subgrp in src[src_group].keys():
                        src.copy(src_group + "/" + subgrp, dst, subgrp)
                    dst[dst_group].attrs.update(src[src_group].attrs)
                else:
                    src.copy(src_group, dst, dst_group if dst_group != "/" else None)


def cp(src_uri: str, dst_uri: str, overwrite: bool = False) -> None:
    """Copy a group or dataset from one file to another or within the same file."""
    _copy(src_uri, dst_uri, overwrite, link=False, rename=False, soft_link=False)


def mv(src_uri: str, dst_uri: str, overwrite: bool = False) -> None:
    """Rename a group or dataset within the same file."""
    _copy(src_uri, dst_uri, overwrite, link=False, rename=True, soft_link=False)


def ln(
    src_uri: str, dst_uri: str, soft: bool = False, overwrite: bool = False
) -> None:
    """Create a hard link to a group or dataset in the same file. Also
    supports soft links (in the same file) or external links (different files).
    """
    _copy(src_uri, dst_uri, overwrite, link=not soft, rename=False, soft_link=soft)


######
# Tree rendering. Borrowed from zarr-python.

def tree_widget(group: h5py.Group, expand: bool, level: int | None) -> Any:
    try:
        import ipytree
    except ImportError as error:
        raise ImportError(
            f"{error}: Run `pip install ipytree` or `conda install ipytree`"
            "to get the required ipytree dependency for displaying the tree "
            "widget. If using jupyterlab, you also need to run "
            "`jupyter labextension install ipytree`"
        ) from None

    def _tree_get_icon(
        stype: Literal["Dataset", "Array", "Group", "File"]
    ) -> Literal["table", "folder"]:
        if stype in {"Dataset", "Array"}:
            return "table"
        elif stype in {"Group", "File"}:
            return "folder"
        else:
            raise ValueError(f"Unknown type: {stype}")

    def _tree_widget_sublist(
        node: TreeNode, root: bool = False, expand: bool = False
    ) -> ipytree.Node:

        result = ipytree.Node()
        result.icon = _tree_get_icon(node.get_type())
        if (
            root or (expand is True) or
            (isinstance(expand, int) and node.depth < expand)
        ):
            result.opened = True
        else:
            result.opened = False
        result.name = node.get_text()
        result.nodes = [
            _tree_widget_sublist(c, expand=expand) for c in node.get_children()
        ]
        result.disabled = True

        return result

    result = ipytree.Tree()
    root = TreeNode(group, level=level)
    result.add_node(_tree_widget_sublist(root, root=True, expand=expand))

    return result


class TreeTraversal(Traversal):
    def get_children(self, node: TreeNode) -> list[TreeNode]:
        return node.get_children()

    def get_root(self, tree: TreeNode) -> TreeNode:
        return tree

    def get_text(self, node: TreeNode) -> str:
        return node.get_text()


class TreeViewer:
    """
    Generates ascii- or html-based reprs for "Groupy" objects.
    Borrowed with minor modifications from the zarr project
    (Zarr Developers, MIT-licensed).

    <https://github.com/zarr-developers/zarr>

    See: zarr.util.TreeViewer, zarr.util.tree_html

    """

    def __init__(
        self,
        group: h5py.Group,
        expand: bool = False,
        level: int | None = None,
        node_cls: Any = TreeNode
    ):
        self.group = group
        self.expand = expand
        self.level = level

        self.text_kwargs = {"horiz_len": 2, "label_space": 1, "indent": 1}

        self.bytes_kwargs = {
            "UP_AND_RIGHT": "+",
            "HORIZONTAL": "-",
            "VERTICAL": "|",
            "VERTICAL_AND_RIGHT": "+",
        }

        self.unicode_kwargs = {
            "UP_AND_RIGHT": "\u2514",
            "HORIZONTAL": "\u2500",
            "VERTICAL": "\u2502",
            "VERTICAL_AND_RIGHT": "\u251C",
        }

        self.node_cls = node_cls

    def __bytes__(self) -> bytes:
        drawer = LeftAligned(
            traverse=TreeTraversal(),
            draw=BoxStyle(gfx=self.bytes_kwargs, **self.text_kwargs),
        )
        root = self.node_cls(self.group, level=self.level)
        result = drawer(root)

        # Unicode characters slip in on Python 3.
        # So we need to straighten that out first.
        result = result.encode()

        return result

    def __unicode__(self) -> str:
        drawer = LeftAligned(
            traverse=TreeTraversal(),
            draw=BoxStyle(gfx=self.unicode_kwargs, **self.text_kwargs),
        )
        root = self.node_cls(self.group, level=self.level)
        return drawer(root)

    def __repr__(self) -> str:
        return self.__unicode__()

    def _repr_mimebundle_(self) -> Any:
        tree = tree_widget(self.group, expand=self.expand, level=self.level)
        tree._repr_mimebundle_()
        return tree


def read_attr_tree(group: h5py.Group, level: int | None) -> dict[str, Any]:
    def _getdict(node, root=False):
        attrs = node.obj.attrs
        result = {"@attrs": {k: decode_attr_value(v) for k, v in attrs.items()}}
        children = node.get_children()
        if children:
            for child in children:
                result[child.get_text()] = _getdict(child)
        return result

    return _getdict(AttrNode(group, level=level), root=True)


def pprint_attr_tree(uri: str, level: int | None) -> str:
    from io import StringIO

    import yaml

    path, group = parse_cooler_uri(uri)
    with h5py.File(path, "r") as f:
        grp = f[group]
        s = StringIO()
        yaml.dump(read_attr_tree(grp, level), s)
        return s.getvalue()


def pprint_data_tree(uri: str, level: int | None) -> str:
    path, group = parse_cooler_uri(uri)
    with h5py.File(path, "r") as f:
        grp = f[group]
        return repr(TreeViewer(grp, level=level))


######
