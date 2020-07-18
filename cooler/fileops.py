# -*- coding: utf-8 -*-
from __future__ import absolute_import, print_function, division
from datetime import datetime
# from textwrap import dedent
import warnings
import simplejson as json
import os
from six import PY2
import six

try:
    from simplejson import JSONDecodeError
except ImportError:
    JSONDecodeError = ValueError  # PY35+

from asciitree import BoxStyle, LeftAligned
from asciitree.traversal import Traversal
import numpy as np
import h5py

from .util import parse_cooler_uri, natsorted
from .create import MAGIC, URL, MAGIC_SCOOL

__all__ = ["is_cooler", "is_multires_file", "list_coolers", "cp", "mv", "ln"]


def json_dumps(o):
    """Write JSON in a consistent, human-readable way."""
    return json.dumps(
        o, indent=4, sort_keys=True, ensure_ascii=True, separators=(',', ': ')
    )


def json_loads(s):
    """Read JSON in a consistent way."""
    return json.loads(s)


def decode_attr_value(obj):
    """Decode a numpy object or HDF5 attribute into a JSONable object"""
    if hasattr(obj, "item"):
        o = obj.item()
    elif hasattr(obj, "tolist"):
        o = obj.tolist()
    elif isinstance(obj, six.string_types):
        try:
            o = datetime.strptime(obj, "%Y-%m-%dT%H:%M:%S.%f")
        except ValueError:
            try:
                o = json.loads(obj)
            except JSONDecodeError:
                o = obj
    else:
        o = obj
    return o


class TreeNode(object):
    def __init__(self, obj, depth=0, level=None):
        self.obj = obj
        self.depth = depth
        self.level = level

    def get_type(self):
        return type(self.obj).__name__

    def get_children(self):
        if hasattr(self.obj, "values"):
            if self.level is None or self.depth < self.level:
                depth = self.depth + 1
                children = self.obj.values()
                return [
                    self.__class__(o, depth=depth, level=self.level) for o in children
                ]
        return []

    def get_text(self):
        name = self.obj.name.split("/")[-1] or "/"
        if hasattr(self.obj, "shape"):
            name += " {} {}".format(self.obj.shape, self.obj.dtype)
        return name


class AttrNode(TreeNode):
    def get_text(self):
        return self.obj.name.split("/")[-1] or "/"


def visititems(group, func, level=None):
    """Like :py:method:`h5py.Group.visititems`, but much faster somehow.
    """

    def _visititems(node, func, result=None):
        children = node.get_children()
        if children:
            for child in children:
                result[child.obj.name] = func(child.obj.name, child.obj)
                _visititems(child, func, result)
        return result

    root = TreeNode(group, level=level)
    return _visititems(root, func, {})


def _is_cooler(grp):
    fmt = grp.attrs.get("format", None)
    if fmt == MAGIC:
        keys = ("chroms", "bins", "pixels", "indexes")
        if not all(name in grp.keys() for name in keys):
            warnings.warn("Cooler path {} appears to be corrupt".format(grp.name))
        return True
    return False


def is_cooler(uri):
    """
    Determine if a URI string references a cooler data collection.
    Returns False if the file or group path doesn't exist.
    
    """
    filepath, grouppath = parse_cooler_uri(uri)
    if not h5py.is_hdf5(filepath):
        return False
    with h5py.File(filepath) as f:
        return _is_cooler(f[grouppath])


def is_multires_file(filepath, min_version=1):
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


def is_scool_file(filepath):
    """
    Determine if a file is a single-cell cooler file.
    Returns False if the file doesn't exist.

    """
    if not h5py.is_hdf5(filepath):
        raise OSError("'{}' is not an HDF5 file.".format(filepath))
        return False

    with h5py.File(filepath) as f:
        fmt = f.attrs.get("format", None)
        if fmt == MAGIC_SCOOL:
            keys = ("chroms", "bins", "cells")
            if not all(name in f.keys() for name in keys):
                warnings.warn("Scooler path {} appears to be corrupt".format(grp.name))
                return False
            if "cells" in f.keys() and len(f["cells"].keys()) > 0:
                for cells in f["cells"].keys():
                    if not _is_cooler(f["cells"][cells]):
                        return False
                return True
    return False


def list_coolers(filepath):
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
        raise OSError("'{}' is not an HDF5 file.".format(filepath))

    listing = []

    def _check_cooler(pth, grp):
        if _is_cooler(grp):
            listing.append("/" + pth if not pth.startswith("/") else pth)

    with h5py.File(filepath, "r") as f:
        _check_cooler("/", f)
        visititems(f, _check_cooler)

    return natsorted(listing)


def list_scool_cells(filepath):
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
    if is_scool_file(filepath):

        listing = []
        def _check_cooler(pth, grp):
            if _is_cooler(grp):
                listing.append("/" + pth if not pth.startswith("/") else pth)
        with h5py.File(filepath, "r") as f:
            _check_cooler("/", f)
            visititems(f, _check_cooler)
        if '/' in listing:
            listing.remove('/')
        return natsorted(listing)
    else:
        raise OSError("'{}' is not a scool file.".format(filepath))
    return False


def ls(uri):
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
        raise OSError("'{}' is not an HDF5 file.".format(filepath))

    listing = []

    def _check_all(pth, grp):
        listing.append("/" + pth if not pth.startswith("/") else pth)

    with h5py.File(filepath, "r") as f:
        _check_all(grouppath, f)
        visititems(f[grouppath], _check_all)

    return listing


def _copy(src_uri, dst_uri, overwrite, link, rename, soft_link):
    src_path, src_group = parse_cooler_uri(src_uri)
    dst_path, dst_group = parse_cooler_uri(dst_uri)

    if sum([link, rename, soft_link]) > 1:
        raise ValueError('Must provide at most one of: "link", "rename", "soft_link"')

    if not os.path.isfile(dst_path) or overwrite:
        write_mode = "w"
    else:
        write_mode = "r+"

    with h5py.File(src_path, "r+") as src, h5py.File(
        dst_path, write_mode
    ) as dst:  # noqa

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


def cp(src_uri, dst_uri, overwrite=False):
    """Copy a group or dataset from one file to another or within the same file.
    """
    _copy(src_uri, dst_uri, overwrite, link=False, rename=False, soft_link=False)


def mv(src_uri, dst_uri, overwrite=False):
    """Rename a group or dataset within the same file.
    """
    _copy(src_uri, dst_uri, overwrite, link=False, rename=True, soft_link=False)


def ln(src_uri, dst_uri, soft=False, overwrite=False):
    """Create a hard link to a group or dataset in the same file. Also
    supports soft links (in the same file) or external links (different files).
    """
    _copy(src_uri, dst_uri, overwrite, link=not soft, rename=False, soft_link=soft)


######
# Tree rendering. Borrowed from zarr-python.


def _tree_get_icon(stype):
    if stype in {"Dataset", "Array"}:
        return 'table'
    elif stype in {"Group", "File"}:
        return 'folder'
    else:
        raise ValueError("Unknown type: %s" % stype)


def _tree_widget_sublist(node, root=False, expand=False):
    import ipytree

    result = ipytree.Node()
    result.icon = _tree_get_icon(node.get_type())
    if root or (expand is True) or (isinstance(expand, int) and node.depth < expand):
        result.opened = True
    else:
        result.opened = False
    result.name = node.get_text()
    result.nodes = [_tree_widget_sublist(c, expand=expand) for c in node.get_children()]
    result.disabled = True

    return result


def tree_widget(group, expand, level):
    try:
        import ipytree
    except ImportError as error:
        raise ImportError(
            "{}: Run `pip install ipytree` or `conda install ipytree`"
            "to get the required ipytree dependency for displaying the tree "
            "widget. If using jupyterlab, you also need to run "
            "`jupyter labextension install ipytree`".format(error)
        )

    result = ipytree.Tree()
    root = TreeNode(group, level=level)
    result.add_node(_tree_widget_sublist(root, root=True, expand=expand))

    return result


class TreeTraversal(Traversal):
    def get_children(self, node):
        return node.get_children()

    def get_root(self, tree):
        return tree

    def get_text(self, node):
        return node.get_text()


class TreeViewer(object):
    """
    Generates ascii- or html-based reprs for "Groupy" objects.
    Borrowed with minor modifications from the zarr project
    (Zarr Developers, MIT-licensed).

    <https://github.com/zarr-developers/zarr>

    See: zarr.util.TreeViewer, zarr.util.tree_html

    """

    def __init__(self, group, expand=False, level=None, node_cls=TreeNode):
        self.group = group
        self.expand = expand
        self.level = level

        self.text_kwargs = dict(horiz_len=2, label_space=1, indent=1)

        self.bytes_kwargs = dict(
            UP_AND_RIGHT="+",
            HORIZONTAL="-",
            VERTICAL="|",
            VERTICAL_AND_RIGHT="+"
        )

        self.unicode_kwargs = dict(
            UP_AND_RIGHT=u"\u2514",
            HORIZONTAL=u"\u2500",
            VERTICAL=u"\u2502",
            VERTICAL_AND_RIGHT=u"\u251C",
        )

        self.node_cls = node_cls

    def __bytes__(self):
        drawer = LeftAligned(
            traverse=TreeTraversal(),
            draw=BoxStyle(gfx=self.bytes_kwargs, **self.text_kwargs),
        )
        root = self.node_cls(self.group, level=self.level)
        result = drawer(root)

        # Unicode characters slip in on Python 3.
        # So we need to straighten that out first.
        if not PY2:
            result = result.encode()

        return result

    def __unicode__(self):
        drawer = LeftAligned(
            traverse=TreeTraversal(),
            draw=BoxStyle(gfx=self.unicode_kwargs, **self.text_kwargs),
        )
        root = self.node_cls(self.group, level=self.level)
        return drawer(root)

    def __repr__(self):
        if PY2:  # pragma: py3 no cover
            return self.__bytes__()
        else:  # pragma: py2 no cover
            return self.__unicode__()

    def _ipython_display_(self):
        tree = tree_widget(self.group, expand=self.expand, level=self.level)
        tree._ipython_display_()
        return tree


def read_attr_tree(group, level):
    def _getdict(node, root=False):
        attrs = node.obj.attrs
        result = {"@attrs": {k: decode_attr_value(v) for k, v in attrs.items()}}
        children = node.get_children()
        if children:
            for child in children:
                result[child.get_text()] = _getdict(child)
        return result

    return _getdict(AttrNode(group, level=level), root=True)


def pprint_attr_tree(uri, level):
    import yaml
    from six import StringIO

    path, group = parse_cooler_uri(uri)
    with h5py.File(path, "r") as f:
        grp = f[group]
        s = StringIO()
        yaml.dump(read_attr_tree(grp, level), s)
        return s.getvalue()


def pprint_data_tree(uri, level):
    path, group = parse_cooler_uri(uri)
    with h5py.File(path, "r") as f:
        grp = f[group]
        return repr(TreeViewer(grp, level=level))

######
