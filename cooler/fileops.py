# -*- coding: utf-8 -*-
from __future__ import absolute_import, print_function, division
from datetime import datetime
from textwrap import dedent
import posixpath
import tempfile
import warnings
import uuid
import json
import os

from six.moves import map
from six import PY2
import six
try:
    from json import JSONDecodeError
except ImportError:
    JSONDecodeError = ValueError  # PY35+

from asciitree import BoxStyle, LeftAligned
from asciitree.traversal import Traversal
from pandas.api.types import is_categorical, is_integer
import pandas as pd
import numpy as np
import h5py

from ._logging import get_logger
from .util import parse_cooler_uri, natsorted
from .create import MAGIC, URL

__all__ = [
    'is_cooler', 'is_multires_file', 'list_coolers', 'cp', 'mv', 'ln'
]


class TreeNode(object):

    def __init__(self, obj, depth=0, level=None):
        self.obj = obj
        self.depth = depth
        self.level = level

    def get_type(self):
        return type(self.obj).__name__

    def get_children(self):
        if hasattr(self.obj, 'values'):
            if self.level is None or self.depth < self.level:
                depth = self.depth + 1
                children = self.obj.values()
                return [self.__class__(o, depth=depth, level=self.level)
                            for o in children]
        return []

    def get_text(self):
        name = self.obj.name.split("/")[-1] or "/"
        if hasattr(self.obj, 'shape'):
            name += ' {} {}'.format(self.obj.shape, self.obj.dtype)
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
    fmt = grp.attrs.get('format', None)
    url = grp.attrs.get('format-url', None)
    if fmt == MAGIC or url == URL:
        keys = ('chroms', 'bins', 'pixels', 'indexes')
        if not all(name in grp.keys() for name in keys):
            warnings.warn(
                'Cooler path /{} appears to be corrupt'.format(pth))
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
    filepath, grouppath = parse_cooler_uri(uri)
    if not h5py.is_hdf5(filepath):
        return False

    with h5py.File(filepath) as f:
        fmt = f.attrs.get('format', None)
        if 'resolutions' in f.keys() and len(f['resolutions'].keys()) > 0:
            name = next(list(f['resolutions'].keys()))
            if fmt == 'HDF5::MCOOL' and _is_cooler(f['resolutions'][name]):
                return True
        elif '0' in f.keys() and _is_cooler(f['0']) and min_version < 2:
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
            listing.append('/' + pth if not pth.startswith('/') else pth)

    with h5py.File(filepath, 'r') as f:
        _check_cooler('/', f)
        visititems(f, _check_cooler)

    return natsorted(listing)


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
        listing.append('/' + pth if not pth.startswith('/') else pth)

    with h5py.File(filepath, 'r') as f:
        _check_all(grouppath, f)
        visititems(f[grouppath], _check_all)

    return listing


def _copy(src_uri, dst_uri, overwrite, link, rename, soft_link):
    src_path, src_group = parse_cooler_uri(src_uri)
    dst_path, dst_group = parse_cooler_uri(dst_uri)

    if sum([link, rename, soft_link]) > 1:
        raise ValueError(
            'Must provide at most one of: "link", "rename", "soft_link"')

    if not os.path.isfile(dst_path) or overwrite:
        write_mode = 'w'
    else:
        write_mode = 'r+'

    with h5py.File(src_path, 'r+') as src, \
         h5py.File(dst_path, write_mode) as dst:

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
                if dst_group == '/':
                    for subgrp in src[src_group].keys():
                        src.copy(src_group + '/' + subgrp, dst, subgrp)
                    dst[dst_group].attrs.update(src[src_group].attrs)
                else:
                    src.copy(
                        src_group, dst,
                        dst_group if dst_group != '/' else None)


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


def _tree_html(node, root=False, expand=False):
    result = ''
    data_jstree = '{"type": "%s"}' % node.get_type()
    if root or (expand is True) or (isinstance(expand, int) and node.depth < expand):
        css_class = 'jstree-open'
    else:
        css_class = ''
    result += "<li data-jstree='{}' class='{}'>".format(data_jstree, css_class)
    result += '<span>{}</span>'.format(node.get_text())
    children = node.get_children()
    if children:
        result += '<ul>'
        for child in children:
            result += _tree_html(child, expand=expand)
        result += '</ul>'
    result += '</li>'
    return result


def tree_html(group, expand, level):
    tree_group_icon = 'fa fa-folder'
    tree_array_icon = 'fa fa-table'
    # alternatives...
    # tree_group_icon: 'jstree-folder'
    # tree_array_icon: 'jstree-file'

    result = ''

    # include CSS for jstree default theme
    css_url = '//cdnjs.cloudflare.com/ajax/libs/jstree/3.3.3/themes/default/style.min.css'
    result += '<link rel="stylesheet" href="{}"/>'.format(css_url)

    # construct the tree as HTML nested lists
    node_id = uuid.uuid4()
    result += '<div id="{}" class="zarr-tree">'.format(node_id)
    result += '<ul>'

    root = TreeNode(group, level=level)
    result += _tree_html(root, root=True, expand=expand)

    result += '</ul>'
    result += '</div>'

    # construct javascript
    result += dedent("""
        <script>
            if (!require.defined('jquery')) {
                require.config({
                    paths: {
                        jquery: '//cdnjs.cloudflare.com/ajax/libs/jquery/1.12.1/jquery.min'
                    },
                });
            }
            if (!require.defined('jstree')) {
                require.config({
                    paths: {
                        jstree: '//cdnjs.cloudflare.com/ajax/libs/jstree/3.3.3/jstree.min'
                    },
                });
            }
            require(['jstree'], function() {
                $('#%s').jstree({
                    types: {
                        Group: {
                            icon: "%s"
                        },
                        Array: {
                            icon: "%s"
                        }
                    },
                    plugins: ["types"]
                });
            });
        </script>
    """ % (node_id, tree_group_icon, tree_array_icon))

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

        self.text_kwargs = dict(
            horiz_len=2,
            label_space=1,
            indent=1
        )

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
            VERTICAL_AND_RIGHT=u"\u251C"
        )

        self.node_cls = node_cls

    def __bytes__(self):
        drawer = LeftAligned(
            traverse=TreeTraversal(),
            draw=BoxStyle(gfx=self.bytes_kwargs, **self.text_kwargs)
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
            draw=BoxStyle(gfx=self.unicode_kwargs, **self.text_kwargs)
        )
        root = self.node_cls(self.group, level=self.level)
        return drawer(root)

    def __repr__(self):
        if PY2:  # pragma: py3 no cover
            return self.__bytes__()
        else:  # pragma: py2 no cover
            return self.__unicode__()

    def _repr_html_(self):
        return tree_html(self.group, expand=self.expand, level=self.level)


def pprint_data_tree(uri, level):
    path, group = parse_cooler_uri(uri)
    with h5py.File(path, 'r') as f:
        grp = f[group]
        return repr(TreeViewer(grp, level=level))


def _decode_attr_value(obj):
    if hasattr(obj, 'item'):
        o = np.asscalar(obj)
    elif hasattr(obj, 'tolist'):
        o = obj.tolist()
    elif isinstance(obj, six.string_types):
        try:
            o = datetime.strptime(obj, '%Y-%m-%dT%H:%M:%S.%f')
        except ValueError:
            try:
                o = json.loads(obj)
            except JSONDecodeError:
                o = obj
    else:
        o = obj
    return o


def read_attr_tree(group, level):

    def _getdict(node, root=False):
        attrs = node.obj.attrs
        result = {'@attrs': {k: _decode_attr_value(v)
                                for k, v in attrs.items()}}
        children = node.get_children()
        if children:
            for child in children:
                result[child.get_text()] = _getdict(child)
        return result

    return _getdict(AttrNode(group, level=level), root=True)


def pprint_attr_tree(uri, level):
    import yaml
    from io import StringIO
    path, group = parse_cooler_uri(uri)
    with h5py.File(path, 'r') as f:
        grp = f[group]
        s = StringIO()
        yaml.dump(read_attr_tree(grp, level), s)
        return s.getvalue()



# if not h5py.is_hdf5(filepath):
#     return False
# if group is None:
#     return len(ls(filepath)) > 0
# if not group.startswith('/'):
#     group = '/' + group
# return group in ls(filepath)

# def has_cooler(uri):
#     """
#     Determine if a file contains a valid Cooler data hierarchy.

#     Parameters
#     ----------
#     filepath : str
#     group : str, optional
#         Path to specific group to check. Otherwise returns True
#         if any Cooler paths are detected.

#     Returns
#     -------
#     bool

#     """
#     filepath, grouppath = parse_cooler_uri(uri)
#     if not h5py.is_hdf5(filepath):
#         return False
#     return len(list_coolers(uri)) > 0
