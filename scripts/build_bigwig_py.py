"""Minimal pure-Python bigWig writer (standard library only).

This module exists because compiling upstream libBigWig/pyBigWig with a
MinGW toolchain on Windows produces a writer that segfaults (access
violation inside ``bwAddIntervals`` / at ``close()``) even for trivial
inputs, regardless of zoom levels or any ``io.c`` patch.  The UCSC
``bedGraphToBigWig`` Windows binaries are also no longer published
(``admin/exe/windows.x86_64/`` returns HTTP 404).  This writer ports the
on-disk layout of libBigWig's ``bwWrite.c`` (pyBigWig 0.3.26 sdist) to
Python so the maintainer can regenerate bigWig files on any platform
without a C toolchain.

On-disk layout reproduced here (all values little-endian):

* fixed header (64 bytes): magic ``0x888FFC26``, version 4, zero zoom
  levels; chromosome-tree/data/index offsets are patched in place after
  every section has been written (the same multi-pass approach libBigWig
  takes with ``writeAtPos``);
* chromosome list as a B+ tree (magic ``0x78CA8C91``); keys are padded
  to ``keySize`` bytes like libBigWig's ``writeChromList``;
* data section: 8-byte block-count stub followed by zlib-compressed
  bedGraph sections (24-byte section header + 12 bytes per item);
* R-tree index over the data sections (magic ``0x2468ACE0``) with 64
  children per node, exactly like libBigWig's ``writeIndex``;
* total summary (40 bytes: nBasesCovered/minVal/maxVal/sumData/
  sumSquared) written right after the fixed header;
* trailing bigWig magic, as libBigWig's ``bwFinalize`` emits.

Produced files are plain bigWig files without zoom levels and with
``fieldCount = definedFieldCount = 0`` (the values libBigWig itself
writes); they are readable by pyBigWig, UCSC tools and IGV.

Typical use::

    from build_bigwig_py import write_bigwig

    write_bigwig(
        "out.bw",
        [("chr5", 151834684)],
        [("chr5", 1000, 1025, 0.5), ...],  # half-open, 0-based
    )

CLI form::

    python build_bigwig_py.py input.bedGraph chrom.sizes output.bw
"""

from __future__ import annotations

import argparse
import os
import sys
import zlib
from struct import pack

BIGWIG_MAGIC = 0x888FFC26
CHROM_TREE_MAGIC = 0x78CA8C91  # libBigWig name: CIRTREE_MAGIC
RTREE_INDEX_MAGIC = 0x2468ACE0  # libBigWig name: IDX_MAGIC

#: Raw (uncompressed) section budget.  Readers allocate their inflate
#: buffer from the header's ``uncompressBufSize`` field, which this
#: writer pins to the same value libBigWig uses in ``bwCreateHdr``.
UNCOMPRESSED_BLOCK_SIZE = 32768

#: Maximum children per R-tree node (libBigWig ``writeBuffer.blockSize``).
RTREE_BLOCK_SIZE = 64

_HEADER_SIZE = 64
_SUMMARY_SIZE = 40
_INDEX_HEADER_SIZE = 48
_SECTION_HEADER_SIZE = 24
_BEDGRAPH_ITEM_SIZE = 12
_MAX_ITEMS_PER_SECTION = (
    UNCOMPRESSED_BLOCK_SIZE - _SECTION_HEADER_SIZE
) // _BEDGRAPH_ITEM_SIZE


def write_bigwig(path, chrom_sizes, entries, compresslevel=9):
    """Write ``entries`` to ``path`` as a (zoom-free) bigWig file.

    Parameters
    ----------
    path:
        Output file path.
    chrom_sizes:
        Sequence of ``(chrom, size)`` tuples; chromosome ids are
        assigned in this order and every chromosome present in
        ``entries`` must be listed here.
    entries:
        Sequence of ``(chrom, start, end, value)`` tuples with
        half-open, 0-based coordinates.  Entries are sorted internally
        (stable) by chromosome id and start coordinate; ``None`` values
        are skipped, leaving gaps in the bedGraph coverage.
    compresslevel:
        zlib compression level for data blocks (default 9).

    Returns
    -------
    dict
        Summary with entry/block counts and the output file size.
    """
    chrom_sizes = list(chrom_sizes)
    tid_of = {name: index for index, (name, _size) in enumerate(chrom_sizes)}
    items = []
    for chrom, start, end, value in entries:
        if chrom not in tid_of:
            raise ValueError(f"chromosome {chrom!r} is not listed in chrom_sizes")
        start = int(start)
        end = int(end)
        if start < 0 or end <= start:
            raise ValueError(f"invalid interval {chrom}:{start}-{end}")
        if value is not None:
            items.append((tid_of[chrom], start, end, float(value)))
    items.sort(key=lambda item: (item[0], item[1], item[2]))

    header_bytes = _pack_header()
    chrom_tree_bytes = _pack_chrom_tree(chrom_sizes)
    summary_bytes = _pack_summary(items)
    data_bytes, blocks = _pack_data_sections(items, compresslevel)
    chrom_tree_start = _HEADER_SIZE + _SUMMARY_SIZE
    data_start = chrom_tree_start + len(chrom_tree_bytes)
    # The R-tree stores ABSOLUTE file offsets of every data section; the
    # sections were laid out relative to the block-count stub, so shift
    # them into file coordinates (stub = 8 bytes before the first block).
    for block in blocks:
        block["offset"] += data_start + 8
    index_start = data_start + 8 + len(data_bytes)
    index_bytes = _pack_rtree(blocks, index_start)
    with open(path, "wb") as fh:
        fh.write(header_bytes)
        fh.write(summary_bytes)
        fh.write(chrom_tree_bytes)
        fh.write(pack("<Q", 0))  # data-block count stub; patched below
        fh.write(data_bytes)
        fh.write(index_bytes)
        fh.write(pack("<I", BIGWIG_MAGIC))  # trailing magic
        fh.flush()

        # ---- patch header offsets now that every size is known ----
        _patch_u16(fh, 0x06, 0)  # zoom levels
        _patch_u64(fh, 0x08, chrom_tree_start)  # chromTreeOffset
        _patch_u64(fh, 0x10, data_start)  # dataOffset
        _patch_u64(fh, 0x18, index_start)  # fullIndexOffset
        _patch_u64(fh, 0x2C, _HEADER_SIZE)  # totalSummaryOffset
        max_raw = max((block["raw_size"] for block in blocks), default=0)
        # uncompressBufSize: >= every raw section (libBigWig pins 32768)
        _patch_u32(
            fh,
            0x34,
            max(UNCOMPRESSED_BLOCK_SIZE if blocks else 0, max_raw),
        )
        _patch_u64(fh, data_start, len(blocks))  # number of data blocks

    return {
        "n_entries": len(items),
        "n_blocks": len(blocks),
        "n_chroms": len(chrom_sizes),
        "bytes": os.path.getsize(path),
    }


# ---------------------------------------------------------------------------
# fixed header + chromosome B+ tree ("circles" tree)
# ---------------------------------------------------------------------------

def _pack_header():
    header = bytearray(_HEADER_SIZE)
    header[0:4] = pack("<I", BIGWIG_MAGIC)
    header[4:6] = pack("<H", 4)  # version
    header[6:8] = pack("<H", 0)  # zoom levels (no zoom headers follow)
    header[32:34] = pack("<H", 0)  # fieldCount (int, like libBigWig)
    header[34:36] = pack("<H", 0)  # definedFieldCount
    # 0x08 chromTreeOffset, 0x10 dataOffset, 0x18 fullIndexOffset,
    # 0x2c totalSummaryOffset, 0x34 uncompressBufSize patched in place.
    return bytes(header)


def _pack_chrom_tree(chrom_sizes):
    names = [name for name, _size in chrom_sizes]
    key_size = max(len(name) for name in names)
    n_keys = len(names)
    # libBigWig rule: items per leaf block is min(nKeys, 0x7FFF).
    n_per_block = min(n_keys, 0x7FFF)
    n_blocks = (n_keys + n_per_block - 1) // n_per_block

    tree = bytearray()
    tree += pack("<IIII", CHROM_TREE_MAGIC, n_per_block, key_size, 8)
    tree += pack("<QQ", n_keys, 0)  # itemCount + 8 reserved bytes
    if n_blocks > 1:
        # Internal node listing the first key of every leaf block.
        leaf_size = n_per_block * (key_size + 8) + 4
        # nonLeafEnd = end of this node's entry array (leaves follow).
        non_leaf_end = len(tree) + 4 + n_per_block * (key_size + 8)
        tree += bytes((0, 0)) + pack("<H", n_blocks)
        for i in range(n_per_block):
            if i < n_blocks:
                key = names[i * n_per_block].encode().ljust(key_size, b"\0")
                tree += key + pack("<Q", non_leaf_end + i * leaf_size)
            else:
                tree += b"\0" * key_size + pack("<Q", 0)
    # Leaf blocks (the actual key/value slots), indices in chrom order.
    j = 0
    for _block in range(n_blocks):
        n_here = min(n_keys - j, n_per_block)
        tree += bytes((1, 0)) + pack("<H", n_here)
        for _k in range(n_per_block):
            if j >= n_keys:
                tree += b"\0" * key_size + pack("<Q", 0)
            else:
                name, size_value = chrom_sizes[j]
                tree += name.encode().ljust(key_size, b"\0")
                tree += pack("<II", j, size_value)
                j += 1
    return bytes(tree)


# ---------------------------------------------------------------------------
# data sections: 24-byte section header + bedGraph items, zlib compressed
# ---------------------------------------------------------------------------

def _pack_data_sections(items, compresslevel):
    """Pack bedGraph items into zlib-compressed data sections."""
    sections = []
    section_items = []  # (start, end, value) triples of the open section
    section_tid = None

    def flush_section():
        nonlocal section_items
        if not section_items:
            return
        raw = bytearray(_SECTION_HEADER_SIZE)
        raw[0:4] = pack("<I", section_tid)
        raw[4:8] = pack("<I", section_items[0][0])  # start
        raw[8:12] = pack("<I", section_items[-1][1])  # end
        raw[12:16] = pack("<I", 0)  # itemStep
        raw[16:20] = pack("<I", 0)  # itemSpan
        raw[20] = 1  # bedGraph (type 1)
        raw[21] = 0  # padding byte
        raw[22:24] = pack("<H", len(section_items))
        for start, end, value in section_items:
            raw += pack("<IIf", start, end, value)
        sections.append(bytes(raw))
        section_items = []

    prev_tid = None
    for chrom_tid, start, end, value in items:
        if prev_tid is not None and chrom_tid != prev_tid:
            flush_section()
        if (
            section_items
            and _SECTION_HEADER_SIZE
            + (len(section_items) + 1) * _BEDGRAPH_ITEM_SIZE
            > UNCOMPRESSED_BLOCK_SIZE
        ):
            flush_section()
        section_tid = chrom_tid if not section_items else section_tid
        section_items.append((start, end, value))
        prev_tid = chrom_tid
    flush_section()

    data = bytearray()
    blocks = []
    for raw in sections:
        compressed = zlib.compress(raw, compresslevel)
        blocks.append(
            {
                "tid": int.from_bytes(raw[0:4], "little"),
                "start": int.from_bytes(raw[4:8], "little"),
                "end": int.from_bytes(raw[8:12], "little"),
                "raw_size": len(raw),
                "n_items": int.from_bytes(raw[22:24], "little"),
                "offset": len(data),
                "size": len(compressed),
            }
        )
        data += compressed
    return bytes(data), blocks


# ---------------------------------------------------------------------------
# total summary + R-tree index
# ---------------------------------------------------------------------------

def _pack_summary(items):
    n_bases = sum(end - start for _tid, start, end, _value in items)
    values = [value for _tid, _start, _end, value in items]
    if values:
        min_val, max_val = min(values), max(values)
    else:
        min_val = max_val = 0.0
    sum_data = sum((end - start) * value for _tid, start, end, value in items)
    sum_squared = sum(
        (end - start) * value * value for _tid, start, end, value in items
    )
    return pack("<Qdddd", n_bases, min_val, max_val, sum_data, sum_squared)


class _Node:
    __slots__ = ("is_leaf", "children", "entries")

    def __init__(self, is_leaf=False, children=None, entries=None):
        self.is_leaf = is_leaf
        self.children = children if children is not None else []
        self.entries = entries if entries is not None else []

    def bounds(self):
        if self.is_leaf:
            first, last = self.entries[0], self.entries[-1]
            return first["tid"], first["start"], last["tid"], last["end"]
        first = self.children[0].bounds()
        last = self.children[-1].bounds()
        return first[0], first[1], last[2], last[3]


def _node_size(node):
    if node.is_leaf:
        return 4 + 32 * len(node.entries)
    return 4 + 24 * len(node.children)


def _pack_rtree(blocks, index_start):
    """Build the R-tree index exactly like libBigWig's ``writeIndex``.

    ``index_start`` is the absolute file position where the index begins;
    internal-node child offsets must be absolute (libBigWig stores
    ``ftell`` results), while data-section offsets were already shifted
    into absolute coordinates by the caller.
    """
    if not blocks:
        return b""
    # Group data sections into leaf nodes of <= RTREE_BLOCK_SIZE entries
    # each, then keep grouping internal nodes until a single root
    # remains (same hierarchy libBigWig's addLeaves builds).
    level = [
        _Node(is_leaf=True, entries=blocks[i : i + RTREE_BLOCK_SIZE])
        for i in range(0, len(blocks), RTREE_BLOCK_SIZE)
    ]
    while len(level) > 1:
        level = [
            _Node(is_leaf=False, children=level[i : i + RTREE_BLOCK_SIZE])
            for i in range(0, len(level), RTREE_BLOCK_SIZE)
        ]
    root = level[0]

    # Pass 1: assign sequential positions to every node.  The root is
    # stored first (right after the 48-byte index header), then its
    # descendants level by level (BFS), matching writeIndexTreeNode.
    positions = {}
    current_level = [root]
    pos = _INDEX_HEADER_SIZE
    while current_level:
        next_level = []
        for node in current_level:
            positions[id(node)] = pos
            pos += _node_size(node)
            next_level.extend(node.children)
        current_level = next_level

    # Pass 2: serialize the index header and every node.
    index = bytearray()
    start_tid, start_base, end_tid, end_base = root.bounds()
    index += pack(
        "<IIQIIIIQII",
        RTREE_INDEX_MAGIC,
        RTREE_BLOCK_SIZE,  # blockSize (children per node)
        len(blocks),  # number of data sections
        start_tid,
        start_base,
        end_tid,
        end_base,
        pos,  # idxSize: total byte length of the index
        1,  # itemsPerSlot (libBigWig writes 1)
        0,  # padding u32
    )

    def make_node_bytes(node):
        out = bytearray()
        out += bytes((1 if node.is_leaf else 0, 0))
        count = len(node.entries) if node.is_leaf else len(node.children)
        out += pack("<H", count)
        if node.is_leaf:
            for entry in node.entries:
                out += pack(
                    "<IIIIQQ",
                    entry["tid"],
                    entry["start"],
                    entry["tid"],  # endTid == tid for single-chrom blocks
                    entry["end"],
                    entry["offset"],
                    entry["size"],
                )
        else:
            for child in node.children:
                child_tid, child_start, child_tid_end, child_end = child.bounds()
                out += pack(
                    "<IIII",
                    child_tid,
                    child_start,
                    child_tid_end,
                    child_end,
                )
                # dataOffset u64 completes the 24-byte internal entry and
                # must be an ABSOLUTE file position (libBigWig uses ftell)
                out += pack("<Q", index_start + positions[id(child)])
        return bytes(out)

    index += make_node_bytes(root)
    current_level = [root] if not root.is_leaf else []
    while current_level:
        next_level = []
        for node in current_level:
            for child in node.children:
                index += make_node_bytes(child)
                if not child.is_leaf:
                    next_level.append(child)
        current_level = next_level
    return bytes(index)


# ---------------------------------------------------------------------------
# in-place header patches
# ---------------------------------------------------------------------------

def _patch_u16(fh, offset, value):
    fh.seek(offset)
    fh.write(pack("<H", value))


def _patch_u32(fh, offset, value):
    fh.seek(offset)
    fh.write(pack("<I", value))


def _patch_u64(fh, offset, value):
    fh.seek(offset)
    fh.write(pack("<Q", value))


# ---------------------------------------------------------------------------
# CLI: bedGraph text + chrom sizes -> bigWig
# ---------------------------------------------------------------------------

def _main(argv=None):
    parser = argparse.ArgumentParser(
        description="Convert a bedGraph file (+ chrom sizes) into a bigWig "
        "file using only the Python standard library."
    )
    parser.add_argument("bedgraph", help="0-based 'chrom start end value' text file")
    parser.add_argument("chrom_sizes", help="one 'chrom size' pair per line")
    parser.add_argument("output", help="output .bw path")
    parser.add_argument(
        "--compress-level",
        type=int,
        default=9,
        help="zlib level for data blocks (default %(default)s)",
    )
    args = parser.parse_args(argv)

    chrom_sizes = []
    with open(args.chrom_sizes, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            name, size = line.split()[:2]
            chrom_sizes.append((name, int(size)))

    entries = []
    with open(args.bedgraph, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(("track", "browser", "#")):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            entries.append((parts[0], int(parts[1]), int(parts[2]), float(parts[3])))

    summary = write_bigwig(args.output, chrom_sizes, entries, args.compress_level)
    print(
        f"wrote {args.output}: {summary['n_entries']} entries in "
        f"{summary['n_blocks']} blocks ({summary['bytes']:,} bytes)"
    )
    return 0


if __name__ == "__main__":
    sys.exit(_main())
