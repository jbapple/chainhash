#!/usr/bin/env python3
"""Recompute test/128v2/vectors.txt with the independent Python oracle (pyref_2l.py) and compare.
Default: the lengths up to 4200 bytes (seconds). --full: every length, and the regenerated file
must equal vectors.txt byte for byte (a few minutes). The archive is never overwritten."""
import io, os, sys
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import pyref_2l

full = '--full' in sys.argv[1:]
archive = open(os.path.join(here, 'vectors.txt')).read()
if full:
    out = io.StringIO()
    key = bytes((i * 73 + 11) & 255 for i in range(192))
    lens = pyref_2l.lengths(16)
    msg = bytes((i * 137 + 29) & 255 for i in range(max(lens)))
    print('# CH-128/P v2 (M=16: 2048 B blocks, 8 blocks per region, level-2 regions of R=8). key[i]=(i*73+11)&255 (192 bytes), msg[i]=(i*137+29)&255', file=out)
    print('# len digest(16 bytes, little-endian hex)', file=out)
    for L in lens:
        print('%d %s' % (L, pyref_2l.hash_2l(16, key, msg[:L]).to_bytes(16, 'little').hex()), file=out)
    assert out.getvalue() == archive, 'the Python oracle does not reproduce test/128v2/vectors.txt'
    print('PASS Python oracle reproduces test/128v2/vectors.txt (%d lengths, byte for byte)' % len(lens))
else:
    key = bytes((i * 73 + 11) & 255 for i in range(192))
    rows = [l.split() for l in archive.splitlines() if l and not l.startswith('#')]
    rows = [(int(a), b) for a, b in rows if int(a) <= 4200]
    msg = bytes((i * 137 + 29) & 255 for i in range(4200))
    for L, want in rows:
        got = pyref_2l.hash_2l(16, key, msg[:L]).to_bytes(16, 'little').hex()
        assert got == want, 'length %d: oracle %s, archive %s' % (L, got, want)
    print('PASS Python oracle matches test/128v2/vectors.txt at %d lengths <= 4200 B (--full: all)' % len(rows))
