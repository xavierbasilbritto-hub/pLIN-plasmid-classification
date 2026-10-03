# Copyright (C) 2025 Basil Xavier Britto. License: GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
FracMinHash sketches of canonical k-mers (numpy only, deterministic).

A sketch keeps every canonical k-mer whose 64-bit hash falls below
2**64 / scaled, so the containment of A in B is estimated by
|sketch(A) ∩ sketch(B)| / |sketch(A)| (Irber et al. 2022, sourmash's
FracMinHash). k-mers containing non-ACGT characters are skipped. The hash is
splitmix64 of the 2-bit packed canonical k-mer, so sketches are identical on
every platform and need no external tool.
"""

import numpy as np

K = 21
SCALED = 200
_CODE = np.full(256, 255, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    _CODE[_b] = _i
    _CODE[ord(chr(_b).lower())] = _i
_MASK = np.uint64((1 << (2 * K)) - 1)


def _cutoff(scaled):
    """Largest kept hash at this scale (scale 1 keeps every k-mer)."""
    return np.uint64((2 ** 64 - 1) // scaled)


def _splitmix64(x):
    x = x + np.uint64(0x9E3779B97F4A7C15)
    x = (x ^ (x >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
    x = (x ^ (x >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
    return x ^ (x >> np.uint64(31))


def sketch(sequence, k=K, scaled=SCALED):
    """Sorted unique uint64 hashes of the kept canonical k-mers."""
    c = _CODE[np.frombuffer(sequence.encode("ascii", "replace"), dtype=np.uint8)]
    n = len(c) - k + 1
    if n <= 0:
        return np.zeros(0, dtype=np.uint64)
    valid = c != 255
    # a window is usable only if all k bases are ACGT
    bad = np.concatenate([[0], np.cumsum(~valid)])
    ok = (bad[k:] - bad[:-k]) == 0
    c64 = np.where(valid, c, 0).astype(np.uint64)
    fwd = np.zeros(n, dtype=np.uint64)
    rev = np.zeros(n, dtype=np.uint64)
    with np.errstate(over="ignore"):
        for j in range(k):
            fwd = (fwd << np.uint64(2)) | c64[j:j + n]
            rev = rev | ((np.uint64(3) - c64[j:j + n]) << np.uint64(2 * j))
        canon = np.minimum(fwd & _MASK, rev & _MASK)[ok]
        h = _splitmix64(canon)
    keep = h <= _cutoff(scaled)
    return np.unique(h[keep])


TARGET = 400        # aim for >= 400 kept k-mers per plasmid
MAX_SCALED = 256


def adaptive_scaled(length):
    """Power-of-two scale giving ~TARGET kept k-mers: small plasmids keep proportionally more."""
    s = 1
    while s * 2 <= MAX_SCALED and length / (s * 2) >= TARGET:
        s *= 2
    return s


def adaptive_sketch(sequence):
    """(hashes, scaled) with the adaptive scale for this sequence length."""
    s = adaptive_scaled(len(sequence))
    return sketch(sequence, scaled=s), s


def downsample(h, scaled):
    """Exact FracMinHash downsampling to a coarser scale."""
    return h[h <= _cutoff(scaled)]


def _common(a, sa, b, sb):
    s = max(sa, sb)
    a, b = (downsample(a, s) if sa < s else a), (downsample(b, s) if sb < s else b)
    return a, b, np.intersect1d(a, b, assume_unique=True).size


def containment(a, b, sa=SCALED, sb=SCALED):
    """Estimated fraction of A's k-mers present in B (compared at the coarser of the two scales)."""
    a, b, n = _common(a, sa, b, sb)
    return n / a.size if a.size else 0.0


def min_containment(a, b, sa=SCALED, sb=SCALED):
    """Shared k-mers over the larger sketch: high only if both plasmids are mostly shared (symmetric)."""
    a, b, n = _common(a, sa, b, sb)
    m = max(a.size, b.size)
    return n / m if m else 0.0
