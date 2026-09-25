"""Minimal interval index over BED-like records (per-chromosome sorted starts + bisect)."""
import bisect, gzip
from collections import defaultdict

class IntervalIndex:
    def __init__(self):
        self._starts = defaultdict(list); self._ends = defaultdict(list); self._names = defaultdict(list)
        self._maxlen = defaultdict(int); self._built = False
    def add(self, chrom, start, end, name=None):
        self._starts[chrom].append(start); self._ends[chrom].append(end); self._names[chrom].append(name); self._built = False
    def build(self):
        for c in self._starts:
            order = sorted(range(len(self._starts[c])), key=lambda i: self._starts[c][i])
            self._starts[c] = [self._starts[c][i] for i in order]; self._ends[c] = [self._ends[c][i] for i in order]
            self._names[c] = [self._names[c][i] for i in order]
            self._maxlen[c] = max((e - s for s, e in zip(self._starts[c], self._ends[c])), default=0)
        self._built = True
    def overlaps(self, chrom, start, end):
        """Yield (start, end, name) of records overlapping [start, end) (0-based half-open)."""
        if not self._built: self.build()
        S = self._starts.get(chrom); 
        if not S: return
        lo = bisect.bisect_left(S, start - self._maxlen[chrom]); hi = bisect.bisect_right(S, end)
        for i in range(lo, hi):
            if self._ends[chrom][i] > start and S[i] < end:
                yield S[i], self._ends[chrom][i], self._names[chrom][i]
    def any(self, chrom, start, end):
        return next(self.overlaps(chrom, start, end), None) is not None
    def distance_to_edge(self, chrom, pos):
        """Distance from a 0-based position to the nearest boundary of the record containing it (None if outside)."""
        best = None
        for s, e, _ in self.overlaps(chrom, pos, pos + 1):
            d = min(pos - s, e - 1 - pos); best = d if best is None else min(best, d)
        return best
    @classmethod
    def from_bed(cls, path, name_col=3):
        ix = cls(); op = gzip.open if str(path).endswith(".gz") else open
        with op(path, "rt") as fh:
            for line in fh:
                if not line.strip() or line.startswith(("#", "track", "browser")): continue
                f = line.rstrip("\n").split("\t")
                ix.add(f[0], int(f[1]), int(f[2]), f[name_col] if len(f) > name_col else None)
        ix.build(); return ix
