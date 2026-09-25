"""Clone tree, per-haplotype copy-number state, and expected-VAF arithmetic for one dataset."""
from .intervals import IntervalIndex


class CloneModel:
    def __init__(self, ds_cfg, arms_bed):
        self.purity = float(ds_cfg["purity"])
        self.wgd = bool(ds_cfg.get("wgd", False))
        self.clones = ds_cfg["clones"]
        self.base = (2, 2) if self.wgd else (1, 1)
        self.arms = {}
        for line in open(arms_bed):
            c, s, e, arm = line.split()[:4]
            self.arms[f"{c}:{arm}"] = (c, int(s), int(e))
        children = {c: [k for k, v in self.clones.items() if v["parent"] == c] for c in self.clones}
        # exclusive fraction: tumor cells whose most-derived clone is c
        self.excl = {c: self.clones[c]["ccf"] - sum(self.clones[k]["ccf"] for k in children[c]) for c in self.clones}
        if any(v < -1e-9 for v in self.excl.values()):
            raise ValueError("child CCFs exceed their parent's CCF")
        self.seg = {}
        for c, evs in (ds_cfg.get("cna") or {}).items():
            ix = IntervalIndex()
            for ev in evs:
                chrom, s, e = self.region(ev["region"])
                ix.add(chrom, s, e, (tuple(ev["cn"]), ev.get("label", "")))
            ix.build()
            self.seg[c] = ix
        self.lineage = {c: self._lineage(c) for c in self.clones}

    def region(self, r):
        if r in self.arms:
            return self.arms[r]
        chrom, span = r.split(":")
        s, e = span.split("-")
        return chrom, int(s), int(e)

    def _lineage(self, c):
        out = []
        while c is not None:
            out.append(c)
            c = self.clones[c]["parent"]
        return out[::-1]

    def is_descendant(self, c, anc):
        return anc in self.lineage[c]

    def cn(self, clone, chrom, pos1):
        """(cnA, cnB, label) at a 1-based position in a clone; the most-derived event in the lineage wins."""
        state, label = self.base, "baseline"
        for c in self.lineage[clone]:
            ix = self.seg.get(c)
            if ix is None:
                continue
            for _s, _e, (cnab, lab) in ix.overlaps(chrom, pos1 - 1, pos1):
                state, label = cnab, lab
        return state[0], state[1], label

    def vaf(self, chrom, pos1, hap, clone, pre_cna=True):
        """Expected bulk-tumor VAF for an event on haplotype `hap` acquired in `clone`.
        Truncal events with pre_cna=True predate the copy-number changes (multiplicity = current copies of
        that haplotype, which also models pre-WGD timing); pre_cna=False or subclonal events have multiplicity 1
        where the haplotype is retained, 0 where it is lost.
        Returns (vaf, {clone: multiplicity}, mean tumor copy number at the locus)."""
        num = 0.0
        den_t = 0.0
        mults = {}
        for c in self.clones:
            f = self.excl[c]
            if f <= 0:
                continue
            a, b, _ = self.cn(c, chrom, pos1)
            cnh = (a, b)[hap]
            den_t += f * (a + b)
            if self.is_descendant(c, clone):
                m = cnh if (clone == "T" and pre_cna) else (1 if cnh >= 1 else 0)
                mults[c] = m
                num += f * m
        p = self.purity
        den = p * den_t + (1 - p) * 2
        return (p * num / den if den > 0 else 0.0), mults, den_t

    def ccf(self, clone):
        return self.clones[clone]["ccf"]

    def clonality_tier(self, clone, chrom, pos1, hap):
        if clone != "T":
            return clone
        a, b, _ = self.cn("T", chrom, pos1)
        cnh = (a, b)[hap]
        if cnh >= 4:
            return "T_amp"
        if (a == 0) != (b == 0):
            return "T_LOH"
        return "T_het"

    def retained_haplotypes(self, clone, chrom, pos1):
        """Haplotype indices with at least one copy in `clone` at the locus."""
        a, b, _ = self.cn(clone, chrom, pos1)
        return [h for h, c in enumerate((a, b)) if c >= 1]

    def regions_with(self, predicate, clone="T"):
        """List (chrom, start0, end0, label) of truncal-lineage CN segments satisfying predicate(cnA, cnB)."""
        out = []
        for c in self.lineage[clone]:
            ix = self.seg.get(c)
            if ix is None:
                continue
            for chrom in ix._starts:
                for s, e, (cnab, lab) in zip(ix._starts[chrom], ix._ends[chrom], ix._names[chrom]):
                    if predicate(*cnab):
                        out.append((chrom, s, e, lab))
        return out
