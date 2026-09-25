"""Sequence-context strata for a candidate site."""
from .intervals import IntervalIndex


class ContextAnnotator:
    def __init__(self, genome, germline, exome_bed, rmsk_bed, segdups_bed):
        self.g = genome
        self.gl = germline
        import pysam
        self.exome = IntervalIndex.from_bed(exome_bed)
        self.rmsk = pysam.TabixFile(rmsk_bed)   # 5.6 M records: query through the tabix index instead of loading
        self.segdup = IntervalIndex.from_bed(segdups_bed)

    def _rmsk_names(self, chrom, p0):
        try:
            return [row.split("\t")[3] for row in self.rmsk.fetch(chrom, p0, p0 + 1)]
        except ValueError:
            return []

    def features(self, chrom, pos1, cm=None):
        """Raw context features at a 1-based position (cm: CodingModel of the transcript, optional)."""
        p0 = pos1 - 1
        hp = self.g.homopolymer_run(chrom, p0)
        hets = self.gl.hets_near(chrom, pos1, 30)
        edge = self.exome.distance_to_edge(chrom, p0)
        rep = self._rmsk_names(chrom, p0)
        segdup = self.segdup.any(chrom, p0, p0 + 1)
        cds_edge = None
        if cm is not None:
            for s, e in cm.t.cds:
                if s <= pos1 <= e:
                    cds_edge = min(pos1 - s, e - pos1)
        return {
            "homopolymer_run": hp,
            "germline_hets_within_30bp": len(hets),
            "nearest_germline_het_bp": min((abs(h[0] - pos1) for h in hets), default=None),
            "capture_edge_bp": edge,           # None when outside the capture BED
            "on_target": edge is not None,
            "repeat": rep[0] if rep else None,
            "segdup": segdup,
            "cds_edge_bp": cds_edge,
            "gc_100bp": round(self.g.gc(chrom, p0 - 50, p0 + 50), 3),
        }

    @staticmethod
    def stratum(f):
        """Assign the single context stratum a site belongs to (priority order)."""
        if not f["on_target"]:
            return "deep_intronic"
        if f["segdup"]:
            return "segdup"
        if f["homopolymer_run"] >= 6:
            return "homopolymer"
        if f["germline_hets_within_30bp"] > 0:
            return "near_germline_het"
        if f["cds_edge_bp"] is not None and f["cds_edge_bp"] <= 10:
            return "exon_edge"
        if f["repeat"] is not None:
            return "repeat_other"
        return "clean"
