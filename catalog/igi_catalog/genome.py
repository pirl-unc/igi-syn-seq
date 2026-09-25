"""Reference genome access and sequence-context helpers."""
import pysam

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def revcomp(s): return s.translate(COMP)[::-1]

class Genome:
    def __init__(self, fasta):
        self.fa = pysam.FastaFile(fasta)
        self.lengths = dict(zip(self.fa.references, self.fa.lengths))
    def seq(self, chrom, start, end):
        """0-based half-open, uppercase."""
        return self.fa.fetch(chrom, max(0, start), min(end, self.lengths[chrom])).upper()
    def homopolymer_run(self, chrom, pos0):
        """Length of the longest homopolymer run touching 0-based position pos0 (run of the base at pos0 or its neighbours)."""
        s = self.seq(chrom, pos0 - 12, pos0 + 13); c = 12
        best = 1
        for anchor in (c - 1, c, c + 1):
            if anchor < 0 or anchor >= len(s): continue
            b = s[anchor]; i = anchor; j = anchor
            while i > 0 and s[i - 1] == b: i -= 1
            while j < len(s) - 1 and s[j + 1] == b: j += 1
            best = max(best, j - i + 1)
        return best
    def gc(self, chrom, start, end):
        s = self.seq(chrom, start, end); n = len(s)
        return (s.count("G") + s.count("C")) / n if n else 0.0
