"""GENCODE GTF parsing to a compact transcript model, with CDS/protein extraction and coordinate mapping."""
import gzip, json, os, re
from dataclasses import dataclass, field, asdict
from .genome import revcomp

CODON = {a + b + c: aa for (a, b, c), aa in zip(
    [(x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"],
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}
def translate(cds):
    return "".join(CODON.get(cds[i:i + 3], "X") for i in range(0, len(cds) - len(cds) % 3, 3))

@dataclass
class Transcript:
    tid: str; gene_id: str; gene_name: str; gene_type: str; chrom: str; strand: str
    start: int; end: int  # 1-based inclusive transcript span
    exons: list = field(default_factory=list)   # [(start,end)] 1-based inclusive, genomic order
    cds: list = field(default_factory=list)     # [(start,end)] 1-based inclusive, genomic order (no stop codon in GENCODE CDS)
    tags: list = field(default_factory=list)
    level: int = 0
    def cds_len(self): return sum(e - s + 1 for s, e in self.cds)
    def n_exons(self): return len(self.exons)

_ATTR = re.compile(r'(\S+) "([^"]*)"')
def parse_gtf(path, keep_types=("protein_coding", "lncRNA")):
    """Return {tid: Transcript} for genes of keep_types (None = all)."""
    op = gzip.open if path.endswith(".gz") else open
    tx = {}
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"): continue
            f = line.rstrip("\n").split("\t")
            if f[2] not in ("transcript", "exon", "CDS"): continue
            a = dict(_ATTR.findall(f[8]))
            gt = a.get("gene_type", "")
            if keep_types and gt not in keep_types: continue
            tid = a["transcript_id"]
            if f[2] == "transcript":
                tags = re.findall(r'tag "([^"]+)"', f[8])
                tx[tid] = Transcript(tid, a["gene_id"], a.get("gene_name", ""), gt, f[0], f[6], int(f[3]), int(f[4]), tags=tags, level=int(a.get("level", 0)))
            elif tid in tx:
                (tx[tid].exons if f[2] == "exon" else tx[tid].cds).append((int(f[3]), int(f[4])))
    for t in tx.values():
        t.exons.sort(); t.cds.sort()
    return tx

def load_or_build(gtf, cache):
    if cache and os.path.exists(cache):
        with gzip.open(cache, "rt") as fh:
            return {k: Transcript(**v) for k, v in json.load(fh).items()}
    tx = parse_gtf(gtf)
    if cache:
        with gzip.open(cache, "wt") as fh:
            json.dump({k: asdict(v) for k, v in tx.items()}, fh)
    return tx

def representative_transcripts(tx):
    """One coding transcript per protein-coding gene: MANE_Select, else Ensembl_canonical, else longest CDS among 'basic'."""
    by_gene = {}
    for t in tx.values():
        if t.gene_type != "protein_coding" or not t.cds or "cds_start_NF" in t.tags or "cds_end_NF" in t.tags: continue
        by_gene.setdefault(t.gene_id, []).append(t)
    rep = {}
    for g, ts in by_gene.items():
        def key(t): return (("MANE_Select" in t.tags), ("Ensembl_canonical" in t.tags), ("basic" in t.tags), t.cds_len())
        rep[g] = max(ts, key=key)
    return rep

class CodingModel:
    """CDS sequence, protein, and genomic<->CDS coordinate mapping for one transcript."""
    def __init__(self, t: Transcript, genome):
        self.t = t; self.g = genome
        segs = t.cds if t.strand == "+" else t.cds[::-1]
        self.map = []  # genomic 1-based positions in CDS order
        for s, e in segs:
            rng = range(s, e + 1) if t.strand == "+" else range(e, s - 1, -1)
            self.map.extend(rng)
        seq = "".join(genome.seq(t.chrom, s - 1, e) for s, e in t.cds)
        self.cds = seq if t.strand == "+" else revcomp(seq)
        # append the stop codon (GENCODE CDS excludes it) when in frame
        self.protein = translate(self.cds)
        self.pos_index = {p: i for i, p in enumerate(self.map)}
    def cds_index(self, gpos): return self.pos_index.get(gpos)
    def codon_number(self, gpos):
        i = self.cds_index(gpos); return None if i is None else i // 3
    def exon_number(self, gpos):
        for n, (s, e) in enumerate(self.t.exons if self.t.strand == "+" else self.t.exons[::-1], start=1):
            if s <= gpos <= e: return n
        return None
    def is_last_coding_exon(self, gpos):
        segs = self.t.cds if self.t.strand == "+" else self.t.cds[::-1]
        s, e = segs[-1]; return s <= gpos <= e
    def mutate_protein(self, gpos, ref, alt):
        """Apply a simple substitution/indel at genomic position (1-based, ref/alt on + strand) and return
        (mutant_protein, first_changed_aa_index, consequence). Handles SNV, insertion, deletion within CDS."""
        i = self.cds_index(gpos)
        if i is None: return None, None, "non_coding"
        cds = self.cds
        if self.t.strand == "-":
            ref_t, alt_t = revcomp(ref), revcomp(alt)
            # on minus strand the anchor base is the last base of the reversed allele
            i = i - (len(ref_t) - 1)
        else:
            ref_t, alt_t = ref, alt
        if cds[i:i + len(ref_t)] != ref_t: return None, None, "ref_mismatch"
        mut = cds[:i] + alt_t + cds[i + len(ref_t):]
        # extend with downstream genomic sequence so frameshifts can run into the 3'UTR
        tail = self._downstream(3000)
        mprot = translate(mut + tail)
        wprot = self.protein
        stop = mprot.find("*"); mprot = mprot[:stop] if stop >= 0 else mprot
        k = 0
        while k < min(len(mprot), len(wprot)) and mprot[k] == wprot[k]: k += 1
        d = len(alt_t) - len(ref_t)
        if d == 0 and len(ref_t) == 1:
            if k == len(wprot) and len(mprot) == len(wprot): cons = "synonymous"
            elif k < len(wprot) and len(mprot) == k: cons = "nonsense" if k < len(wprot) else "synonymous"
            elif k == 0 and mprot[:1] != "M": cons = "start_loss"
            elif len(mprot) > len(wprot): cons = "stop_loss"
            else: cons = "missense" if len(mprot) == len(wprot) else "nonsense"
        elif d % 3 == 0: cons = "inframe_insertion" if d > 0 else "inframe_deletion"
        else: cons = "frameshift"
        return mprot, k, cons
    def _downstream(self, n):
        t = self.t
        if t.strand == "+":
            last = t.cds[-1][1]; return self.g.seq(t.chrom, last, last + n)
        first = t.cds[0][0]; return revcomp(self.g.seq(t.chrom, first - 1 - n, first - 1))
