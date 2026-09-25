"""Expression baseline: per-transcript TPM table (GENCODE v23 IDs from Xena) mapped onto the working annotation."""
import gzip, math

class Expression:
    def __init__(self, tsv_gz, tx):
        self.tx_tpm = {}   # versionless ENST -> median TPM
        with gzip.open(tsv_gz, "rt") as fh:
            hdr = fh.readline().rstrip("\n").split("\t"); i_med = hdr.index("median_tpm")
            for line in fh:
                f = line.rstrip("\n").split("\t"); self.tx_tpm[f[0].split(".")[0]] = float(f[i_med])
        # gene TPM = sum over that gene's transcripts present in the table; transcripts absent inherit the gene mean
        self.gene_tpm = {}; per_gene = {}
        for t in tx.values():
            v = self.tx_tpm.get(t.tid.split(".")[0])
            per_gene.setdefault(t.gene_id, []).append(v)
        for g, vs in per_gene.items():
            known = [v for v in vs if v is not None]
            self.gene_tpm[g] = sum(known) if known else 0.0
        self._per_gene = per_gene
    def transcript(self, tid, gene_id=None):
        v = self.tx_tpm.get(tid.split(".")[0])
        if v is not None: return v
        if gene_id and gene_id in self._per_gene:
            known = [x for x in self._per_gene[gene_id] if x is not None]
            return (sum(known) / len(known)) if known else 0.0
        return 0.0
    def gene(self, gene_id): return self.gene_tpm.get(gene_id, 0.0)
    @staticmethod
    def tier(tpm):
        """Expression tiers used by the design: 0 (<0.5), 1 (~1), 10, 100, 1000 TPM."""
        if tpm < 0.5: return "T0"
        if tpm < 3: return "T1"
        if tpm < 30: return "T10"
        if tpm < 300: return "T100"
        return "T1000"
