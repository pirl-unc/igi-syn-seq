"""Human-readable 'diff cards': the patient's two germline haplotypes at a locus, the mutant haplotype with the
somatic change marked, and the wild-type vs mutant protein window. One card per designed event."""

FLANK_NT = 20
FLANK_AA = 12


def _mark(width, offset, length=1):
    return " " * offset + "^" * max(1, length)


def dna_card(genome, germline, chrom, pos1, ref, alt, hap):
    """Return (lines, wt_window, mut_window) for a substitution/indel applied on haplotype `hap`."""
    s1, e1 = pos1 - FLANK_NT, pos1 + len(ref) - 1 + FLANK_NT
    h0, h1, m0, m1 = germline.haplotype_seqs(genome, chrom, s1, e1, extra=(pos1, ref, alt, hap))
    refseq = genome.seq(chrom, s1 - 1, e1)
    off = pos1 - s1
    lines = [
        f"  window    {chrom}:{s1}-{e1}  somatic {ref}>{alt} at {pos1} on hap{hap}",
        f"  ref       {refseq}",
        f"  hap0      {h0}" + ("   <- germline differs" if h0 != refseq else ""),
        f"  hap1      {h1}" + ("   <- germline differs" if h1 != refseq else ""),
        f"  hap{hap}+mut  {(m0, m1)[hap]}",
        f"            {_mark(len(refseq), off, len(alt))}",
    ]
    return lines, (h0, h1)[hap], (m0, m1)[hap]


def protein_card(wt_prot, mut_prot, k, cons):
    """Aligned wild-type vs mutant protein windows around the first changed residue k."""
    if k is None or wt_prot is None:
        return ["  protein   n/a"], "", ""
    a, b = max(0, k - FLANK_AA), k + FLANK_AA + 1
    wt_w = wt_prot[a:b]
    if cons == "frameshift":
        mut_w = mut_prot[a:k] + mut_prot[k:k + 2 * FLANK_AA + 1]  # show the whole neo-ORF start
    else:
        mut_w = mut_prot[a:b]
    lines = [
        f"  protein   aa {a + 1}-{b}  first change at aa {k + 1} ({cons})",
        f"  WT        {wt_w}",
        f"  MUT       {mut_w}",
        f"            {_mark(len(wt_w), k - a)}",
    ]
    return lines, wt_w, mut_w


def render(event_id, title, dna_lines, prot_lines, extra=None):
    out = [f"## {event_id}  {title}"] + dna_lines + prot_lines
    if extra:
        out += [f"  {k:9s} {v}" for k, v in extra.items()]
    return "\n".join(out) + "\n"
