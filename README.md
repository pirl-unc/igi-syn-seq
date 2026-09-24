# igi-syn-seq

Generator and documentation for **IGI-SYN-SEQ** (datasets IGI-SYN-SEQ-01, HG002 baseline, and IGI-SYN-SEQ-02, IPISRC044 baseline), two fully synthetic, known-truth multi-assay
tumor/normal datasets (Illumina WES/RNA/10x 5' GEX + TCR, PacBio HiFi WGS, Kinnex bulk and
single-cell) built to benchmark LENS across every antigen class.

- `docs/IGI-SYN-SEQ-design.md`: design specification (baselines, TNBC tumor architecture,
  variant catalog and evidence tiers, assays, truth bundle, chr1to6 subset).
- `docs/data-sources.md`: external inputs and where they live.
- `docs/baseline-references.md`: how the phased germline baselines are built (DeepVariant, WhatsHap,
  SHAPEIT5, sniffles2, GIAB Q100), written with placeholders so it can be followed on any cluster.

Generator code will be added once the design settles; for now this repository holds the design and
procedures only.
