# External data sources for IGI-SYN-SEQ

Path placeholders: `$RAFT` = RAFT workspace root, `$DATASETS` = shared lab datasets directory.

| Purpose | Source | Location on disk | Status (2026-09-24) |
|---|---|---|---|
| HG002 phased germline baseline (IGI-SYN-SEQ-01) | GIAB HG2-T2TQ100-V1.1 GRCh38 smvar + stvar (defrabb V0.020) | `$DATASETS/HG002/GIAB_Q100_benchmark_defrabbV0.020/` | downloaded 2026-09-24 |
| HiFi WGS error model | GIAB HG002 Revio HiFi raw movies (m84039_230928_213653_s3, m84039_231005_222902_s1) | `$DATASETS/HG002/PacBio_HiFi-Revio_20231031/` | downloaded 2026-09-24 |
| Kinnex sc error / MAS-seq model | PacBio DATA-Revio-Kinnex-HG002-10x5p (pre-skera BAM + matched 10x Illumina) | `$DATASETS/HG002/Kinnex_scRNA_HG002_10x5p/` | blocked: downloads.pacbcloud.com unreachable from cluster; see `datasets/HG002/NEEDS_EXTERNAL_DOWNLOAD.md` |
| Kinnex bulk error / MAS-seq array model | PacBio Kinnex-full-length-RNA DATA-Revio-HG002-1 (post-skera `segmented.bam`, 38.7 M reads, with skera `ds`/`di`/`dl`/`zm` tags that let the 8-segment arrays be reconstructed per ZMW; `flnc.bam`, 37.2 M reads) | `$DATASETS/HG002/Kinnex_fullLengthRNA_HG002_Revio/` | copied 2026-09-24 from a local institutional mirror of downloads.pacbcloud.com |
| Extra bulk long-read RNA (HG002) | GIAB RNA effort MAS-seq, GM24385 refined FLNC FASTQ + lima/refine reports | `$DATASETS/HG002/GIAB_MASseq_RNA_GM24385/` | copied 2026-09-24, same mirror |
| IPISRC044 germline baseline (IGI-SYN-SEQ-02) | DeepVariant WGS on blood-normal BAM + WhatsHap (short reads + T1 ONT) | `$RAFT/inputs/vcfs/IPISRC044/germline_wgs/` | in progress (2026-09-24) |
| IPISRC044 short-read error models | real IPISRC044 WES / RNA / 10x FASTQs | `$RAFT/inputs/fastqs/IPISRC044/` | present |
| Immune/stromal single-cell profiles | IPISRC044 UCSF SCG1 (10x 5') | `$RAFT/inputs/fastqs/IPISRC044/assays/ucsf/T*/FASTQ/` | present |
| Tumor expression baseline | UCSC Xena Toil `tcga_Kallisto_tpm` (transcript log2 TPM, all TCGA) + PanCanAtlas subtypes `tcga_subtypes_syn8402849.csv` (BRCA PAM50 in `Subtype_mRNA`) | `$DATASETS/ucsc_xena/`, `$DATASETS/TCGA/tcga_subtypes_syn8402849/` | present |
| ERV / CTA pan-normal expression reference | lens-v2.0.0-dev-alt `pan_normal_and_mtec_exp.95th_perc.homo_sapiens.quant.sf`, `erv_pep_exp_norm_and_mtec.homo_sapiens.tsv` | `$RAFT/references/homo_sapiens/{expression,erv}/` | present (to be regenerated later per issue 31) |
| Statistical phasing panel | 1000 Genomes 30x GRCh38 phased panel (for SHAPEIT5) | `$DATASETS/1000G_GRCh38_phased/` | to download |

Note: LENS's own germline VCFs for IPISRC044 (lens-v2.0.0-dev, IPISRC044-WGS project) are
exome-restricted (`--regions hg38_exome.50bpflank.bed`, 156k records) and phased from
short reads only (block N50 138 bp), so they are not usable as a genome-wide phased baseline.

## Containers (all pulled into `$IMAGE_DIR` with `singularity pull docker://...`)

| Tool | Image |
|---|---|
| DeepVariant 1.8.0 | `google/deepvariant:1.8.0` |
| WhatsHap 2.4 | `quay.io/biocontainers/whatshap:2.4--py310h184ae93_0` |
| SHAPEIT5 5.1.1 | `quay.io/biocontainers/shapeit5:5.1.1--h34261f4_2` |
| sniffles 2.8.0 | `quay.io/biocontainers/sniffles:2.8.0--pyhdfd78af_1` |
| minimap2 2.28 / samtools 1.21 / bcftools 1.21 | `quay.io/biocontainers/{minimap2:2.28--h577a1d6_4,samtools:1.21--h50ea8bc_0,bcftools:1.21--h8b25389_0}` |
| pbsim3 3.0.5 | `quay.io/biocontainers/pbsim3:3.0.5--h9948957_2` |
| ART 2016.06.05 | `quay.io/biocontainers/art:2016.06.05--h0704011_13` |
| skera 1.4.0 / pbtk 3.5.0 / lima 2.13 / isoseq 4.3 / pbmm2 1.17 / pigeon 1.4 | `quay.io/biocontainers/{pbskera:1.4.0--hdfd78af_0,pbtk:3.5.0--h9ee0642_0,lima:2.13.0--h9ee0642_0,isoseq:4.3.0--h9ee0642_0,pbmm2:1.17.0--h9ee0642_0,pbpigeon:1.4.0--h9948957_0}` |
| Cell Ranger 9.0.1 | `quay.io/nf-core/cellranger:9.0.1` (community build of the 10x tarball; 10x EULA applies) |
| netMHCpan 4.1b / mhcflurry 2.1.1 | LENS images `spvensko/netmhcpan:4.1b`, `spvensko/mhcflurry:2.1.1` |
