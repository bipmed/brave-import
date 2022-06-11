Import variant data from VCF file into BraVE server

This is the newest implementation of brave-import that solves issues regarding missing FORMAT columns. It is faster than the earlier implementations and uses an official VCF parser (htslib).

Download and install Rust: https://www.rust-lang.org/tools/install
Clone this repository and then build the binary file.

```bash
cargo build --release
```

Binary file will be in `target/release/brave-import`, you might want to copy to somewhere in your PATH.

There are INFO columns required to build BraVE database properly, make sure that your VCF file has the following columns:

- `NS` - Number of samples with data
- `AF` - Allele frequency
- `ANN` - Standard annotation format. Added by snpEff or ClinEff
- `CLNSIG` - Variant clinical significance. Added by by snpEff or ClinEff

It accepts VCF files (v4.2) as input and submit variants to server instance. No genotype (FORMAT column) data is sent to server. FORMAT/DP and FORMAT/GQ are used to calculate distribution (min, q25, median, q75, max and average) of every variant. By default only variant that passed all filters are imported to database (FILTER = PASS or .). Use `--dont-filter` option to import all variants, regardless of FILTER column.

```bash
brave import \
    [--dont-filter] \
    [--dryrun] \
    [--verbose] \
    [--host http://localhost:8080] \
    [--username admin] \
    [--password secret] \
    --assembly hg38 \
    --dataset bipmed \
    bipmed.hg38.vcf.gz
```