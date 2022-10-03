# Sample Contamination Estimate from VCF (SCE-VCF)

SCE-VCF is a simple tool to evalute sample contamination from sequencing experiments. It uses WGS/WES VCF files to compute few metrics useful to detect potential sample contamination:

- CHARR: ref read fraction in hom alt vars as proposed by Wenham Lu, Broad Institute
- Mean ref allele AB observed in hom alt genotypes
- heterozygosity ratio (Nhet / Nhom)
- rate on inconsistent het calls (het call with AB outside threshold)

The tool is better suited for cohort analysis so one can compare value distributions to identify outliers. Even a small cohort of 20 samples work fine in our test.

## Usage

```bash
Usage:
  sceVCF [options] [vcf ...]

Arguments:
  [vcf ...]        input VCF/BCF file(s). Glob pattern allowed

Options:
  -h, --help
  -o, --out=OUT              Output file (TSV). If not provided output to stdout
  -f, --ad_field=AD_FIELD    FORMAT field containing the AD values (default: AD)
  -a, --af_field=AF_FIELD    INFO field containing the allele frequency (default: AF)
  -l, --het_ab_limit=HET_AB_LIMIT
                             Comma separated min and max allele balance accepted for het calls (default: 0.25,0.75)
  -q, --min_GQ=MIN_GQ        Min GQ for hom var to be included in charr computation (default: 20)
  -d, --min_DP=MIN_DP        Min DP for var to be included in computation (default: 20)
  -s, --samples=SAMPLES      Restrict analysis to the given samples. Comma separated list or file wiht one sample per line
  -r, --region=REGION        Specify genomic region for processing. Format is chr[:start-end]. Comma-separated list of regions or file with 1 region per line.
```

## How it works

1. any variants with DP less than `min_DP` will be skipped.
2. variants with GQ less than `min_GQ` will be used only to compute inconsistent het rate, but discarded when computed het ratio and CHARR
3. CHARR is computed as `ADref / (AFref * ADtot)` considering only high-quality hom sites with ALT AF < 0.95, to avoid issues with very low AFref

## Notes

- **AF data:** The tool needs an annotated AF field in the input VCF to compute CHARR values. Ideally, this should be AF computed from the samples in the VCF, so the tool is better suited for cohort analysis. Even a small cohort of 20 samples work fine in our test.

- **Multiple input VCFs:** The tool can read multiple VCFs and statistics are then aggregated by sample IDs if there are matching sample IDs across the input VCFs. In this way it's easy to analyze cohort VCFs splitted by chromosome for example.

- **Subset samples:** User can specify a subset of samples to analyze in the cohort VCF, this can be conveniently used also to parallelize the analysis for large cohorts. Note that some time is required anyway to parse sample initially, so run time will increase as more samples are present in the input VCF

- **Subset regions** Analysis can be limited to specific chromosome(s) or region(s) by using `--region` option. If you are analyzing VCF from WGS of a mid/large cohort a realiable estimation can be obtained also analyzing single autosome.
