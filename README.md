# Sample Contamination Estimate from VCF (SCE-VCF)

SCE-VCF is a simple tool to evalute sample contamination from sequencing experiments. It uses WGS/WES VCF files to compute few metrics useful to detect potential sample contamination:

- CHARR: ref read fraction in hom alt vars as proposed by Wenham Lu, Broad Institute
- Mean ref allele AB observed in hom alt genotypes
- heterozygosity ratio (N HQ het / (N HQ het + N HQ hom))
- rate on inconsistent het calls (het call with AB outside threshold)

The ideal scenario is a cohort analysis so one can compare value distributions to identify outliers. Even a small cohort of 50 samples works fine in our test.

Example distributions for contamination values based on 1000G 30X WGS samples and variants identified using either GATK or DeepVariant are provided in the example folder.

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
  -i, --is_refAF             AF in the af_field is the REF AF
  -t, --refaf_limit=REFAF_LIMIT
                             Limits of REF AF. Comma-separated lower and upper limit (default: 0.1,0.9)
  -l, --het_ab_limit=HET_AB_LIMIT
                             Comma separated min and max allele balance accepted for het calls (default: 0.25,0.75)
  -q, --min_GQ=MIN_GQ        Min GQ for hom var to be included in charr computation (default: 20)
  -d, --dp_limit=DP_LIMIT    Limits of genotype DP. Comma-separated lower and upper limit (default: 20,100)
  -s, --samples=SAMPLES      Restrict analysis to the given samples. Comma separated list or file wiht one sample per line
  -r, --region=REGION        Specify genomic region for processing. Format is chr[:start-end]. Comma-separated list of regions or file with 1 region per line.
```

## How it works

1. Multi-allelic variants and variants with REF AF outside the configured limits (`refaf_limit`) are skipped.
2. Any genotype with DP outside the configure DP limits `dp_limit` is skipped.
3. Genotypes with GQ less than `min_GQ` are used only to compute inconsistent het rate, but discarded when computed het rate and CHARR.
4. CHARR is computed as `ADref / (AFref * ADtot)` considering only high-quality hom alt genotypes.

## Output columns

| Column | Description |
|--------|-------------|
| SAMPLE | Sample ID as read from input VCF file |
| HQ_HOM | Number of high quality homozygous sites used in CHARR computation |
| HQ_HOM_RATE | Fraction of HQ hom sites across all hom sites with enough coverage |
| HQ_HET | Number of high quality heterozygous sites used to compute het rate |
| HQ_HET_RATE | Fraction of HQ het sites across all het sites with enough coverage |
| CHARR | CHARR value computed as explained above using only HQ hom alt sites |
| MEAN_REF_AB_HOM_ALT | Mean ref allele frequency observed in HQ hom alt sites |
| HETEROZYGOSITY_RATE | Heterozygoisty ratio (Nhet / (Nhet + Nhom)) computed using only HQ genotypes |
| INCONSISTENT_AB_HET_RATE | Fraction of het sites with AB outside threashold across all het sites with enough coverage (this includes all het sites with DP > threashold, disregarding GQ) |

## Interpretation of results

If you have analyzed a large enough cohort or you have a cohort processed with a similar pipeline for variant-calling, we suggest to consider the distribution of CHARR, HETEROZYGOSITY_RATE, and INCONSISTENT_AB_HET_RATE values across all individuals and mark as possibly contaminated the outliers with values larger than mean + 3SD.

Based on our test, a CHARR value of 0.02-0.03 can be considered a warning signal and a value larger than 0.03 is strongly suggestive of contamination. However, we suggest to jointly evaluate all the 3 metrics when determining contaminated sample. Usually contaminated samples show also high value for HETEROZYGOSITY_RATE (>= 0.7) and/or INCONSISTENT_AB_HET_RATE (>= 0.05).

## Suggested workflow

If you are processing a cohort VCF that contains enough samples to reliably estimate the cohort allele frequency and this value is annoteated in the INFO column (usually AF), then you can use your file directly to estimate contamination in each sample.

If you are working with a single sample VCF or a small cohort we suggest to proceed as follows:

1. Select only PASS bi-allelic SNPs using `bcftools view -f PASS,. -m2 -M2 -v snps`
2. Annotate the resulting file with AF from large reference population. For example, you can rapidly annotate global AF from gnomAD using echtvar (see [the official repo](https://github.com/brentp/echtvar))
3. Perform contamination estimation

**NB.** Decomposing the file with `bcftools norm` may lead to unexpected results due to variants derived from multi-allelic split being included in the computation. We suggest to instead select only biallelic sites before normalize (see point 1 above).

## Running time

Running time on the full 1000G cohort VCF containing 2504 samples and about 120M variants is ~20h.

## Notes

- **AF data:** The tool needs an annotated AF field in the input VCF to compute REF AF and the CHARR values. Ideally, this should be AF computed from the samples in the VCF. Even a small cohort of 50 samples works fine in our test. When analyzing a single sample VCF or a small cohort, we suggest to annotate your file with AF from large external populations (like gnomAD) and then pass the relevant INFO field to `--af-field`.

- **Multiple input VCFs:** The tool can read multiple VCFs and statistics are then aggregated by sample IDs if there are matching sample IDs across the input VCFs. In this way it's easy to analyze cohort VCFs splitted by chromosome for example.

- **Subset samples:** User can specify a subset of samples to analyze in the cohort VCF, this can be conveniently used also to parallelize the analysis for large cohorts. Note that some time is required anyway to parse sample initially, so run time will increase as more samples are present in the input VCF

- **Subset regions** Analysis can be limited to specific chromosome(s) or region(s) by using `--region` option. If you are analyzing VCF from WGS of a mid/large cohort a realiable estimation can be obtained also analyzing chr1 only for example.
