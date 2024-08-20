# Sample Contamination Estimate from VCF (SCE-VCF)

SCE-VCF is a simple tool to evalute sample contamination from sequencing experiments. It uses WGS/WES VCF files to compute few metrics useful to detect potential sample contamination:

- CHARR: ref read fraction in hom alt vars as proposed by Wenham Lu, Broad Institute
- Mean ref allele AB observed in hom alt genotypes
- heterozygosity ratio (N HQ het / (N HQ het + N HQ hom))
- rate on inconsistent het calls (het call with AB outside threshold)

The ideal scenario is a cohort analysis so one can compare value distributions to identify outliers. Even a small cohort of 50 samples works fine in our test.

Example distributions for contamination values based on 1000G 30X WGS samples and variants identified using either GATK or DeepVariant are provided in the example folder.

## Installation

On most Linux based systems, you can just download the pre-compiled executable from the latest release and it should work out of the box.

### Compile from source

Compiling the executable with nim has the following requirements:

- nim >= 1.4.8
- hts >= 0.3.21
- argparse >= 3.0.0

You can usually install them with `nimble install hts argparse`.

Then you can clone the repository and compile the executable with `nim c -d:release src/sceVCF.nim`.

### Generate a statically linked executable

If you have singularity installed, you can generate a statically linked executable with the following command:

```bash
bash nim_compile.sh
```

See [here](See: https://github.com/brentp/hts-nim#static-builds) for more details on how this works. 

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

## Suggested workflow

If you are processing a cohort VCF that contains enough samples to reliably estimate the cohort allele frequency and this value is annoteated in the INFO column (usually AF), then you can use your file directly to estimate contamination in each sample.

**NB.** Decomposing the file with `bcftools norm` may lead to unexpected results due to variants derived from multi-allelic split being included in the computation. We suggest to instead select only biallelic SNPs (for example `bcftools view -f PASS,. -m2 -M2 -v snps`) and then compute contamination metrics on this subset.

If you are working with a single sample VCF or a small cohort we suggest to proceed as follows:

1. Select only PASS bi-allelic SNPs using `bcftools view -f PASS,. -m2 -M2 -v snps`
2. Annotate the resulting file with AF from large reference population. For example, you can rapidly annotate global AF from gnomAD using echtvar (see [the official repo](https://github.com/brentp/echtvar))
3. Perform contamination estimation

## Output columns

| Column | Description |
|--------|-------------|
| SAMPLE | Sample ID as read from input VCF file |
| HOM_TOTAL | Total number of ALT homozygous sites observed for this sample before DP/GQ filters |
| HOM_COVERED | Number of ALT homozygous sites within the configured DP limits |
| HQ_HOM | Number of high quality ALT homozygous sites used in CHARR computation |
| HQ_HOM_RATE | Fraction of HQ hom sites across all hom sites with enough coverage (HQ_HOM / HOM_COVERED) |
| HET_TOTAL | Total number of heterozygous sites observed for this sample before DP/GQ filters |
| HET_COVERED | Number of heterozygous sites within the configured DP limits |
| HQ_HET | Number of high quality heterozygous sites used to compute het rate |
| HQ_HET_RATE | Fraction of HQ het sites across all het sites with enough coverage (HQ_HET / HET_COVERED) |
| CHARR | CHARR value computed as explained above using only HQ hom alt sites |
| MEAN_REF_AB_HOM_ALT | Mean ref allele frequency observed in HQ hom alt sites |
| HETEROZYGOSITY_RATE | Heterozygoisty ratio (Nhet / (Nhet + Nhom)) computed using only HQ genotypes |
| INCONSISTENT_AB_HET_RATE | Fraction of het sites with AB outside threshold across all het sites with enough coverage (this includes all het sites with DP within threshold, disregarding GQ) |

**NB.** All values in result table are based only on variant passing the global variants filter (PASS, biallelic, SNV only, REF AF within the configured threshold). The number of variants considered for the analysis is stated in the header.

## Interpretation of results

If you have analyzed a large enough cohort or you have a cohort processed with a similar pipeline for variant-calling, we suggest to identify contaminated samples considering the distribution of CHARR, HETEROZYGOSITY_RATE, and INCONSISTENT_AB_HET_RATE values across all individuals. For example, one can mark as possibly contaminated the outliers with values larger than mean + 2 or 3 SD.

When evaluating a single sample, a CHARR value of 0.02-0.03 may indicate a subtle contamination (1-5%) and a value larger than 0.03 is strongly suggestive of contamination (contamination > 5%). However, we suggest to jointly evaluate all the 3 metrics when determining contaminated sample. When a sample is heavily contaminated (> 20%), CHARR metric may look good, but you can observe high value for HETEROZYGOSITY_RATE (>= 0.7) and INCONSISTENT_AB_HET_RATE (>= 0.1).

In the end we suggest the following thresholds for contamined samples:

- WARNING: CHARR > 0.02 | INCONSISTENT_AB_HET_RATE > 0.1
- FAIL: CHARR > 0.03 | INCONSISTENT_AB_HET_RATE > 0.15

This approach classify samples similarly to traditional freemix approach, but it is more robust to contamination of different magnitude and it is more sensitive to subtle contamination.

### Simulation results

To simulate contamination we generated mixed samples using 2 samples from the GIAB reference datasets. Briefly, we generated a 30X BAM file of NA24631 GIAB sample and then mixed this with various proportion of reads from NA24385. See results below and plots in the example folder.

|SAMPLE         |HQ_HOM |HQ_HOM_RATE|HQ_HET |HQ_HET_RATE|CHARR  |MEAN_REF_AB_HOM_ALT|HETEROZYGOSITY_RATE|INCONSISTENT_AB_HET_RATE|
|---------------|-------|-----------|-------|-----------|-------|-------------------|-------------------|------------------------|
|original_sample|1101063|1.00000    |1516921|0.96191    |0.01144|0.00371            |0.57942            |0.00999                 |
|cont_0.99_0.01 |1092683|1.00000    |1515890|0.96113    |0.01669|0.00598            |0.58112            |0.01068                 |
|cont_0.98_0.02 |1079498|1.00000    |1514989|0.95938    |0.02134|0.00794            |0.58393            |0.01249                 |
|cont_0.97_0.03 |1062034|1.00000    |1514340|0.95648    |0.02543|0.00963            |0.58778            |0.01580                 |
|cont_0.96_0.04 |1040781|1.00000    |1514274|0.95240    |0.02899|0.01107            |0.59266            |0.02079                 |
|cont_0.95_0.05 |1027437|1.00000    |1512849|0.94977    |0.03044|0.01164            |0.59554            |0.02412                 |
|cont_0.85_0.15 |741203 |1.00000    |1578589|0.86162    |0.04422|0.01578            |0.68049            |0.15782                 |
|cont_0.9_0.1   |874280 |1.00000    |1532181|0.90654    |0.04179|0.01559            |0.63669            |0.08385                 |
|cont_0.8_0.2   |635342 |1.00000    |1652251|0.82522    |0.04192|0.01431            |0.72227            |0.22597                 |
|cont_0.7_0.3   |505642 |1.00000    |1832859|0.78729    |0.03092|0.00972            |0.78378            |0.29380                 |
|cont_0.5_0.5   |409987 |1.00000    |2103835|0.78209    |0.01383|0.00371            |0.83691            |0.29484                 |

## Running time

Running time on a single sample WGS VCF containing about 4M variants is ~25 secs. Analysis of the full 1000G cohort VCF containing 2504 samples and about 120M variants takes ~20h.

## Notes

- **AF data:** The tool needs an annotated AF field in the input VCF to compute REF AF and the CHARR values. Ideally, this should be AF computed from the samples in the VCF. Even a small cohort of 50 samples works fine in our test. When analyzing a single sample VCF or a small cohort, we suggest to annotate your file with AF from large external populations (like gnomAD) and then pass the relevant INFO field to `--af-field`.

- **Multiple input VCFs:** The tool can read multiple VCFs and statistics are then aggregated by sample IDs if there are matching sample IDs across the input VCFs. In this way it's easy to analyze cohort VCFs splitted by chromosome for example.

- **Subset samples:** User can specify a subset of samples to analyze in the input VCF. Note that some time is required anyway to parse samples initially, so run time will increase as more samples are present in the input VCF

- **Subset regions** Analysis can be limited to specific chromosome(s) or region(s) by using `--region` option.
