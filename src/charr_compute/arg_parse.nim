import argparse
import strformat
from ./utils import log

var p = newParser("CSQ Selector"):
    help("Select gene consequences from snpEff, VEP or bcftools annotations in VCF file")
    arg("vcf", nargs = -1, help="input VCF/BCF file(s). Glob pattern allowed")
    option("-o", "--out", help="Output file (TSV). If not provided output to stdout")
    option("-f", "--ad_field", help="FORMAT field containing the AD values", default = some("AD"))
    option("-a", "--af_field", help="INFO field containing the allele frequency", default = some("AF"))
    option("-l", "--het_ab_limit", help="Comma separated min and max allele balance accepted for het calls", default = some("0.25,0.75"))
    option("-q", "--min_GQ", help="Min GQ for hom var to be included in charr computation", default = some("20"))
    option("-d", "--min_DP", help="Min DP for var to be included in computation", default = some("20"))
    option("-s", "--samples", help="Restrict analysis to the given samples. Comma separated list or file wiht one sample per line")
    option("-r", "--region", help="Specify genomic region for processing. Format is chr[:start-end]. Comma-separated list of regions or file with 1 region per line.")                

proc parseCmdLine*(): ref =
    try:
        result = p.parse() 
    except ShortCircuit as e:
        if e.flag == "argparse_help":
            echo p.help
            quit QuitSuccess
    except UsageError:
        stderr.writeLine getCurrentExceptionMsg() 
        echo "Use --help for usage information"
        quit QuitSuccess

proc logArgs*(opts: ref) {.discardable.} =
    log("ARG", fmt"Input VCF: {opts.vcf}")
    let out_stream = (if opts.out != "": opts.out else: "stdout")
    log("ARG", fmt"Output: {out_stream}")
    log("ARG", fmt"AD field: {opts.ad_field}")
    log("ARG", fmt"AF field: {opts.af_field}")
    log("ARG", fmt"Het AB limit: {opts.het_ab_limit}")
    log("ARG", fmt"Min GQ: {opts.min_GQ}")
    log("ARG", fmt"Min DP: {opts.min_DP}")
    if opts.samples != "":
        log("ARG", fmt"Samples {opts.samples}")
