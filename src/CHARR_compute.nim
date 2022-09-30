# This is just an example to get you started. A typical binary package
# uses this file as the main entry point of the application.
import hts
import times
import strformat
import strutils
import streams
import tables
import math
#import sets
from os import fileExists
import charr_compute/arg_parse
import charr_compute/utils

const VERSION = "0.1.1"
const TSV_HEADER = "#SAMPLE\tCHARR\tVALID_HOM_FOR_CHARR\tMEAN_REF_AB_HOM_ALT\tN_HETS\tHET_RATE\tINCONSISTENT_HET_RATE"

type Contamination_data = object
  charrs: float
  abnormal_charrs: int
  ref_ab: float
  n_het: int
  n_hom: int
  abnormal_hets: int

iterator readvar(v: VCF, regions: seq[string]): Variant =
  if regions.len == 0:
    for variant in v: yield variant
  else:
    for r in regions:
      for variant in v.query(r): yield variant

proc update_values(sdata: var Table[string, Contamination_data], v: Variant, ad: string, af: string, samples: seq[string], het_ab_limit: (float,float), minGQ: int, minDP: int) =
  var
    ads: seq[int32]
    gts: seq[int32]
    afs: seq[float32]
    gqs: seq[int32]
    dps: seq[int32]

  doAssert v.format.get(ad, ads) == Status.OK
  doAssert v.format.get("GQ", gqs) == Status.OK
  doAssert v.format.get("DP", dps) == Status.OK
  doAssert v.info.get(af, afs) == Status.OK
  let genos = v.format.genotypes(gts)
   
  for i in 0..samples.high:
    if dps[i] < minDP: continue
    let sample_ads = @[ads[i*2],ads[i*2+1]]
    var x = sdata.getOrDefault(samples[i])

    case genos[i].alts:
      of 2: #compute charr for hom alt vars
        if gqs[i] < minGQ: continue
        let ref_af = 1 - (if afs[0] < 0: 0.0 else: afs[0])
        let ref_ab = sample_ads[0] / sum(sample_ads)
        let charr = ref_ab.float * (1/ref_af)
        #echo fmt"{$v} - {$genos[i].alts} - ADS: {sample_ads} - AF: {afs[0]} - REF_AF: {ref_af} - CHARR: {charr}"
        if charr.classify == fcNan or charr.classify == fcInf:
          x.abnormal_charrs += 1
        else:
          x.charrs += charr
        x.n_hom += 1
        x.ref_ab += ref_ab
      of 1: #check 
        x.n_het += 1
        let het_ab = sample_ads[1] / sum(sample_ads)
        if het_ab < het_ab_limit[0] or het_ab > het_ab_limit[1]:
          x.abnormal_hets += 1
      else:
        discard

    sdata[samples[i]] = x 
    
proc main* () =
  log("INFO", fmt"SCE-VCF v{VERSION}")

  #Parse command line arguments
  var opts = parseCmdLine()
  opts.logArgs()

  let 
    regions = read_list(opts.region, "region")
    minGQ = parseInt(opts.min_GQ)
    minDP = parseInt(opts.min_DP)
  var samples = read_list(opts.samples, "samples")
  var het_ab_limit: (float, float) 
  try:
    het_ab_limit = (
      parseFloat(opts.het_ab_limit.split(",")[0]),
      parseFloat(opts.het_ab_limit.split(",")[1])
    )
  except ValueError:
    log("FATAL", fmt"Cannot parse float values from het_ab_limit '{opts.het_ab_limit}'")

  #Open output file and write header
  var write_to_file = false
  if opts.out != "": write_to_file = true

  var out_tsv:FileStream
  if write_to_file:
    out_tsv = newFileStream(opts.out, fmWrite)
    if isNil(out_tsv):
        log("FATAL", "Could not open output file " & opts.out)
        quit "", QuitFailure
    else:
        out_tsv.writeLine(TSV_HEADER)
  else:
    stdout.writeLine(TSV_HEADER)

  for f in opts.vcf:
    if not fileExists(f):
      log("FATAL", "Could not open input VCF file " & f)
      quit "", QuitFailure

  #Process input files
  log("INFO", "Variant processing started")
  var
    start_time = cpuTime()
    t0 = cpuTime()
    n_multiallele = 0
    interval = 50000
    n = 0
    sample_data: Table[string, Contamination_data]

  for f in opts.vcf:
    log("INFO", fmt"Processing file {f}")
    var vcf:VCF
    discard open(vcf, f, samples=samples)
      
    var desc: string
    try:
      desc = vcf.header.get(opts.af_field, BCF_HEADER_TYPE.BCF_HL_INFO)["Description"]
    except:
      log("FATAL", fmt"Didn't find {opts.af_field} INFO in header in {vcf.fname}")

    try:
      desc = vcf.header.get(opts.ad_field, BCF_HEADER_TYPE.BCF_HL_FMT)["Description"]
    except:
      log("FATAL", fmt"Didn't find {opts.ad_field} FORMAT in header in {vcf.fname}")

    var empty_val: Contamination_data
    empty_val.n_het = 0
    empty_val.n_hom = 0
    empty_val.charrs = 0
    empty_val.abnormal_hets = 0
    for s in vcf.samples:
      discard sample_data.hasKeyOrPut(s, empty_val)

    for v in vcf.readvar(regions):
      n = n + 1
      var (dolog, log_msg) = progress_counter(n, interval, t0)
      if dolog: log("INFO", log_msg)

      var chrom = $v.CHROM
      if chrom.replace("chr","") in @["X","Y","M","MT"]: continue  
      if len(v.ALT) > 1: #consider only biallelic sites
        n_multiallele += 1
        continue

      sample_data.update_values(v, opts.ad_field, opts.af_field, vcf.samples, het_ab_limit, minGQ, minDP)
    log("INFO", fmt"{n} variants processed, {n_multiallele} multiallelic vars ignored")
    close(vcf)
  

  log("INFO", "Computing contamination values")
  var written_samples = 0
  for sample_id, values in sample_data.pairs:
    written_samples += 1
    let
      valid_charrs =  values.n_hom - values.abnormal_charrs
      charr_value = values.charrs / valid_charrs.float
      het_rate = values.n_het / values.n_hom
      bad_het_rate = values.abnormal_hets / values.n_het
      mean_ref_ab = values.ref_ab / values.n_hom.float

    let result_line = &"{sample_id}\t{charr_value}\t{valid_charrs}\t{mean_ref_ab}\t{values.n_het}\t{het_rate}\t{bad_het_rate}"
    if write_to_file:
      out_tsv.writeLine(result_line)
    else:
      stdout.writeLine(result_line)

  if write_to_file: close(out_tsv)
  log("INFO", fmt"Computed contamination values for {written_samples} sample(s) in {elapsed_time(start_time)}")

when isMainModule:
  main()
