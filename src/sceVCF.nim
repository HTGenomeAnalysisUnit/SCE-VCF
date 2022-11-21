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
import scevcf/arg_parse
import scevcf/utils

const VERSION = "0.1.1"
const TSV_HEADER = "#SAMPLE\tHQ_HOM\tHQ_HOM_RATE\tHQ_HET\tHQ_HET_RATE\tCHARR\tMEAN_REF_AB_HOM_ALT\tHET_RATE\tINCONSISTENT_AB_HET_RATE"

type Contamination_data = object
  charr: float
  ref_ab: float
  het: tuple[n: int, hq: int, bad: int]
  hom: tuple[n: int, hq: int]

iterator readvar(v: VCF, regions: seq[string]): Variant =
  if regions.len == 0:
    for variant in v: yield variant
  else:
    for r in regions:
      for variant in v.query(r): yield variant

proc update_values(sdata: var Table[string, Contamination_data], genos: Genotypes, ads: seq[int32], afs: seq[float32], gqs: seq[int32], dps: seq[int32], samples: seq[string], het_ab_limit: (float,float), minGQ: int, minDP: int) =
  for i in 0..samples.high:
    if dps[i] < minDP: continue
    let 
      ref_ad = ads[i*2]
      alt_ad = ads[i*2+1]
      tot_ad = ref_ad + alt_ad
    var x = sdata.getOrDefault(samples[i])

    case genos[i].alts:
      of 2: #compute charr for hom alt vars
        x.hom.n += 1
        if gqs[i] < minGQ: continue
        let 
          ref_af = 1 - (if afs[0] < 0: 0.0 else: afs[0])
          ref_ab = ref_ad / tot_ad
          charr = ref_ab * (1/ref_af)
        if charr.classify != fcNan and charr.classify != fcInf:
          #echo fmt"{$genos[i].alts} - TOT_AD: {tot_ad} - REF_AD: {ref_ad} - REF_AB: {ref_ab} - REF_AF: {ref_af} - CHARR: {charr}"
          x.charr += charr
        x.hom.hq += 1
        if ref_ab.classify != fcNan and ref_ab.classify != fcInf:
          x.ref_ab += ref_ab
      of 1: #check 
        x.het.n += 1
        if gqs[i] >= minGQ: x.het.hq += 1
        let het_ab = alt_ad / tot_ad
        if het_ab < het_ab_limit[0] or het_ab > het_ab_limit[1]:
          x.het.bad += 1
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
    n_large_af = 0
    interval = 50000
    n = 0
    sample_data: Table[string, Contamination_data]

  var
    ads: seq[int32]
    gts: seq[int32]
    afs: seq[float32]
    gqs: seq[int32]
    dps: seq[int32]

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

    for v in vcf.readvar(regions):
      n = n + 1
      var (dolog, log_msg) = progress_counter(n, interval, t0)
      if dolog: log("INFO", log_msg)

      var chrom = $v.CHROM
      if chrom.replace("chr","") in @["X","Y","M","MT"]: continue
      
      if len(v.ALT) > 1: #consider only biallelic sites
        n_multiallele += 1
        continue
      
      doAssert v.format.get(opts.ad_field, ads) == Status.OK
      doAssert v.format.get("GQ", gqs) == Status.OK
      doAssert v.format.get("DP", dps) == Status.OK
      doAssert v.info.get(opts.af_field, afs) == Status.OK
      let genos = v.format.genotypes(gts)

      if afs[0] > 0.95: 
        n_large_af += 1
        continue #skip vars where ref is too rare
      
      sample_data.update_values(genos, ads, afs, gqs, dps, vcf.samples, het_ab_limit, minGQ, minDP)
    
    log("INFO", fmt"{n} variants processed, {n_multiallele + n_large_af} vars ignored: {n_multiallele} multiallelic, {n_large_af} AF > 0.95")
    close(vcf)
  
  log("INFO", "Computing contamination values")
  var written_samples = 0
  for sample_id, values in sample_data.pairs:
    written_samples += 1
    let
      charr_value = (values.charr / values.hom.hq.float).formatFloat(ffDecimal, 5)
      het_rate = (values.het.hq / values.hom.hq).formatFloat(ffDecimal, 5)
      bad_het_rate = (values.het.bad / values.het.n).formatFloat(ffDecimal, 5)
      mean_ref_ab = (values.ref_ab / values.hom.hq.float).formatFloat(ffDecimal, 5)
      hom_hq_rate = (values.hom.hq / values.hom.n).formatFloat(ffDecimal, 5)
      het_hq_rate = (values.het.hq / values.het.n).formatFloat(ffDecimal, 5)

    let result_line = &"{sample_id}\t{values.hom.hq}\t{hom_hq_rate}\t{values.het.hq}\t{het_hq_rate}\t{charr_value}\t{mean_ref_ab}\t{het_rate}\t{bad_het_rate}"
    if write_to_file:
      out_tsv.writeLine(result_line)
    else:
      stdout.writeLine(result_line)

  if write_to_file: close(out_tsv)
  log("INFO", fmt"Computed contamination values for {written_samples} sample(s) in {elapsed_time(start_time)}")

when isMainModule:
  main()