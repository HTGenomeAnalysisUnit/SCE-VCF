# This is just an example to get you started. A typical binary package
# uses this file as the main entry point of the application.
import hts
import times
import strformat
import strutils
import sequtils
import streams
import tables
import math
#import sets
from os import fileExists
import scevcf/arg_parse
import scevcf/utils

const VERSION = "0.1.2"
const TSV_HEADER = "#SAMPLE\tHQ_HOM\tHQ_HOM_RATE\tHQ_HET\tHQ_HET_RATE\tCHARR\tMEAN_REF_AB_HOM_ALT\tHETEROZYGOSITY_RATE\tINCONSISTENT_AB_HET_RATE"
const AUTOSOMES = map(to_seq(1..22), proc(x: int): string = $x) 

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

proc update_values(sdata: var Table[string, Contamination_data], genos: Genotypes, ads: seq[int32], afs: seq[float32], gqs: seq[int32], dps: seq[int32], samples: seq[string], het_ab_limit: (float,float), minGQ: int, dp_limit: seq[int]) =
  for i in 0..samples.high:
    if dps[i] < dp_limit[0] or dps[i] > dp_limit[1]: continue
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
    dp_lims_opts = opts.dp_limit.split(",")
    dp_lims = @[parseInt(dp_lims_opts[0]), parseInt(dp_lims_opts[1])]
    refAF_lims_opts = opts.refaf_limit.split(",")
    refAF_lims = @[parseFloat(refAF_lims_opts[0]), parseFloat(refAF_lims_opts[1])]
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
    interval = 50000
    n = 0
    sample_data: Table[string, Contamination_data]
    n_multiallele = 0
    n_indels = 0
    n_large_af = 0
    n_no_aftag = 0

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

      #consider only PASS vars from autosomes
      if v.FILTER != "PASS" and v.FILTER != ".": continue
      var chrom = $v.CHROM
      if chrom.replace("chr","") notin AUTOSOMES: continue
      
      #consider only biallelic sites
      if len(v.ALT) > 1: 
        n_multiallele += 1
        continue

      #skip indels
      if v.REF.len > 1 or v.ALT[0].len > 1:
        n_indels += 1
        continue
      
      doAssert v.format.get(opts.ad_field, ads) == Status.OK
      doAssert v.format.get("GQ", gqs) == Status.OK
      doAssert v.format.get("DP", dps) == Status.OK 
      let genos = v.format.genotypes(gts)

      #skip variant if the selected AF field is not found
      if v.info.get(opts.af_field, afs) != Status.OK:
        n_no_aftag += 1
        continue

      #skip var where ref AF is too low or high
      #because this will results in outlier CHARR values
      var refAF: float
      if (opts.is_refAF):
        refAF = afs[0]
      else:
        refAF = 1-afs[0]
        
      #echo fmt"{$v} - AF: {afs[0]} - refAF {refAF} - AF lims low: {refAF_lims[0]} - AF lims high: {refAF_lims[1]}"
      if refAF < refAF_lims[0] or refAF > refAF_lims[1]: 
        n_large_af += 1
        continue 
      
      sample_data.update_values(genos, ads, afs, gqs, dps, vcf.samples, het_ab_limit, minGQ, dp_lims)
    
    log("INFO", fmt"{n} variants processed, {n_multiallele + n_large_af + n_no_aftag + n_indels} vars ignored: {n_multiallele} multiallelic, {n_indels} indels, {n_large_af} outside ref AF limits, {n_no_aftag} missing AF tag")
    close(vcf)
  
  log("INFO", "Computing contamination values")
  var written_samples = 0
  for sample_id, values in sample_data.pairs:
    written_samples += 1
    let
      charr_value = (values.charr / values.hom.hq.float).formatFloat(ffDecimal, 5)
      het_rate = (values.het.hq / (values.hom.hq + values.het.hq)).formatFloat(ffDecimal, 5)
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