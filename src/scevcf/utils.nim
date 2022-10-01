import times
import strutils
import strformat
import system
import math
import os

proc log* (level: string = "INFO", msg: string) {.discardable.} = 
    let t = now()
    let f = initTimeFormat("HH:mm:ss-fff")
    var time_string = format(t, f)
    let log_msg = fmt"[{time_string}] - {level}: {msg}"
    stderr.write(log_msg & "\n")

proc elapsed_time* (start_time: float): string =
    let interval = cpuTime() - start_time
    let s = floor(interval)
    let m = floor(((interval - s) * 1000))
    let time_interval = initDuration(seconds = int(s), milliseconds = int(m))
    result = $time_interval

proc progress_counter* (n:int, interval: var int, t0: var float): (bool, string) {.discardable.} =
    result = (false, "")

    case n
        #of 50000: interval = 25000
        of 150000: interval = 50000
        of 500000: interval = 100000
        of 1000000: interval = 500000
        else: discard
    
    if floorMod(n, interval) == 0:
            result = (true, fmt"{n} variants processed. Last batch: {interval} in {elapsed_time(t0)}")
            t0 = cpuTime()

proc read_list* (list: string, name: string = ""): seq[string] =
    if list == "":
        result = @[]
    else:
        if os.fileExists(list):
            log("INFO", fmt"Reading {name} list from file: {list}")
            for line in lines(list):
                result.add(line)
        else:
            log("INFO", fmt"Reading {name} list from comma-separated string")
            result = list.split(",")

        log("INFO", fmt"Imported {result.len} {name}s")