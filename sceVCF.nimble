# Package

version       = "0.1.0"
author        = "edoardo.giacopuzzi"
description   = "Fast calculation of metric useful to estimate sample contamination"
license       = "MIT"
srcDir        = "src"
bin           = @["bin/sceVCF"]
skipDirs      = @["test"]


# Dependencies

requires "nim >= 1.4.8", "hts >= 0.3.21", "argparse >= 3.0.0"
