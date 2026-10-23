"""Exercise autosome filtering against a compiled sceVCF executable."""

import csv
import io
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


BINARY = str(Path(sys.argv.pop(1)).resolve()) if len(sys.argv) > 1 else "sceVCF"


class ContigFilterTest(unittest.TestCase):
    def run_vcf(self, contigs):
        header = [
            "##fileformat=VCFv4.2",
            '##FILTER=<ID=PASS,Description="All filters passed">',
            '##INFO=<ID=AF,Number=A,Type=Float,Description="ALT frequency">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
            '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Quality">',
        ]
        header += [f"##contig=<ID={chrom}>" for chrom in dict.fromkeys(contigs)]
        header += ["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample"]
        records = []
        for chrom in contigs:
            for pos, genotype in [(1, "1/1:2,38:40:99"), (2, "0/1:20,20:40:99")]:
                records.append(
                    f"{chrom}\t{pos}\t.\tA\tC\t.\tPASS\tAF=0.5\tGT:AD:DP:GQ\t{genotype}"
                )
        with tempfile.TemporaryDirectory(prefix="scevcf-contigs-") as directory:
            vcf = Path(directory) / "input.vcf"
            vcf.write_text("\n".join(header + records) + "\n")
            return subprocess.run(
                [BINARY, str(vcf)], capture_output=True, text=True, check=False
            )

    def check_autosomes(self, contigs, kept):
        result = self.run_vcf(contigs)
        self.assertEqual(result.returncode, 0, result.stderr)
        skipped = 2 * (len(contigs) - kept)
        lines = result.stdout.splitlines()
        self.assertEqual(
            lines[0],
            f"## {2 * kept} variants considered for analysis, "
            f"{skipped}/{2 * len(contigs)} vars ignored",
        )
        rows = list(csv.DictReader(io.StringIO("\n".join(lines[1:])), delimiter="\t"))
        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertEqual(row["#SAMPLE"], "sample")
        for field in ("HOM_TOTAL", "HOM_COVERED", "HQ_HOM", "HET_TOTAL", "HET_COVERED", "HQ_HET"):
            self.assertEqual(int(row[field]), kept, field)
        self.assertAlmostEqual(float(row["CHARR"]), 0.1)
        self.assertAlmostEqual(float(row["MEAN_REF_AB_HOM_ALT"]), 0.05)
        self.assertAlmostEqual(float(row["HETEROZYGOSITY_RATE"]), 0.5)
        self.assertAlmostEqual(float(row["INCONSISTENT_AB_HET_RATE"]), 0.0)

    def test_all_autosomes(self):
        for prefix in ("", "chr"):
            with self.subTest(prefix=prefix):
                self.check_autosomes([f"{prefix}{i}" for i in range(1, 23)], 22)

    def test_consecutive_excluded_contigs(self):
        self.check_autosomes(
            ["chrX", "chrY", "chrM", "1", "2", "chrUn", "GL000220.1", "chr21", "chr22", "MT"],
            4,
        )

    def test_only_excluded_contigs(self):
        result = self.run_vcf(["X", "Y", "MT", "chrX", "chrY", "chrM", "chrUn"])
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("14 not in autosomes", result.stderr)
        self.assertIn("No variants remaining for analysis", result.stderr)
        self.assertEqual(result.stdout, "")


if __name__ == "__main__":
    unittest.main()
