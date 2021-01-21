#!/usr/bin/env python

import pysam
import sys
import csv
import os
import pkg_resources
from pathlib import Path


class Variant:
    def __init__(self, contig, position, reference, alt):
        self.contig = contig
        self.position = position
        self.reference = reference
        self.alt = alt
        self.name = None

    def key(self):
        return ",".join([str(x) for x in [self.contig, self.position, self.reference, self.alt]])


def load_vcf(filename):
    variants = list()
    f = pysam.VariantFile(filename, 'r')
    for record in f:
        if len(record.alts) > 1:
            sys.stderr.write("Multi-allelic VCF not supported\n")
            sys.exit(1)

        v = Variant(record.chrom, record.pos, record.ref, record.alts[0])
        if "Name" in record.info:
            v.name = record.info["Name"]
        variants.append(v)

    return variants


# noinspection PyBroadException
def load_ivar_variants(filename):
    variants = list()
    try:
        with open(filename, 'r') as ifh:
            reader = csv.DictReader(ifh, delimiter='\t')
            for record in reader:
                ref = record["REF"]
                alt = record["ALT"]
                if alt[0] == "-":
                    ref += alt[1:]
                    alt = ref[0]
                elif alt[0] == "+":
                    alt = ref + alt[1:]

                variants.append(Variant(record["REGION"], record["POS"], ref, alt))
    except:
        pass
    return variants


def get_from_directory(directory):
    def _glob(pattern):
        return [path for path in Path(directory).rglob(pattern)]

    def _inner():
        files = _glob("*pass.vcf") + _glob("*pass.vcf.gz") + _glob("*variants.tsv")
        for f in files:
            yield str(f)
    return _inner


def get_from_stdin():
    for line in sys.stdin:
        yield line.rstrip()


def main():
    mutation_sets = pkg_resources.resource_listdir(__name__, 'watchlists')
    mutation_sets = [Path(mutation_set).stem for mutation_set in mutation_sets]

    import argparse
    parser = argparse.ArgumentParser(description='Report samples containing a variant in the watchlist')
    parser.add_argument('-m', '--mutation-set', required=True,
                        choices=mutation_sets,
                        help='Mutation set to screen variants against')
    parser.add_argument('-o', '--output',
                        type=argparse.FileType('w'),
                        default=sys.stdout,
                        help='Write output to file, defaults to stdout')
    parser.add_argument('input', nargs='?',
                        help='Root of directories holding variant files, blank for stdin')

    args = parser.parse_args()
    out = csv.writer(args.output, delimiter='\t')

    mutation_set_path = Path("watchlists") / Path(args.mutation_set + ".vcf")
    mutation_set = pkg_resources.resource_filename(__name__, str(mutation_set_path))

    watch_variants = load_vcf(mutation_set)
    watch_dict = {v.key(): v.name for v in watch_variants}

    gen_func = get_from_stdin
    if args.input:
        gen_func = get_from_directory(args.input)

    out.writerow(["sample", "mutation", "contig", "position", "reference", "alt"])
    for f in gen_func():
        if f.find("variants.tsv") >= 0:
            variants = load_ivar_variants(f)
        else:
            variants = load_vcf(f)

        for v in variants:
            if v.key() in watch_dict:
                name = watch_dict[v.key()]
                if name is None:
                    name = "not annotated"
                out.writerow([os.path.basename(f), name, v.contig, str(v.position), v.reference, v.alt])


if __name__ == "__main__":
    main()
