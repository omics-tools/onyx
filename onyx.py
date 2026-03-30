#!/usr/bin/env python3

import argparse
import json
import re
import shutil
import subprocess
import sys
import tempfile
import tarfile
import hashlib
from pathlib import Path
import requests
from tqdm import tqdm, trange

VERSION = "0.0.1"

DB_LIST_URL = "https://raw.githubusercontent.com/omics-tools/onyx/main/databases.json"

class CustomFormatter(
    argparse.ArgumentDefaultsHelpFormatter,
    argparse.RawDescriptionHelpFormatter
):
    pass

#Print ONYX title

def print_title():
    return f"""

     ██████╗ ███╗   ██╗██╗   ██╗██╗  ██╗
    ██╔═══██╗████╗  ██║╚██╗ ██╔╝╚██╗██╔╝
    ██║   ██║██╔██╗ ██║ ╚████╔╝  ╚███╔╝
    ██║   ██║██║╚██╗██║  ╚██╔╝   ██╔██╗
    ╚██████╔╝██║ ╚████║   ██║   ██╔╝ ██╗
     ╚═════╝ ╚═╝  ╚═══╝   ╚═╝   ╚═╝  ╚═╝  ver.{VERSION}

An alignment-free biological sex inference from high-throughput sequencing data
Developer: Koji Ishiya
"""


#External command info
# Run external command
# Suppresses stdout/stderr unless verbose=True

def run(cmd, verbose=False):
    if verbose:
        print("[CMD]", " ".join(map(str, cmd)))
        subprocess.run(cmd, check=True)
    else:
        subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )


#ONYX environment check

_ENV_DONE = False

def ensure_env(tools):
    global _ENV_DONE
    if _ENV_DONE:
        return
    for tool in tools:
        if shutil.which(tool) is None:
            print(f"ERROR: {tool} not found in PATH", file=sys.stderr)
            sys.exit(1)
    _ENV_DONE = True


#Database metadata utility

def fetch_db_list():
    try:
        r = requests.get(DB_LIST_URL)
        r.raise_for_status()
        return r.json()
    except Exception as e:
        raise RuntimeError(f"Failed to fetch database list from {DB_LIST_URL}\n{e}")


def list_databases():
    db_meta = fetch_db_list()

    print(print_title())
    print("Available ONYX preset databases:\n")

    for name, meta in db_meta.items():
        species = meta.get("species", "unknown")
        ref = meta.get("reference", "unknown")
        print(f"{name:<10} {species} ({ref})")


#Download preset databases

def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def download_db(args):

    if args.list:
        list_databases()
        return

    db_meta = fetch_db_list()

    if args.name not in db_meta:
        raise RuntimeError(
            f"Unknown database: {args.name}\nAvailable: {', '.join(db_meta.keys())}"
        )

    meta = db_meta[args.name]
    url = meta["url"]

    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    tar_path = outdir / Path(url).name

    print(print_title())
    print(f"[INFO] Downloading database: {args.name}")

    existing = tar_path.stat().st_size if tar_path.exists() else 0
    headers = {"Range": f"bytes={existing}-"} if existing > 0 else {}

    response = requests.get(url, stream=True, headers=headers)
    response.raise_for_status()

    downloaded_size = int(response.headers.get("content-length", 0))
    total_size = existing + downloaded_size
    mode = "ab" if existing > 0 else "wb"

    with open(tar_path, mode) as f, tqdm(
        total=total_size,
        initial=existing,
        unit="B",
        unit_scale=True,
        unit_divisor=1024,
        desc="Downloading",
    ) as bar:
        for chunk in response.iter_content(chunk_size=8192):
            if chunk:
                f.write(chunk)
                bar.update(len(chunk))

    print("[INFO] Verifying checksum")
    observed_sha256 = sha256_file(tar_path)

    if meta["sha256"] != "dummy" and observed_sha256 != meta["sha256"]:
        raise RuntimeError(
            f"Checksum mismatch for {tar_path.name}\n"
            f"expected: {meta['sha256']}\n"
            f"observed: {observed_sha256}"
        )

    print("[INFO] Extracting database")
    with tarfile.open(tar_path, "r:gz") as tar:
        tar.extractall(outdir, filter="data")

    tar_path.unlink()

    print(f"[INFO] Installed to {outdir}")


#Detect input file formats

def detect_flag(path: Path) -> str:
#Detect input format for KMC
#FASTQ -> -fq
#FASTA -> -fm
#BAM   -> -fbam
    s = str(path).lower()

    if s.endswith((".fq", ".fastq", ".fq.gz", ".fastq.gz")):
        return "-fq"

    if s.endswith((".fa", ".fasta", ".fa.gz", ".fasta.gz")):
        return "-fm"

    if s.endswith(".bam"):
        return "-fbam"

    raise RuntimeError(f"Unsupported input format: {path}")


def validate_seq_formats(seqs):
    flags = {detect_flag(Path(x)) for x in seqs}
    if len(flags) != 1:
        raise RuntimeError("All input files in --seqs must have the same format.")
    return next(iter(flags))


#Helper functions

def calc_ratio(n, d):
    return n / d if d > 0 else 0.0


def norm_ratios(rh, re):
    s = rh + re
    return (rh / s, re / s) if s > 0 else (0.0, 0.0)


def classify_sex(kr_het, threshold):
    return "HET" if kr_het >= threshold else "HOM"


def format_sex(cls, sex_system):
    if sex_system == "XY":
        return "XY" if cls == "HET" else "XX"
    if sex_system == "ZW":
        return "ZW" if cls == "HET" else "ZZ"
    raise RuntimeError(f"Unsupported sex system: {sex_system}")


def parse_k_values(k):
    return [int(x) for x in str(k).split(",")]


def parse_threshold_values(threshold_arg, k_values):
    if threshold_arg is None:
        return None

    vals = [float(x) for x in str(threshold_arg).split(",")]
    if len(vals) == 1:
        return {str(k): vals[0] for k in k_values}

    if len(vals) != len(k_values):
        raise RuntimeError(
            "--threshold must be either one value or the same number of values as --k"
        )

    return {str(k): v for k, v in zip(k_values, vals)}


def get_total_kmers(db_prefix: Path) -> int:
    res = subprocess.run(
        ["kmc_tools", "info", str(db_prefix)],
        capture_output=True,
        text=True,
        check=True,
    )

    m = re.search(r"total k-mers\s*:\s*(\d+)", res.stdout)
    if not m:
        raise RuntimeError(f"Could not parse 'total k-mers' from kmc_tools info for {db_prefix}")

    return int(m.group(1))


def intersect_kmers(query_db, ref_db, out_db, verbose=False):
    run(
        ["kmc_tools", "simple", str(query_db), str(ref_db), "intersect", str(out_db)],
        verbose=verbose,
    )
    return get_total_kmers(out_db)


def write_filelist(paths, out_path: Path):
    with open(out_path, "w") as fh:
        for p in paths:
            fh.write(str(p) + "\n")
    return f"@{out_path}"


def remove_kmc_db(prefix):
    for ext in [".kmc_pre", ".kmc_suf"]:
        p = Path(str(prefix) + ext)
        if p.exists():
            p.unlink()

#Build ONYX database

def build(args):
    # Build ONYX k-mer database
    #
    # Overview:
    #   1. Split reference genome into individual contigs
    #   2. Identify homologous and heterologous sex chromosomes
    #   3. Construct autosome + sex chromosome references
    #   4. Count k-mers using KMC
    #   5. Extract unique sex-linked k-mers

    ensure_env(["seqkit", "kmc", "kmc_tools"])
    print(print_title())
    print("[INFO] Build started")

    ref = Path(args.ref)
    out = Path(args.out)
    k_values = parse_k_values(args.k)
    threshold_map = parse_threshold_values(args.threshold, k_values)

    contigs = out / "contigs"
    fasta_dir = out / "build_fasta"
    kmc_dir = out / "kmc"
    tmp_dir = out / "tmp"

    for d in [contigs, fasta_dir, kmc_dir, tmp_dir]:
        d.mkdir(parents=True, exist_ok=True)

    print("[INFO] Splitting reference")
    run(["seqkit", "split", "--by-id", str(ref), "-O", str(contigs)], verbose=args.verbose)

    hom_matches = list(contigs.glob(f"*.part_{args.sex_hom}.fa"))
    het_matches = list(contigs.glob(f"*.part_{args.sex_het}.fa"))

    if len(hom_matches) != 1:
        raise RuntimeError(f"{args.sex_hom} not uniquely found in split reference")
    if len(het_matches) != 1:
        raise RuntimeError(f"{args.sex_het} not uniquely found in split reference")

    hom = hom_matches[0]
    het = het_matches[0]

    autos = [f for f in contigs.glob("*.fa") if f not in (hom, het)]

    auto = fasta_dir / "autosome.fa"
    with open(auto, "wb") as o:
        for f in autos:
            o.write(f.read_bytes())

    auto_hom = fasta_dir / f"autosome_{args.sex_hom}.fa"
    auto_het = fasta_dir / f"autosome_{args.sex_het}.fa"

    auto_hom.write_bytes(auto.read_bytes() + hom.read_bytes())
    auto_het.write_bytes(auto.read_bytes() + het.read_bytes())

    meta = {}

    for k in k_values:
        print(f"[INFO] Building k={k}")

        kdir = kmc_dir / f"k{k}"
        ktmp = tmp_dir / f"k{k}"
        kdir.mkdir(parents=True, exist_ok=True)
        ktmp.mkdir(parents=True, exist_ok=True)

        hom_union = kdir / "HOM_union"
        het_union = kdir / "HET_union"
        hom_db = kdir / "hom_kmers"
        het_db = kdir / "het_kmers"

        run(
            ["kmc", f"-k{k}", f"-t{args.threads}", "-ci1", "-fm", str(auto_hom), str(hom_union), str(ktmp)],
            verbose=args.verbose,
        )
        run(
            ["kmc", f"-k{k}", f"-t{args.threads}", "-ci1", "-fm", str(auto_het), str(het_union), str(ktmp)],
            verbose=args.verbose,
        )

        run(
            ["kmc_tools", "simple", str(hom_union), str(het_union), "kmers_subtract", str(hom_db), "-ci1", "-cx1"],
            verbose=args.verbose,
        )
        run(
            ["kmc_tools", "simple", str(het_union), str(hom_union), "kmers_subtract", str(het_db), "-ci1", "-cx1"],
            verbose=args.verbose,
        )

        entry = {
            "hom_total_kmers": get_total_kmers(hom_db),
            "het_total_kmers": get_total_kmers(het_db),
            "hom_db": str(Path("kmc") / f"k{k}" / "hom_kmers"),
            "het_db": str(Path("kmc") / f"k{k}" / "het_kmers"),
        }

        if threshold_map is not None:
            entry["threshold"] = threshold_map[str(k)]

        meta[str(k)] = entry

        remove_kmc_db(hom_union)
        remove_kmc_db(het_union)

    build_info = {
        "reference": str(ref),
        "hom_name": args.sex_hom,
        "het_name": args.sex_het,
        "k_values": k_values,
        "preset": bool(args.preset),
        "kmer_metadata": meta,
    }

    with open(out / "build_info.json", "w") as fh:
        json.dump(build_info, fh, indent=4)

    print("[INFO] Cleaning temporary files")
    for d in [contigs, fasta_dir, tmp_dir]:
        if d.exists():
            shutil.rmtree(d)

    print("[INFO] Build finished")


#Run boostrap resampling

def bootstrap_resample_fastx(seqs, out_dir, fraction, seed, verbose=False):
    # Bootstrap resampling of sequencing reads from FASTA or FASTQ
    # Used to estimate confidence intervals of KR_het
    out_paths = []

    for idx, seq in enumerate(seqs):
        seq = Path(seq)
        name = seq.name
        out_path = out_dir / f"boot_{idx}_{name}"

        run(
            [
                "seqkit", "sample",
                "-p", str(fraction),
                "-s", str(seed),
                str(seq),
                "-o", str(out_path),
            ],
            verbose=verbose,
        )
        out_paths.append(out_path)

    return out_paths

def bootstrap_resample_bam(seqs, out_dir, fraction, seed, verbose=False):
    #Bootstrap resampling of BAM files using samtools.
    out_paths = []

    for idx, seq in enumerate(seqs):
        seq = Path(seq)
        name = seq.name
        out_path = out_dir / f"boot_{idx}_{name}"

        run(
            [
                "samtools",
                "view",
                "--subsample",
                str(fraction),
                "--subsample-seed",
                str(seed),
                "-b",
                str(seq),
                "-o",
                str(out_path),
            ],
            verbose=verbose,
        )

        out_paths.append(out_path)

    return out_paths

#Sex classification

def classify(args):
    # Classify biological sex using precomputed sex-linked k-mers
    #
    # Overview:
    #   1. Count k-mers from sequencing reads using KMC
    #   2. Intersect read k-mers with sex-specific k-mer databases
    #   3. Compute KR_hom and KR_het statistics
    #   4. Determine sex based on KR_het threshold
    ensure_env(["kmc", "kmc_tools", "seqkit","samtools"])
    print(print_title())

    db = Path(args.db)
    with open(db / "build_info.json") as fh:
        meta = json.load(fh)

    seqs = [Path(x) for x in args.seqs]
    flag = validate_seq_formats(seqs)

    seqs_str = ",".join([str(x) for x in seqs])
    sample = args.sample_id if args.sample_id else seqs[0].name

    print(f"[INFO] Inputs: {seqs_str}")

    temp = Path(tempfile.mkdtemp())

    try:
        with open(args.out, "w") as f:
            use_bootstrap_cols = args.bootstrap > 0

            if use_bootstrap_cols:
                header = (
                    "sample\tinputs\tk\tKR_hom\tKR_het\tclass\tsex\tsex_system\t"
                    "ci_low\tci_high\tconfidence\tbootstrap_n\tbootstrap_fraction\tbootstrap_seed\tthreshold\n"
                )
            else:
                header = "sample\tinputs\tk\tKR_hom\tKR_het\tclass\tsex\tsex_system\n"

            f.write(header)

            original_list = temp / "original_reads.txt"
            kmc_input = write_filelist(seqs, original_list) if len(seqs) > 1 else str(seqs[0])

            for k in meta["k_values"]:
                print(f"[INFO] Processing k={k}")

                q = temp / f"query_k{k}"
                run(
                    [
                        "kmc",
                        f"-k{k}",
                        f"-t{args.threads}",
                        f"-m{args.memory}",
                        "-ci1",
                        flag,
                        str(kmc_input),
                        str(q),
                        str(temp),
                    ],
                    verbose=args.verbose,
                )

                hom = intersect_kmers(q, db / meta["kmer_metadata"][str(k)]["hom_db"], temp / f"hom_{k}", verbose=args.verbose)
                het = intersect_kmers(q, db / meta["kmer_metadata"][str(k)]["het_db"], temp / f"het_{k}", verbose=args.verbose)

                hom_total = meta["kmer_metadata"][str(k)]["hom_total_kmers"]
                het_total = meta["kmer_metadata"][str(k)]["het_total_kmers"]

                rh = calc_ratio(hom, hom_total)
                re = calc_ratio(het, het_total)

                KR_hom, KR_het = norm_ratios(rh, re)

                if args.threshold is not None:
                    thr = args.threshold
                else:
                    thr = meta["kmer_metadata"][str(k)].get("threshold")

                if thr is not None:
                    cls = classify_sex(KR_het, thr)
                    sex = format_sex(cls, args.system)
                else:
                    cls = "NA"
                    sex = "NA"

                if use_bootstrap_cols:
                    if thr is None:
                        raise RuntimeError(
                            "Bootstrap requires a threshold. Provide --threshold or use a preset DB with threshold."
                        )

                    seed = args.bootstrap_seed if args.bootstrap_seed is not None else 42

                    print("[INFO] Running bootstrap")

                    KR_list = []
                    final = cls
                    match = 0

                    for i in trange(args.bootstrap, desc="Bootstrap"):
                        boot_dir = temp / f"boot_iter_{i}"
                        boot_dir.mkdir(parents=True, exist_ok=True)

                        if flag == "-fbam":
                            boot_paths = bootstrap_resample_bam(
                                seqs=seqs,
                                out_dir=boot_dir,
                                fraction=args.bootstrap_fraction,
                                seed=seed + i,
                                verbose=args.verbose,
                            )
                        else:
                            boot_paths = bootstrap_resample_fastx(
                                seqs=seqs,
                                out_dir=boot_dir,
                                fraction=args.bootstrap_fraction,
                                seed=seed + i,
                                verbose=args.verbose,
                            )


                        boot_list = boot_dir / "boot_reads.txt"
                        boot_input = write_filelist(boot_paths, boot_list) if len(boot_paths) > 1 else str(boot_paths[0])

                        qb = temp / f"qb_k{k}_{i}"
                        run(
                            [
                                "kmc",
                                f"-k{k}",
                                f"-t{args.threads}",
                                f"-m{args.memory}",
                                "-ci1",
                                flag,
                                str(boot_input),
                                str(qb),
                                str(temp),
                            ],
                            verbose=args.verbose,
                        )

                        hb = intersect_kmers(
                            qb,
                            db / meta["kmer_metadata"][str(k)]["hom_db"],
                            temp / f"hb_k{k}_{i}",
                            verbose=args.verbose,
                        )
                        eb = intersect_kmers(
                            qb,
                            db / meta["kmer_metadata"][str(k)]["het_db"],
                            temp / f"eb_k{k}_{i}",
                            verbose=args.verbose,
                        )

                        rhb = calc_ratio(hb, hom_total)
                        reb = calc_ratio(eb, het_total)

                        KRb = reb / (rhb + reb + 1e-12)
                        KR_list.append(KRb)

                        if classify_sex(KRb, thr) == final:
                            match += 1

                    KR_list.sort()
                    ci_low = KR_list[int(0.025 * len(KR_list))]
                    ci_high = KR_list[int(0.975 * len(KR_list))]
                    confidence = match / len(KR_list)

                    result = [
                        sample,
                        seqs_str,
                        k,
                        KR_hom,
                        KR_het,
                        cls,
                        sex,
                        args.system,
                        ci_low,
                        ci_high,
                        confidence,
                        args.bootstrap,
                        args.bootstrap_fraction,
                        seed,
                        thr,
                    ]

                    result_header = (
                        "sample", "inputs", "k", "KR_hom", "KR_het", "class", "sex", "sex_system",
                        "ci_low", "ci_high", "confidence", "bootstrap_n",
                        "bootstrap_fraction", "bootstrap_seed", "threshold"
                    )

                else:
                    result = [
                        sample,
                        seqs_str,
                        k,
                        KR_hom,
                        KR_het,
                        cls,
                        sex,
                        args.system,
                    ]

                    result_header = (
                        "sample", "inputs", "k", "KR_hom", "KR_het", "class", "sex", "sex_system"
                    )

                f.write("\t".join(map(str, result)) + "\n")

                print("[RESULT]")
                for h, v in zip(result_header, result):
                    print(f"  {h}: {v}")


    finally:
        shutil.rmtree(temp)


#Print command line inference of ONYX
#
# Subcommands:
#   build        build ONYX database
#   classify     infer biological sex
#   download-db  download preset databases

def main():
    if any(x in sys.argv for x in ["-h", "--help"]):
        print(print_title())

    p = argparse.ArgumentParser(
    prog="onyx",
    description="ONYX: An alignment-free biological sex inference from high-throughput sequencing data",
    formatter_class=CustomFormatter,
    epilog="""
Examples

Build database
  onyx build --ref genome.fa --sex_hom chrX --sex_het chrY

Classify sequencing reads
  onyx classify --seqs reads_R1.fq.gz reads_R2.fq.gz --db onyx_db --system XY --out result.tsv

Download preset database
  onyx download-db human --outdir db

List available preset databases
  onyx download-db --list
"""
)

    p.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"ONYX {VERSION}"
    )

    #Commands parser

    s = p.add_subparsers(
        dest="cmd",
        title="commands",
        metavar=""
    )

    #Build parser

    b = s.add_parser(
        "build",
        help="Build a k-mer database for sex inference",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    b.add_argument(
        "--ref",
        required=True,
        help="Reference genome FASTA containing autosomes and sex chromosomes"
    )

    b.add_argument(
        "--out",
        default="onyx_db",
        help="Output directory for the database"
    )

    b.add_argument(
        "--k",
        default="33",
        help="k-mer size(s). Multiple values can be specified with commas (e.g. 21,31,33)"
    )

    b.add_argument(
        "--sex_hom",
        required=True,
        help="Homologous sex chromosome name in the reference (e.g. chrX or chrZ)"
    )

    b.add_argument(
        "--sex_het",
        required=True,
        help="Heterologous sex chromosome name in the reference (e.g. chrY or chrW)"
    )

    b.add_argument(
        "--threads",
        type=int,
        default=8,
        help="Number of CPU threads used for k-mer counting"
    )

    b.add_argument(
        "--preset",
        action="store_true",
        help="Mark the database as a preset database in build_info.json"
    )

    b.add_argument(
        "--threshold",
        help="Sex classification threshold(s) embedded in the database. "
             "Provide a single value or comma-separated values matching --k"
    )

    b.add_argument(
        "--verbose",
        action="store_true",
        help="Print external tool commands and logs"
    )

    b.set_defaults(func=build)

    #Classification parser

    c = s.add_parser(
        "classify",
        help="Infer biological sex from sequencing reads",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    c.add_argument(
        "--seqs",
        nargs="+",
        required=True,
        help="Input sequencing files (FASTQ, FASTA, or BAM). "
             "Multiple files can be provided (e.g. paired-end reads)"
    )

    c.add_argument(
        "--sample-id",
        help="Sample identifier used in output. Default: first input filename"
    )

    c.add_argument(
        "--db",
        required=True,
        help="Path to ONYX database directory"
    )

    c.add_argument(
        "--out",
        required=True,
        help="Output TSV file for classification results"
    )

    c.add_argument(
        "--system",
        required=True,
        choices=["XY", "ZW"],
        help="Sex determination system of the organism"
    )

    c.add_argument(
        "--threshold",
        type=float,
        help="Override classification threshold instead of using database value"
    )

    c.add_argument(
        "--threads",
        type=int,
        default=8,
        help="Number of CPU threads used for k-mer counting"
    )

    c.add_argument(
        "--memory",
        type=int,
        default=12,
        help="Maximum RAM in GB used by KMC"
    )

    c.add_argument(
        "--bootstrap",
        type=int,
        default=0,
        help="Number of bootstrap iterations for confidence estimation"
    )

    c.add_argument(
        "--bootstrap-fraction",
        type=float,
        default=0.7,
        help="Fraction of reads sampled in each bootstrap iteration"
    )

    c.add_argument(
        "--bootstrap-seed",
        type=int,
        help="Base random seed used for bootstrap sampling"
    )

    c.add_argument(
        "--verbose",
        action="store_true",
        help="Print external tool commands and logs"
    )

    c.set_defaults(func=classify)

    #Download parser

    d = s.add_parser(
        "download-db",
        help="Download preset ONYX databases",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    d.add_argument(
        "name",
        nargs="?",
        help="Database name (e.g. human or chicken)"
    )

    d.add_argument(
        "--outdir",
        help="Directory where the database will be installed"
    )

    d.add_argument(
        "--list",
        action="store_true",
        help="List available preset databases"
    )

    d.set_defaults(func=download_db)

    args = p.parse_args()

    #Command validation

    if not hasattr(args, "func"):
        p.print_help()
        sys.exit(1)

    if args.cmd == "download-db":
        if not args.list and args.name is None:
            print("ERROR: database name required\n", file=sys.stderr)
            print("Try:\n  onyx download-db --list\n", file=sys.stderr)
            sys.exit(1)

    args.func(args)

if __name__ == "__main__":
    main()


    
