#!/usr/bin/env python3
"""
pdb2fasta.py — convert PDB files to FASTA (one record per chain). No dependencies.

By default the sequence is read from ATOM/HETATM records (what is actually
modelled). Use --seqres to read the deposited SEQRES sequence instead.

Examples:
    python pdb2fasta.py design.pdb                 # prints FASTA to stdout
    python pdb2fasta.py *.pdb -o all.fasta         # many files -> one FASTA
    python pdb2fasta.py model.pdb --chains A B     # only chains A and B
    python pdb2fasta.py model.pdb --merge          # all chains joined with '/'
    python pdb2fasta.py 1abc.pdb --seqres          # use SEQRES records
    python pdb2fasta.py model.pdb --model all      # every MODEL in NMR/ensemble files
"""
import argparse
import gzip
import sys
from pathlib import Path

THREE_TO_ONE = {
    # standard
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "SEC": "U", "PYL": "O",
    # ambiguous
    "ASX": "B", "GLX": "Z", "UNK": "X",
    # common modified residues -> parent
    "MSE": "M", "SEP": "S", "TPO": "T", "PTR": "Y", "CSO": "C",
    "CSD": "C", "CME": "C", "OCS": "C", "HYP": "P", "MLY": "K",
    "M3L": "K", "KCX": "K", "LLP": "K", "PCA": "E", "CGU": "E",
    "NLE": "L", "FME": "M", "SAC": "S", "DAL": "A",
    # protonation-state / force-field names (Amber, CHARMM, Rosetta)
    "HID": "H", "HIE": "H", "HIP": "H", "HSD": "H", "HSE": "H",
    "HSP": "H", "CYX": "C", "CYM": "C", "ASH": "D", "GLH": "E",
    "LYN": "K",
}
NUCLEOTIDES = {
    "A": "A", "C": "C", "G": "G", "U": "U", "T": "T", "I": "I",
    "DA": "A", "DC": "C", "DG": "G", "DT": "T", "DU": "U", "DI": "I",
}


def to_one(resname, unknown="X"):
    resname = resname.strip().upper()
    return THREE_TO_ONE.get(resname) or NUCLEOTIDES.get(resname) or unknown


def open_any(path):
    path = str(path)
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def parse_atoms(lines, model="first", keep_hetatm=True, fill_gaps=False):
    """Return {model_id: {chain: sequence}} from ATOM/HETATM records."""
    models = {}
    current_model = "1"
    seen_model_record = False
    chains = None
    last_res = {}   # chain -> (resseq, icode)
    last_num = {}   # chain -> resseq int, for gap filling

    def start_model(mid):
        nonlocal chains
        chains = models.setdefault(mid, {})
        last_res.clear()
        last_num.clear()

    start_model(current_model)
    for line in lines:
        rec = line[:6].strip()
        if rec == "MODEL":
            if seen_model_record and model == "first":
                break
            seen_model_record = True
            current_model = line[10:14].strip() or str(len(models) + 1)
            if not models.get("1"):
                models.pop("1", None)
            start_model(current_model)
            continue
        if rec == "ENDMDL" and model == "first":
            break
        if rec not in ("ATOM", "HETATM"):
            continue
        resname = line[17:20].strip()
        if rec == "HETATM":
            if not keep_hetatm or resname not in THREE_TO_ONE:
                continue  # skip waters, ligands, ions, unknown hets
        altloc = line[16:17]
        if altloc not in (" ", "A", "1"):
            continue
        chain = line[21:22].strip() or "_"
        try:
            resseq = int(line[22:26])
        except ValueError:
            continue
        icode = line[26:27]
        key = (resseq, icode)
        if last_res.get(chain) == key:
            continue  # same residue, next atom
        seq = chains.setdefault(chain, [])
        if fill_gaps and chain in last_num and icode == " ":
            gap = resseq - last_num[chain] - 1
            if 0 < gap < 1000:
                seq.append("X" * gap)
        seq.append(to_one(resname))
        last_res[chain] = key
        last_num[chain] = resseq

    return {m: {c: "".join(s) for c, s in ch.items()} for m, ch in models.items() if ch}


def parse_seqres(lines):
    chains = {}
    for line in lines:
        if line.startswith("SEQRES"):
            chain = line[11:12].strip() or "_"
            chains.setdefault(chain, []).extend(to_one(r) for r in line[19:].split())
    return {"1": {c: "".join(s) for c, s in chains.items()}}


def wrap(seq, width):
    if width <= 0:
        return seq
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def main():
    ap = argparse.ArgumentParser(description="Convert PDB file(s) to FASTA.")
    ap.add_argument("pdb", nargs="+", help="PDB file(s); .gz ok; '-' for stdin")
    ap.add_argument("-o", "--output", help="output FASTA (default: stdout)")
    ap.add_argument("--seqres", action="store_true", help="use SEQRES instead of coordinates")
    ap.add_argument("--chains", nargs="+", help="only these chain IDs")
    ap.add_argument("--merge", action="store_true", help="one record per file, chains joined by --sep")
    ap.add_argument("--sep", default="/", help="chain separator for --merge (default '/')")
    ap.add_argument("--model", choices=["first", "all"], default="first", help="which MODEL(s) to read")
    ap.add_argument("--no-hetatm", action="store_true", help="ignore HETATM (drops MSE etc.)")
    ap.add_argument("--fill-gaps", action="store_true", help="insert X for missing residue numbers")
    ap.add_argument("--min-length", type=int, default=1, help="skip chains shorter than this")
    ap.add_argument("--width", type=int, default=0, help="line width (0 = single line)")
    args = ap.parse_args()

    out = open(args.output, "w") if args.output else sys.stdout
    n_records = 0
    for path in args.pdb:
        try:
            with open_any(path) as fh:
                lines = fh.readlines()
        except OSError as e:
            print(f"warning: {e}", file=sys.stderr)
            continue
        name = "stdin" if path == "-" else Path(path).name
        for ext in (".gz", ".pdb", ".ent"):
            name = name[: -len(ext)] if name.endswith(ext) else name

        if args.seqres:
            models = parse_seqres(lines)
            if not models["1"]:
                print(f"warning: no SEQRES in {path}", file=sys.stderr)
                continue
        else:
            models = parse_atoms(lines, args.model, not args.no_hetatm, args.fill_gaps)

        for mid, chains in models.items():
            if args.chains:
                chains = {c: s for c, s in chains.items() if c in args.chains}
            chains = {c: s for c, s in chains.items() if len(s) >= args.min_length}
            if not chains:
                continue
            tag = f"_model{mid}" if len(models) > 1 else ""
            if args.merge:
                ids = "".join(chains)
                out.write(f">{name}{tag}_{ids}\n{wrap(args.sep.join(chains.values()), args.width)}\n")
                n_records += 1
            else:
                for c, s in chains.items():
                    out.write(f">{name}{tag}_{c} length={len(s)}\n{wrap(s, args.width)}\n")
                    n_records += 1

    if args.output:
        out.close()
        print(f"wrote {n_records} sequence(s) to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
