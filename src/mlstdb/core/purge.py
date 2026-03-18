import os
import shutil
from pathlib import Path

import yaml

from mlstdb.utils import error, info, success


def parse_profile_file(profile_path: str):
    """Parse a scheme profile TSV file.

    Returns (header, locus_columns, rows) where:
      - header: list of all column names
      - locus_columns: list of locus column names (between ST and clonal_complex)
      - rows: list of dicts mapping column name → value for each ST row
    """
    with open(profile_path, "r") as f:
        header_line = f.readline().rstrip("\n")
        header = header_line.split("\t")

        # Determine locus columns: everything after ST, excluding clonal_complex
        st_idx = header.index("ST")
        if "clonal_complex" in header:
            cc_idx = header.index("clonal_complex")
            locus_columns = header[st_idx + 1 : cc_idx]
        else:
            locus_columns = header[st_idx + 1 :]

        rows = []
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            values = line.split("\t")
            row = dict(zip(header, values))
            rows.append(row)

    return header, locus_columns, rows


def write_profile_file(profile_path: str, header: list, rows: list):
    """Write rows back to a profile TSV file."""
    with open(profile_path, "w") as f:
        f.write("\t".join(header) + "\n")
        for row in rows:
            f.write("\t".join(row.get(col, "") for col in header) + "\n")


def remove_alleles_from_fasta(fasta_path: str, alleles_to_remove: set):
    """Remove specified allele entries from a FASTA file.

    alleles_to_remove: set of full header names (e.g. {"aroC_1", "aroC_3"})
    """
    entries = []
    current_header = None
    current_seq_lines = []

    with open(fasta_path, "r") as f:
        for line in f:
            if line.startswith(">"):
                if current_header is not None:
                    entries.append((current_header, current_seq_lines))
                current_header = line[1:].strip()
                current_seq_lines = []
            else:
                current_seq_lines.append(line)
        if current_header is not None:
            entries.append((current_header, current_seq_lines))

    with open(fasta_path, "w") as f:
        for header, seq_lines in entries:
            if header not in alleles_to_remove:
                f.write(f">{header}\n")
                for seq_line in seq_lines:
                    f.write(seq_line)


def parse_purge_config(config_path: str):
    """Parse and validate a YAML purge configuration file.

    Returns (entries, global_opts) where:
      - entries: list of dicts with keys 'scheme', 'st' (list), 'alleles' (list)
      - global_opts: dict with optional 'force', 'verbose', 'directory', 'blast_directory'
    """
    with open(config_path, "r") as f:
        data = yaml.safe_load(f)

    if not isinstance(data, dict) or "purge" not in data:
        raise ValueError("Config file must contain a top-level 'purge' key with a list of entries.")

    raw_entries = data["purge"]
    if not isinstance(raw_entries, list) or len(raw_entries) == 0:
        raise ValueError("'purge' must be a non-empty list of scheme entries.")

    entries = []
    for i, entry in enumerate(raw_entries):
        if not isinstance(entry, dict) or "scheme" not in entry:
            raise ValueError(f"Entry {i + 1} must be a dict with at least a 'scheme' key.")

        parsed = {"scheme": str(entry["scheme"]), "st": [], "alleles": []}

        if "st" in entry:
            st_val = entry["st"]
            if not isinstance(st_val, list):
                st_val = [st_val]
            parsed["st"] = [str(s) for s in st_val]

        if "alleles" in entry:
            alleles_val = entry["alleles"]
            if not isinstance(alleles_val, list):
                alleles_val = [alleles_val]
            for a in alleles_val:
                a_str = str(a)
                if ":" not in a_str:
                    raise ValueError(
                        f"Entry {i + 1}: invalid allele format '{a_str}'. "
                        "Expected locus:number (e.g. aroC:3)."
                    )
            parsed["alleles"] = [str(a) for a in alleles_val]

        entries.append(parsed)

    global_opts = {
        "force": bool(data.get("force", False)),
        "verbose": bool(data.get("verbose", False)),
        "directory": data.get("directory"),
        "blast_directory": data.get("blast_directory"),
    }

    return entries, global_opts


def _format_st_list(st_numbers: list) -> str:
    """Format a list of ST numbers for display: show first 3, then 'and N other STs'."""
    sorted_sts = sorted(st_numbers, key=int)
    if len(sorted_sts) <= 3:
        return ", ".join(f"ST {s}" for s in sorted_sts)
    shown = ", ".join(f"ST {s}" for s in sorted_sts[:3])
    remaining = len(sorted_sts) - 3
    return f"{shown} and {remaining} other STs"


def purge_scheme(scheme_dir: str, verbose: bool = False):
    """Purge an entire scheme directory (Mode 1)."""
    shutil.rmtree(scheme_dir)
    success(f"Scheme directory removed: {scheme_dir}")


def purge_st(scheme_dir: str, scheme_name: str, st_number: str,
             force: bool = False, verbose: bool = False):
    """Purge a single ST from a scheme (Mode 2).

    Returns True if changes were made, False otherwise.
    """
    profile_path = os.path.join(scheme_dir, f"{scheme_name}.txt")
    if not os.path.exists(profile_path):
        error(f"Profile file not found: {profile_path}")
        return False

    header, locus_columns, rows = parse_profile_file(profile_path)

    # Find the ST row to remove
    target_row = None
    remaining_rows = []
    for row in rows:
        if row["ST"] == st_number:
            target_row = row
        else:
            remaining_rows.append(row)

    if target_row is None:
        error(f"ST {st_number} not found in {scheme_name}.")
        return False

    if verbose:
        info(f"Removing ST {st_number} from {scheme_name}.")

    # Build set of allele numbers still in use per locus (from remaining rows)
    alleles_in_use = {}
    for locus in locus_columns:
        alleles_in_use[locus] = set()
        for row in remaining_rows:
            alleles_in_use[locus].add(row[locus])

    # Check each allele in the target ST row
    for locus in locus_columns:
        allele_num = target_row[locus]
        allele_name = f"{locus}_{allele_num}"
        fasta_path = os.path.join(scheme_dir, f"{locus}.tfa")

        if allele_num in alleles_in_use[locus]:
            # Allele is still used by other STs
            using_sts = [r["ST"] for r in remaining_rows if r[locus] == allele_num]
            msg = f"Allele {allele_name} is still used by {_format_st_list(using_sts)} \u2014 skipping."
            if force:
                info(f"Allele {allele_name} is still used by {_format_st_list(using_sts)} \u2014 force deleting.")
                if os.path.exists(fasta_path):
                    remove_alleles_from_fasta(fasta_path, {allele_name})
                    if verbose:
                        info(f"Removed allele {allele_name} from {fasta_path}")
            else:
                info(msg)
        else:
            # Allele is orphaned — remove it
            if os.path.exists(fasta_path):
                remove_alleles_from_fasta(fasta_path, {allele_name})
                if verbose:
                    info(f"Removed orphaned allele {allele_name} from {fasta_path}")

    # Write updated profile file
    write_profile_file(profile_path, header, remaining_rows)
    success(f"ST {st_number} removed from {scheme_name}.")
    return True


def purge_allele(scheme_dir: str, scheme_name: str, locus: str, allele_num: str,
                 force: bool = False, verbose: bool = False):
    """Purge an allele and all STs that reference it (Mode 3).

    Returns True if changes were made, False otherwise.
    """
    profile_path = os.path.join(scheme_dir, f"{scheme_name}.txt")
    fasta_path = os.path.join(scheme_dir, f"{locus}.tfa")

    if not os.path.exists(profile_path):
        error(f"Profile file not found: {profile_path}")
        return False
    if not os.path.exists(fasta_path):
        error(f"Allele file not found: {fasta_path}")
        return False

    header, locus_columns, rows = parse_profile_file(profile_path)

    if locus not in locus_columns:
        error(f"Locus '{locus}' not found in scheme {scheme_name}.")
        return False

    # Find affected STs
    affected_sts = [row["ST"] for row in rows if row[locus] == allele_num]
    remaining_rows = [row for row in rows if row[locus] != allele_num]

    allele_name = f"{locus}_{allele_num}"

    if not affected_sts:
        info(f"No STs reference allele {allele_name}.")
    else:
        info(f"Allele {allele_name} is used by {len(affected_sts)} STs: {_format_st_list(affected_sts)}.")

    if not force:
        import click
        if not click.confirm(
            f"Remove allele {allele_name} and {len(affected_sts)} affected ST(s)?",
            default=False,
        ):
            info("Purge cancelled.")
            return False

    # Remove allele from FASTA
    remove_alleles_from_fasta(fasta_path, {allele_name})
    if verbose:
        info(f"Removed allele {allele_name} from {fasta_path}")

    # Remove affected ST rows
    write_profile_file(profile_path, header, remaining_rows)
    if affected_sts:
        success(f"Removed {len(affected_sts)} ST(s) referencing allele {allele_name}.")
    success(f"Allele {allele_name} purged from {scheme_name}.")
    return True


def purge_st_and_allele(scheme_dir: str, scheme_name: str, st_number: str,
                        locus: str, allele_num: str,
                        force: bool = False, verbose: bool = False):
    """Purge a specific ST and check a specific allele (Mode 4).

    Returns True if changes were made, False otherwise.
    """
    profile_path = os.path.join(scheme_dir, f"{scheme_name}.txt")
    fasta_path = os.path.join(scheme_dir, f"{locus}.tfa")

    if not os.path.exists(profile_path):
        error(f"Profile file not found: {profile_path}")
        return False

    header, locus_columns, rows = parse_profile_file(profile_path)

    if locus not in locus_columns:
        error(f"Locus '{locus}' not found in scheme {scheme_name}.")
        return False

    # Find and remove the ST row
    target_row = None
    remaining_rows = []
    for row in rows:
        if row["ST"] == st_number:
            target_row = row
        else:
            remaining_rows.append(row)

    if target_row is None:
        error(f"ST {st_number} not found in {scheme_name}.")
        return False

    if verbose:
        info(f"Removing ST {st_number} from {scheme_name}.")

    # Write updated profile (without target ST)
    write_profile_file(profile_path, header, remaining_rows)
    success(f"ST {st_number} removed from {scheme_name}.")

    # Check if the specified allele is still used
    allele_name = f"{locus}_{allele_num}"
    using_sts = [r["ST"] for r in remaining_rows if r[locus] == allele_num]

    if using_sts:
        msg = f"Allele {allele_name} is still used by {_format_st_list(using_sts)} \u2014 skipping."
        if force:
            info(f"Allele {allele_name} is still used by {_format_st_list(using_sts)} \u2014 force deleting.")
            if os.path.exists(fasta_path):
                remove_alleles_from_fasta(fasta_path, {allele_name})
                if verbose:
                    info(f"Removed allele {allele_name} from {fasta_path}")
        else:
            info(msg)
    else:
        # Orphaned — remove
        if os.path.exists(fasta_path):
            remove_alleles_from_fasta(fasta_path, {allele_name})
            if verbose:
                info(f"Removed orphaned allele {allele_name} from {fasta_path}")
            success(f"Orphaned allele {allele_name} removed.")

    return True
