# Update

The `update` command downloads MLST schemes and creates a BLAST database. This is the main command for keeping your `mlst` databases current.

---

## Basic Usage

```sh
mlstdb update
```

By default, this uses the built-in curated list of ~300 MLST schemes from both PubMLST and Pasteur.

---

## Options

| Option | Description | Default |
|--------|-------------|---------|
| `--input`, `-i` | Path to custom scheme list file | Built-in `mlst_schemes_all.tab` |
| `--directory`, `-d` | Output directory for scheme data | `pubmlst` |
| `--blast-directory`, `-b` | Output directory for BLAST database | `blast` |
| `--no-auth` | Skip OAuth; use unauthenticated access | Off |
| `--resume`, `-r` | Skip already-downloaded schemes | Off |
| `--threads`, `-t` | Parallel download threads (max recommended: 4) | `1` |
| `--verbose`, `-v` | Show detailed debug output | Off |
| `-h`, `--help` | Show help message | |

---

## Examples

### Download all curated schemes

```sh
mlstdb update
```

Uses the built-in curated scheme list. Requires prior setup with `mlstdb connect` for both PubMLST and Pasteur.

### Resume an interrupted download

```sh
mlstdb update --resume
```

Skips schemes whose profile file already exists in the output directory.

### Speed up with parallel downloads

```sh
mlstdb update --threads 4
```

Downloads up to 4 schemes simultaneously. Keep this at 4 or below to avoid overwhelming the API servers.

### Custom output directories

```sh
mlstdb update --directory my_schemes --blast-directory my_blast
```

### Use a custom scheme list

```sh
mlstdb update --input my_filtered_schemes.tab
```

The input file must be tab-delimited with these columns:

```
database	species	scheme_description	scheme	URI
pubmlst	Klebsiella pneumoniae	MLST	klebsiella	https://rest.pubmlst.org/db/pubmlst_klebsiella_seqdef/schemes/1
pasteur	Bordetella	MLST	bordetella_3	https://bigsdb.pasteur.fr/api/db/pubmlst_bordetella_seqdef/schemes/3
```

!!! tip "Generating a custom scheme list"
    Use `mlstdb fetch` with filters to generate a scheme list, then edit it to keep only the schemes you want. See the [Fetch guide](fetch.md) for details.

### Unauthenticated access

```sh
mlstdb update --no-auth
```

!!! warning "Limited availability"
    Unauthenticated access only works for data created before 2024-12-31. Newer schemes require authentication. Use `mlstdb connect` for full access.

---

## What gets downloaded

For each scheme, `mlstdb update` creates a directory containing:

| File | Description |
|------|-------------|
| `<scheme>.txt` | ST profile definitions (tab-delimited) |
| `<locus>.tfa` | Allele sequences for each locus (FASTA) |
| `<scheme>_info.json` | Metadata: source, download date, locus count |
| `database_version.txt` | Database version number |

### Profile file format

The profile file maps sequence type (ST) numbers to allele combinations:

```
ST	gapA	infB	mdh	pgi	phoE	rpoB	tonB
1	4	4	1	1	7	4	10
2	3	4	1	1	9	4	17
```

### Allele file format

Each `.tfa` file contains FASTA-formatted allele sequences:

```
>gapA_1
AACCTGAAGTGGGACGAAGTTGGTGTTGACGTTGTTGCTGAAG...
>gapA_2
AACCTGAAGTGGGACGAAGTTGGTGTTGACGTTGTTGCTGAAG...
```

---

## BLAST Database

After downloading all schemes, `update` concatenates the allele sequences and runs `makeblastdb` to produce the BLAST database:

```
blast/
├── mlst.fa           # Combined allele sequences
├── mlst.fa.ndb       # BLAST index files
├── mlst.fa.nhr
├── mlst.fa.nin
└── ...
```

This is the database you pass to `mlst` via `--blastdb blast/mlst.fa`.

---

## Handling authentication errors

If some schemes fail to download due to authentication issues, `update` will:

1. Report which schemes were skipped
2. Ask if you want to continue and build the BLAST database with the schemes that succeeded

Common causes:

- **Not registered** with a specific database: register via the database's web interface
- **Expired tokens**: re-run `mlstdb connect` to register new credentials

---

## Troubleshooting

### "Input file not found"

Run `mlstdb connect` first, and then `mlstdb update`. Or provide a valid scheme file with `--input`.

### "No credentials available"

Run `mlstdb connect --db pubmlst` and/or `mlstdb connect --db pasteur` to set up authentication.

### "No schemes were successfully downloaded"

- Check your internet connection
- Verify your credentials are valid: `mlstdb connect --db pubmlst`
- Try with `--verbose` to see detailed error messages

### Download is slow

- Use `--threads 4` for parallel downloads
- Use `--resume` if you need to restart, it won't re-download completed schemes
