# Purge

Sometimes a sequence type or allele ends up in your local database that really shouldn't be there. Maybe it's been flagged as contaminated, misassigned, or simply erroneous. Rather than wiping an entire scheme and re-downloading everything, `mlstdb purge` lets you surgically remove just the offending entry and rebuild the BLAST database in one go.

---

## Basic Usage

```sh
# Remove an entire scheme
mlstdb purge --scheme salmonella

# Remove a single ST from a scheme
mlstdb purge --scheme salmonella --st 3

# Remove a specific allele (and all STs that reference it)
mlstdb purge --scheme salmonella --allele aroC:1

# Batch purge across multiple schemes from a YAML config
mlstdb purge --config purge_config.yaml
```

---

## Options

| Option | Short | Description | Default |
|--------|-------|-------------|---------|
| `--scheme` | `-s` | Scheme name to purge (e.g. `salmonella`) | — |
| `--st` | | ST number to remove | — |
| `--allele` | `-a` | Allele to remove, format `locus:number` (e.g. `aroC:3`) | — |
| `--config` | `-c` | Path to a YAML config file for batch purging | — |
| `--force` | `-f` | Skip confirmation prompts and force-delete shared alleles | Off |
| `--verbose` | `-v` | Show detailed output | Off |
| `--directory` | `-d` | Directory containing downloaded schemes | `pubmlst` |
| `--blast-directory` | `-b` | Directory for the BLAST database | `blast` |
| `-h`, `--help` | | Show help | |

!!! note
    Either `--scheme` or `--config` is required. You cannot mix `--config` with `--scheme`, `--st`, or `--allele`. Put those in the config file instead.

---

## Purge modes

### Purge an entire scheme

```sh
mlstdb purge --scheme salmonella
```

Deletes the entire `pubmlst/salmonella/` directory (allele files, profile file, metadata, everything), and then rebuilds the BLAST database from the remaining schemes. You'll be asked to confirm before anything is deleted.

Use `--force` to skip the confirmation:

```sh
mlstdb purge --scheme salmonella --force
```

---

### Purge a specific ST

```sh
mlstdb purge --scheme salmonella --st 3
```

Removes the row for ST 3 from `salmonella.txt`. For each allele in that row, `mlstdb purge` then checks whether the allele is referenced by any other ST:

- **Not used elsewhere** → the allele entry is removed from the corresponding `.tfa` file (it's now an orphan. There's no point keeping it around)
- **Still used by other STs** → a warning is printed and the allele is left untouched:

```
Allele aroC_1 is still used by ST 1, ST 2, ST 8 and 5 other STs — skipping.
```

If you want to force-delete shared alleles regardless, add `--force`:

```sh
mlstdb purge --scheme salmonella --st 3 --force
```

The BLAST database is rebuilt automatically at the end.

---

### Purge a specific allele

```sh
mlstdb purge --scheme salmonella --allele aroC:1
```

Before removing anything, `mlstdb purge` checks which STs reference allele `aroC_1` and reports them:

```
Allele aroC_1 is used by 6 STs: ST 1, ST 2, ST 8 and 3 others.
Remove allele aroC_1 and 6 affected ST(s)? [y/N]:
```

If you confirm (or use `--force`), it will:

1. Remove the `aroC_1` entry from `aroC.tfa`
2. Remove all affected ST rows from `salmonella.txt`
3. Rebuild the BLAST database

This is the most thorough option. If a sequence itself is bad, it makes sense to clean out every ST built on it.

---

### Purge an ST and check a specific allele

```sh
mlstdb purge --scheme salmonella --st 3 --allele aroC:1
```

A targeted combination: removes ST 3, then specifically checks whether allele `aroC_1` has become an orphan as a result. If it has, it's removed too. If other STs still reference it, you'll see a warning (and `--force` overrides this as usual).

---

## Batch purge with a config file

If you need to clean up multiple schemes in one go, use a YAML config file, the BLAST database is only rebuilt once at the end, which is much faster than running `purge` separately for each scheme.

### Config format

```yaml
# purge_config.yaml

purge:
  # Remove two alleles from salmonella
  - scheme: salmonella
    alleles:
      - aroC:1
      - dnaN:5

  # Remove two STs from klebsiella
  - scheme: klebsiella
    st:
      - 3
      - 15

  # Remove an ST and an allele from listeria_2
  - scheme: listeria_2
    st:
      - 42
    alleles:
      - abcZ:12

  # Remove the entire bordetella scheme
  - scheme: bordetella

# Optional global settings (can be overridden by CLI flags)
force: false
verbose: true
directory: pubmlst
blast_directory: blast
```

Then run:

```sh
mlstdb purge --config purge_config.yaml
```

CLI flags override the global settings in the config file. For example, to force all operations without being prompted:

```sh
mlstdb purge --config purge_config.yaml --force
```

### Config keys

| Key | Type | Description |
|-----|------|-------------|
| `purge` | list | Required. List of scheme entries to process. |
| `scheme` | string | Required per entry. Must match a directory under `--directory`. |
| `st` | int or list of ints | ST number(s) to remove. |
| `alleles` | string or list of strings | Allele(s) to remove, in `locus:number` format. |
| `force` | bool | Skip confirmation prompts (default: `false`). |
| `verbose` | bool | Show detailed output (default: `false`). |
| `directory` | string | Schemes directory (default: `pubmlst`). |
| `blast_directory` | string | BLAST output directory (default: `blast`). |

---

## What gets rebuilt

After every purge operation, `mlstdb purge` automatically rebuilds the BLAST database from the remaining schemes in the data directory. This is the same process as `mlstdb update`. It concatenates all `.tfa` files across all scheme directories and runs `makeblastdb`.

!!! warning "Purging is permanent"
    Purged entries are deleted from disk and cannot be recovered. If you're unsure, take a backup of your `pubmlst/` directory before running `purge`.

    ```sh
    cp -r pubmlst/ pubmlst_backup/
    ```

    To restore a purged scheme from scratch, just re-run `mlstdb update`. It will download it again.

---

## Tips

- **Orphaned alleles** are sequences no longer linked to any ST. `mlstdb purge` automatically removes them when you purge an ST (unless they're still in use elsewhere). This keeps your FASTA files tidy.

- **Using `--force` on allele removal** will delete the allele even if other STs reference it, which may leave those STs with a missing allele entry. Only use `--force` if you're confident all the affected STs should go too or follow up with another purge to clean them up.

- **Custom directories**: If your schemes aren't in the default `pubmlst/` directory, use `-d` and `-b` to point `purge` to the right places:

    ```sh
    mlstdb purge --scheme salmonella --st 3 -d my_schemes -b my_blast
    ```
