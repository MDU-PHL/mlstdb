# Getting Started

This guide walks you through the complete setup, from installation to running `mlst` with your freshly updated database.

---

## Step 1: Install mlstdb

```sh
conda create -n mlst -c bioconda mlst
conda activate mlst
pip install mlstdb
```

Verify the installation:

```sh
mlstdb --version
```

---

## Step 2: Register with the databases

Before downloading any schemes, you need to register your OAuth credentials with PubMLST and/or Pasteur. This is a **one-time setup**. Your credentials are saved locally and reused for future updates.

### Connect to PubMLST

```sh
mlstdb connect --db pubmlst
```

### Connect to Pasteur

```sh
mlstdb connect --db pasteur
```

Each `connect` command will:

1. Ask for your **Client ID** (24 characters) and **Client Secret** (42 characters)
2. Open an authorisation URL, visit it in your browser
3. Ask you to paste the **verification code** from the website
4. Save all tokens securely to `~/.config/mlstdb/`

!!! tip "Where do I get my Client ID and Client Secret?"
    See the [Connect guide](usage/connect.md#obtaining-client-credentials) for step-by-step instructions on registering with PubMLST and Pasteur.

!!! info
    If you've already connected before, `mlstdb connect` will test your existing credentials and skip re-registration if they're still valid.

---

## Step 3: Download MLST schemes

```sh
mlstdb update
```

This will:

- Read the built-in curated list of ~170 MLST schemes from both PubMLST and Pasteur
- Download allele sequences (`.tfa` files) and ST profiles (`.txt` files) for each scheme
- Save everything to a `pubmlst/` directory
- Build a BLAST database in `blast/`

!!! note "First run may take a while"
    Downloading hundreds of schemes involves many API calls. You can speed things up with `--threads 4` or download specific schemes by providing a custom input file. See the [Update guide](usage/update.md) for details.

If the download is interrupted, use `--resume` to pick up where you left off:

```sh
mlstdb update --resume
```

---

## Step 4: Verify with mlst

Once the update is complete, test your new database:

```sh
mlst --blastdb blast/mlst.fa --datadir pubmlst your_assembly.fasta
```

Replace `your_assembly.fasta` with the path to any bacterial genome assembly.

---

## What was created?

After a successful update, your directory should look like this:

```
pubmlst/
├── klebsiella/
│   ├── klebsiella.txt          # ST profiles
│   ├── klebsiella_info.json    # Scheme metadata
│   ├── gapA.tfa                # Allele sequences
│   ├── infB.tfa
│   └── ...
├── listeria/
│   └── ...
└── ...

blast/
├── mlst.fa                     # Combined allele sequences
├── mlst.fa.ndb                 # BLAST index files
├── mlst.fa.nhr
└── ...
```

Each scheme gets its own subdirectory containing:

- **Profile file** (`<scheme>.txt`) : maps ST numbers to allele combinations
- **Allele files** (`<locus>.tfa`) : FASTA sequences for each locus
- **Metadata** (`<scheme>_info.json`) : source database, download date, locus count

---

## Keeping your database up to date

Schemes change over time as new STs and alleles are added. To update your database, simply run:

```sh
mlstdb update
```

This will re-download all schemes and rebuild the BLAST database. Use `--resume` if you want to skip schemes that have already been downloaded.

---

## Next Steps

- [Connect — Registration details](usage/connect.md) — How to obtain OAuth credentials
- [Update — All options](usage/update.md) — Custom inputs, parallel downloads, resume
- [Fetch — Advanced](usage/fetch.md) — Explore all available schemes with custom filters
- [Disclaimer](disclaimer.md) — Important safety notes
