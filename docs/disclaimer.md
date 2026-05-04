# Disclaimer

!!! warning "Please read before using mlstdb"

Please take note of the following before using `mlstdb` to update your MLST databases.

## Back up your databases

Always back up your original MLST databases before running any updates. This protects against accidental overwrites or deletions.

```sh
cp -r /path/to/existing/pubmlst /path/to/backup/pubmlst_backup
cp -r /path/to/existing/blast /path/to/backup/blast_backup
```

## Curate your scheme list

Not all MLST schemes available on PubMLST and Pasteur are validated for use with the `mlst` tool. Before updating:

- Review the scheme list if using a custom input file
- Don't blindly download and apply all available schemes
- Overwriting core MLST data with unverified schemes may cause downstream issues

## Bacterial species only

The `mlst` tool is designed for typing **bacterial species only**. If you use the advanced `fetch` command to explore schemes, make sure to filter out non-bacterial schemes before updating.

## Authentication requirements

Some schemes require registration with the respective database. If you encounter authentication errors during `update`, check that you have:

1. Registered with both PubMLST and Pasteur (via `mlstdb connect`)
2. Enrolled in the specific databases you need within each platform

## Network considerations

`mlstdb` makes API calls to external servers (PubMLST and Pasteur). Be aware of:

- Network interruptions — use `--resume` to recover
- Rate limiting — keep `--threads` at 4 or below
- Firewall restrictions — ensure outbound HTTPS access to `rest.pubmlst.org` and `bigsdb.pasteur.fr`
