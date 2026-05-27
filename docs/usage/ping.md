# Ping

The `ping` command sends a single GET request to an API endpoint and displays the HTTP status code and response body. It is useful for verifying that your credentials are working or for exploring raw API responses.

## Basic Usage

```sh
mlstdb ping <URL> [--db pubmlst|pasteur] [--no-auth] [--verbose]
```

## Options

| Option | Description |
|--------|-------------|
| `URL` | The API endpoint to probe (required). |
| `--db`, `-d` | Database whose stored credentials should be used: `pubmlst` or `pasteur`. If omitted and authentication is needed, you will be prompted to choose. |
| `--no-auth` | Skip all authentication and send an unauthenticated request. |
| `--verbose`, `-v` | Print the auth mode, full URL, and status code before the response body. |
| `-h`, `--help` | Show help message. |

## Authentication Priority

Unless `--no-auth` is set, `ping` attempts authentication in the following order:

1. **API key** stored by `mlstdb connect --api-key` (sent as an `X-API-Key` header).
2. **OAuth session token** stored by `mlstdb connect`.
3. **Unauthenticated fallback** if no credentials are found. A warning is printed before sending the request.

This mirrors the same authentication logic used by `mlstdb update`.

## Examples

### Probe a public endpoint without authentication

```sh
mlstdb ping https://rest.pubmlst.org/db --no-auth
```

Expected output:

```
HTTP 200  https://rest.pubmlst.org/db

[
  {
    "description": "Achromobacter spp.",
    "databases": [
      {
        "name": "pubmlst_achromobacter_isolates",
        "description": "Achromobacter spp. isolates",
        "href": "https://rest.pubmlst.org/db/pubmlst_achromobacter_isolates"
      },
      {
        "name": "pubmlst_achromobacter_seqdef",
        "href": "https://rest.pubmlst.org/db/pubmlst_achromobacter_seqdef",
        "description": "Achromobacter spp. sequence/profile definitions"
      }
    ],
    "name": "achromobacter"
  },
    ...
]
```

### Probe an authenticated endpoint using a stored API key

```sh
mlstdb ping https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/schemes/1/profiles_csv --db pubmlst
```

Expected output (first few lines of a profile CSV):

```
HTTP 200  https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/schemes/1/profiles_csv

ST      abcZ    adk     aroE    fumC    gdh     pdhC    pgm     clonal_complex
1       1       3       1       1       1       1       3       ST-1 complex
2       1       3       4       7       1       1       3       ST-1 complex
3       1       3       1       1       1       23      13      ST-1 complex
...
```

### Probe a Pasteur endpoint with verbose output

```sh
mlstdb ping https://bigsdb.pasteur.fr/api/db/pubmlst_bordetella_seqdef/schemes/3 --db pasteur --verbose
```

Expected output:

```
Auth mode:   api-key
Requesting:  https://bigsdb.pasteur.fr/api/db/pubmlst_bordetella_seqdef/schemes/3
Status:      200

HTTP 200  https://bigsdb.pasteur.fr/api/db/pubmlst_bordetella_seqdef/schemes/3

{
  "profiles": "https://bigsdb.pasteur.fr/api/db/pubmlst_bordetella_seqdef/schemes/3/profiles",
  "message": "Please note that you are currently restricted to accessing data that was submitted on or prior to 2024-12-31. Please authenticate to access the full dataset.",
  "description": "MLST",
  "has_primary_key_field": true,
  "loci": [
    "https://bigsdb.pasteur.fr/api/db/pubmlst_bordetella_seqdef/loci/adk",
    "https://bigsdb.pasteur.fr/api/db/pubmlst_bordetella_seqdef/loci/fumC",
    "https://bigsdb.pasteur.fr/api/db/pubmlst_bordetella_seqdef/loci/glyA",
    "https://bigsdb.pasteur.fr/api/db/pubmlst_bordetella_seqdef/loci/tyrB",
    "https://bigsdb.pasteur.fr/api/db/pubmlst_bordetella_seqdef/loci/icd",
    "https://bigsdb.pasteur.fr/api/db/pubmlst_bordetella_seqdef/loci/pepA",
    "https://bigsdb.pasteur.fr/api/db/pubmlst_bordetella_seqdef/loci/pgm"
  ],
  "records": 118,
  ...
}
```

### Unauthenticated fallback when no credentials are stored

If no credentials are stored for the specified database, `ping` falls back to an unauthenticated request and prints a warning:

```sh
mlstdb ping https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/schemes --db pubmlst
```

```
Warning: No credentials found for 'pubmlst', trying unauthenticated...

HTTP 200  https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/schemes

{
  "schemes": [
    {
      "description": "MLST",
      "scheme": "https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/schemes/1"
    },
    ...
  ],
}
```

## Tips

If the server returns a 401 or 403, `ping` prints a short hint:

```
Hint: run 'mlstdb connect --db <db>' to configure credentials, or use --no-auth for unauthenticated access.
```

Use `ping` with `| less` to page through long profile CSV responses:

```sh
mlstdb ping https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/schemes/67/profiles_csv --db pubmlst | less
```
