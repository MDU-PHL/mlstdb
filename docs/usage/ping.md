# Ping

The `ping` command sends a single GET request to an API endpoint and displays the HTTP status code and response body. It is useful for verifying that your credentials are working or for exploring raw API responses before running `mlstdb update`.

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

{
  "databases": [
    {
      "name": "pubmlst_neisseria_seqdef",
      "description": "Neisseria spp. sequence definitions",
      "href": "https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef"
    },
    ...
  ]
}
```

### Probe an authenticated endpoint using a stored API key

```sh
mlstdb ping https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/schemes/67/profiles_csv --db pubmlst
```

Expected output (first few lines of a profile CSV):

```
HTTP 200  https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/schemes/67/profiles_csv

ST,abcZ,adk,aroE,fumC,gdh,pdhC,pgm,clonal_complex
1,1,3,1,1,1,1,3,ST-1 complex
2,1,3,4,2,5,2,1,ST-32 complex
...
```

### Probe a Pasteur endpoint with verbose output

```sh
mlstdb ping https://bigsdb.pasteur.fr/api/db/pubmlst_bordetella_seqdef/schemes/3 --db pasteur --verbose
```

Expected output:

```
Auth mode:   oauth
Requesting:  https://bigsdb.pasteur.fr/api/db/pubmlst_bordetella_seqdef/schemes/3
Status:      200

HTTP 200  https://bigsdb.pasteur.fr/api/db/pubmlst_bordetella_seqdef/schemes/3

{
  "id": 3,
  "description": "Bordetella MLST",
  "loci": [...],
  "profiles_csv": "https://bigsdb.pasteur.fr/api/db/pubmlst_bordetella_seqdef/schemes/3/profiles_csv"
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
  "schemes": [...]
}
```

## Error Responses

### 401 Unauthorised

A 401 response means the server rejected your credentials. `ping` exits with code 1 and suggests remediation steps:

```
HTTP 401  https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/schemes/67/profiles_csv

Unauthorised

✗ 401 Unauthorised — the server rejected your credentials.

This usually means one of the following:
  1. Your session token has expired — run 'mlstdb connect --db <db>' to refresh.
  2. You are not registered for this scheme/database.
     Visit the database website to register your client application.
  3. Your API key is invalid or revoked — run 'mlstdb connect --db <db> --api-key' to update it.
```

### 403 Forbidden

A 403 response means your account does not have permission to access the resource:

```
HTTP 403  ...

✗ 403 Forbidden — your account does not have permission to access this resource.

This usually means:
  1. You are not registered as a curator or user for this scheme.
  2. Contact the database administrator to request access.
```

## Tips

Use `ping` with `| less` to page through long profile CSV responses:

```sh
mlstdb ping https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/schemes/67/profiles_csv --db pubmlst | less
```
