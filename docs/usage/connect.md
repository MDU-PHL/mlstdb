# Connect

The `connect` command registers your credentials with PubMLST or Pasteur databases. This is a **one-time setup** per database. Credentials are saved locally and reused automatically.

Two authentication methods are supported:

| Method | When to use |
|--------|-------------|
| **Personal API Key** *(recommended)* | BIGSdb ≥ v1.53.0 (PubMLST since May 2026). Simpler setup, no OAuth flow, no browser required. Note, only available for PubMLST. |
| **OAuth** *(legacy)* | All BIGSdb versions. Required for Pasteur or any older instance. |

---

## Basic Usage

```sh
# Personal API Key (recommended for PubMLST)
mlstdb connect --db pubmlst --api-key

# OAuth (legacy, or for Pasteur)
mlstdb connect --db pubmlst
mlstdb connect --db pasteur
```

If you omit `--db`, you'll be prompted to choose.

---

## Options

| Option | Description |
|--------|-------------|
| `--db`, `-d` | Database to connect to: `pubmlst` or `pasteur` |
| `--api-key` | Register using a personal API key instead of OAuth (BIGSdb ≥ v1.53.0) |
| `--verbose`, `-v` | Show detailed debug output |
| `-h`, `--help` | Show help message |

---

## Personal API Key

BIGSdb v1.53.0 introduced personal API keys as a simpler alternative to OAuth. A single key is generated from your BIGSdb user profile and passed via the `X-API-Key` request header. There is no browser step, no verification code, and no OAuth token expiry.

### Obtaining a Personal API Key

#### PubMLST

1. Go to [https://pubmlst.org/bigsdb](https://pubmlst.org/bigsdb) and log in
2. Open your **My Account** (top-right menu)
3. Scroll to the **API keys** section and generate a new key for MLST database access (if not already generated).
4. Copy the key.

#### Pasteur

Personal API keys require BIGSdb ≥ v1.53.0. Check with the Pasteur team whether their instance has been upgraded before using `--api-key` for Pasteur.

### Registering your API Key

```sh
mlstdb connect --db pubmlst --api-key
```

You will be prompted:

```
Please enter your personal API key (from your BIGSdb profile page):
API Key: **********************
```

The key is saved to `~/.config/mlstdb/api_keys` with `0600` permissions and used automatically by all subsequent `mlstdb update` and `mlstdb fetch` calls.

!!! tip
    You can regenerate your key at any time from your BIGSdb profile. Run `mlstdb connect --db pubmlst --api-key` again to update the stored key.

---

## OAuth (Legacy)

Use OAuth if you are on an older BIGSdb instance or need access to private data or submission endpoints (API keys only cover public data reads).

### Obtaining Client Credentials

Before running `mlstdb connect`, you need to register as an API client on each database platform to get your **Client ID** and **Client Secret**.

### PubMLST

1. Go to [https://pubmlst.org/bigsdb](https://pubmlst.org/bigsdb)
2. Create an account or log in
3. Under the "Database registrations" section, select "check all" the databases you want access to. You may register for just the ones you need, but you might not be able to access schemes from databases you haven't registered for. You will get an warning while running `mlstdb update` if you haven't registered for a scheme's database.
4. Make sure you register to (e.g., *Neisseria* → `pubmlst_neisseria_seqdef`) as `mlstdb` uses this for connecting to the API and downloading schemes. 
5. Go to **API keys** section and create a new API key for MLST database access. 
6. Copy the **Client ID** (24 characters) and **Client Secret** (42 characters)

### Pasteur (BIGSdb)

1. Go to [https://bigsdb.pasteur.fr/cgi-bin/bigsdb/bigsdb.pl](https://bigsdb.pasteur.fr/cgi-bin/bigsdb/bigsdb.pl)
2. Same as for PubMLST, create an account or log in and register for all the databases you want access to.
3. You will need to email the Pasteur team to request API access and get your Client ID and Client Secret.
4. Use the **Client ID** and **Client Secret** for `mlstdb connect --db pasteur`

!!! tip
    You only need to register your client credentials **once** per platform. The same credentials work across all databases within that platform.

!!! important "Register with individual databases"
    While the client credentials are platform-wide, you may need to **register/enrol** within specific databases on each platform to access their schemes. If you get authentication errors during `mlstdb update`, check that you've registered with the relevant databases.

---

## What Happens During `connect`

### Personal API Key flow

1. **You provide** your personal API key (from your BIGSdb profile)
2. **mlstdb** tests the key against the database API
3. **Key is saved** to `~/.config/mlstdb/api_keys`

### OAuth flow

1. **You provide** your Client ID and Client Secret
2. **mlstdb** requests a temporary token from the API
3. **A URL is displayed** — open it in your browser
4. **Authorise** the application on the website
5. **Enter the verification code** shown on the website back into the terminal
6. **mlstdb** exchanges this for access and session tokens
7. **Credentials are saved** to `~/.config/mlstdb/`

```
~/.config/mlstdb/
├── api_keys              # Personal API key (BIGSdb ≥ v1.53.0)
├── client_credentials    # OAuth Client ID and Secret
├── access_tokens         # OAuth access tokens
└── session_tokens        # OAuth session tokens (used for API calls)
```

!!! note "Security"
    All credential files are stored with restrictive permissions (`0700`). They are only readable by your user account.
---

## Re-connecting

### Personal API Key

If your key is revoked or regenerated on the BIGSdb side, re-run:

```sh
mlstdb connect --db pubmlst --api-key
```

You will be prompted whether to replace the existing key.

### OAuth

If your OAuth credentials expire or become invalid, `mlstdb connect` will detect this:

```sh
mlstdb connect --db pubmlst
# Output: ✓ Credentials found for pubmlst
# Output: ✗ Connection test failed for pubmlst
# Prompt: Do you want to re-register with pubmlst?
```

Answer **yes** to go through the registration flow again.

---

## Environment Variables

For headless execution, CI/CD pipelines, containerised environments, or HPC batch jobs, credentials can be supplied using database-specific environment variables:

| Credential Type | Database | Environment Variable |
| :--- | :--- | :--- |
| **Personal API Key** | PubMLST | `MLSTDB_PUBMLST_API_KEY` |
| **Personal API Key** | Pasteur | `MLSTDB_PASTEUR_API_KEY` |
| **OAuth Client ID** | PubMLST | `MLSTDB_PUBMLST_CLIENT_ID` |
| **OAuth Client Secret** | PubMLST | `MLSTDB_PUBMLST_CLIENT_SECRET` |
| **OAuth Client ID** | Pasteur | `MLSTDB_PASTEUR_CLIENT_ID` |
| **OAuth Client Secret** | Pasteur | `MLSTDB_PASTEUR_CLIENT_SECRET` |

### Precedence

`mlstdb` evaluates credentials in the following order:
1. **Environment Variables** (e.g., `MLSTDB_PUBMLST_API_KEY` or `MLSTDB_PUBMLST_CLIENT_ID` / `MLSTDB_PUBMLST_CLIENT_SECRET`)
2. **Local Configuration Files** (`~/.config/mlstdb/api_keys` or `~/.config/mlstdb/client_credentials`)
3. **Interactive Terminal Prompt** (if credentials are required and not found)



