# Changelog

## [1.2.0] - 2026-05-22

### Added

- Personal API key support (BIGSdb ≥ v1.53.0): `mlstdb connect --api-key` registers a key generated from your BIGSdb profile page. The key is stored in `~/.config/mlstdb/api_keys` (`0600`) and used automatically by `update` and `fetch` in preference to OAuth tokens.
- New `setup_api_key()` and `retrieve_api_key()` functions in `mlstdb.core.auth` for storing and retrieving personal API keys.
- `create_session()` in `mlstdb.core.download` now accepts an `api_key` parameter; when set it creates a plain `requests.Session` with the `X-API-Key` header.
- `remove_db_credentials()` now also removes stored API keys when purging a site's credentials.

### Fixed

- `TypeError: argument of type 'NoneType' is not iterable` in `mlstdb connect --db pubmlst` on Rocky Linux and WSL2 ([#35](https://github.com/MDU-PHL/mlstdb/issues/35)). Root cause: rauth 0.7.3 `_parse_optional_params` calls `req_kwargs.get('params')` without a default, raising `TypeError` when `params` is absent from request kwargs. Fixed by passing `params={}` to every bare `OAuth1Session.get()` call across `auth.py` and `download.py`.
- Unsafe `r.json()['message']` error parsing in `register_tokens()` and throughout the download pipeline replaced by a new `_parse_error_message()` helper, which gracefully handles non-JSON responses (HTML error pages, empty bodies) returned by BIGSdb on
  HTTPS redirects or proxy errors.
- `mlstdb connect` now prints actionable clock-sync instructions when the OAuth request token step fails with a `timestamp more than 600 seconds` error.


## [1.1.1] - 2026-03-31

### Fixed

- `check_dir` now uses an actual write test instead of `os.access`, fixing false "Cannot write to directory" errors on NFS-mounted filesystems common in HPC environments. Thanks to @talasjudit for the detailed report and suggested fix. ([#32](https://github.com/MDU-PHL/mlstdb/issues/32))
- `last-updated` field in scheme info JSON is now capped at `2024-12-31` when running with `--no-auth`, accurately reflecting the data cutoff date for unauthenticated access. ([#31](https://github.com/MDU-PHL/mlstdb/issues/31))
- Removed the legacy `database_version.txt` file from the scheme directory. Scheme metadata is now stored exclusively in `{scheme}_info.json`. ([#11](https://github.com/MDU-PHL/mlstdb/issues/11))

## [1.1.0] - 2026-03-24

### Added

- New `mlstdb purge` command for removing schemes, STs, or individual alleles
  from the local database, with automatic BLAST rebuild afterwards.
  - Purge an entire scheme: `mlstdb purge -s salmonella`
  - Purge a specific ST: `mlstdb purge -s salmonella --st 3`
  - Purge a specific allele: `mlstdb purge -s salmonella -a aroC:1`
  - Batch purge across multiple schemes from a YAML config file: `mlstdb purge -c purge_config.yaml`
  - Before removing a ST, checks whether its alleles are used by other STs and
    warns if so — use `--force` to override
  - Before removing an allele, lists all STs that will be affected and prompts
    for confirmation
  - BLAST database is rebuilt once at the end, minimising redundant work
  - Supports `--force` to skip all confirmation prompts, and `--verbose` for
    detailed logging

### Fixed

- Incomplete scheme directories (missing profiles or allele files) are now removed
  before BLAST database creation when using `--no-auth` or when authentication
  failures leave partial downloads on disk. A warning lists the affected schemes
  and advises re-running with authenticated access. ([#29](https://github.com/MDU-PHL/mlstdb/issues/29))
- Corrected JSON key in `{scheme_name}_info.json` from `"locus"` to `"locii"`
  for accurate scheme metadata. ([#28](https://github.com/MDU-PHL/mlstdb/issues/28))

## [1.0.0] - 2026-03-13

### Added
- New `mlstdb connect` command for streamlined OAuth credential registration ([#25](https://github.com/MDU-PHL/mlstdb/issues/25))
- Curated built-in scheme list (`mlst_schemes_all.tab`) — `mlstdb update` works out of the box without `fetch` ([#10](https://github.com/MDU-PHL/mlstdb/issues/10))
- `--no-auth` flag for unauthenticated access to public APIs on both `fetch` and `update`
- `--resume` flag for `update` to skip already-downloaded schemes
- `--threads` option for parallel downloads on both `fetch` and `update`
- Session reuse and HTTP connection pooling for improved performance
- Restrictive file permissions (`0600`) on stored credential files
- Comprehensive MkDocs documentation site with detailed guides for all commands

### Fixed
- Fetch looping error at 76% when processing databases ([#19](https://github.com/MDU-PHL/mlstdb/issues/19))
- Missing scheme URI resolution errors during fetch ([#18](https://github.com/MDU-PHL/mlstdb/issues/18))
- 401 errors on unregistered databases no longer terminate the process — skipped databases are reported at the end ([#17](https://github.com/MDU-PHL/mlstdb/issues/17))
- Automatic token refresh on expired session tokens

### Changed
- `fetch` command deprecated in favour of `connect` + `update` workflow ([#25](https://github.com/MDU-PHL/mlstdb/issues/25))
- `update` now uses the built-in curated scheme list by default (no `--input` required)
- Simplified README focused on the two-command workflow

## [0.2.0] - 2026-01-05

### Added
- Scheme metadata JSON file for each downloaded scheme ([#11])
- Newline character to `database_version. txt` for Unix tool compatibility ([#20])

### Fixed
- Dependency installation issues when using bioconda ([#16])

### Changed
- Installation instructions to recommend conda-forge channel and pip installation method

## [0.1.7] - 2025-11-18

### Changed
- **License**: Changed from MIT to GPL v3. Original MIT-licensed code is preserved and attributed according to MIT terms.

### Added
- `get_db_type_from_url()` helper function to determine database type from URL
- Acknowledgements section in README.md crediting BIGSdb_downloader and pyMLST projects
- CHANGELOG.md file

### Improved

- Removed redundant `fetch_resources()` function — now using `fetch_json()` directly


[1.0.0]: https://github.com/MDU-PHL/mlstdb/releases/tag/v1.0.0
[0.2.0]: https://github.com/MDU-PHL/mlstdb/releases/tag/v0.2.0
[0.1.7]: https://github.com/MDU-PHL/mlstdb/releases/tag/v0.1.7
[1.1.0]: https://github.com/MDU-PHL/mlstdb/releases/tag/v1.1.0
[1.1.1]: https://github.com/MDU-PHL/mlstdb/releases/tag/v1.1.1
[1.2.0]: https://github.com/MDU-PHL/mlstdb/releases/tag/v1.2.0
[#11]: https://github.com/MDU-PHL/mlstdb/issues/11
[#16]: https://github.com/MDU-PHL/mlstdb/issues/16
[#20]: https://github.com/MDU-PHL/mlstdb/issues/20