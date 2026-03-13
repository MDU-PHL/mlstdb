# Changelog

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
[#11]: https://github.com/MDU-PHL/mlstdb/issues/11
[#16]: https://github.com/MDU-PHL/mlstdb/issues/16
[#20]: https://github.com/MDU-PHL/mlstdb/issues/20