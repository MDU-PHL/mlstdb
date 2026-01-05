# Changelog

## [0.1.7] - 2025-11-18

- **License**: Changed from MIT to GPL v3. Original MIT-licensed code is preserved and attributed according to MIT terms.
- Added `get_db_type_from_url()` helper function to determine database type from URL, eliminating code duplication
- Removed redundant `fetch_resources()` function - now using `fetch_json()` directly
- Added acknowledgements section in README.md crediting BIGSdb_downloader and pyMLST projects
- Added CHANGELOG.md file to document version history

[0.1.7]: https://github.com/himal2007/mlstdb/releases/tag/v0.1.7

## [0.2.0] - 2026-01-05

### Added
- Scheme metadata JSON file for each downloaded scheme ([#11])
- Newline character to `database_version. txt` for Unix tool compatibility ([#20])

### Fixed
- Dependency installation issues when using bioconda ([#16])

### Changed
- Installation instructions to recommend conda-forge channel and pip installation method


[0.1.7]: https://github.com/MDU-PHL/mlstdb/releases/tag/v0.1.7
[0.2.0]: https://github.com/MDU-PHL/mlstdb/releases/tag/v0.2.0
[#11]: https://github.com/MDU-PHL/mlstdb/issues/11
[#16]: https://github.com/MDU-PHL/mlstdb/issues/16
[#20]: https://github.com/MDU-PHL/mlstdb/issues/20