# Changelog

## [0.1.7] - 2025-11-18

- **License**: Changed from MIT to GPL v3. Original MIT-licensed code is preserved and attributed according to MIT terms.
- Added `get_db_type_from_url()` helper function to determine database type from URL, eliminating code duplication
- Removed redundant `fetch_resources()` function - now using `fetch_json()` directly
- Added acknowledgements section in README.md crediting BIGSdb_downloader and pyMLST projects
- Added CHANGELOG.md file to document version history

[0.1.7]: https://github.com/himal2007/mlstdb/releases/tag/v0.1.7

## [0.2.0] - 2026-01-05

- Updated `database_version.txt` to include a trailing newline character (Closes #11).
- Added `<scheme>_info.json` files for each downloaded scheme containing:
   - Scheme name
   - Number of alleles (locus count)
   - Last updated date
   - Source database (pubmlst/pasteur)
   - API endpoint URL
   - (Closes #20)
- Improved installation instructions in README (Closes #16). 

[0.2.0]: https://github.com/himal2007/mlstdb/releases/tag/v0.2.0