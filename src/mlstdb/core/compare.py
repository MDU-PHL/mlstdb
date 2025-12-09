# SPDX-FileCopyrightText: 2025-present Himal Shrestha <stha.himal2007@gmail.com>
#
# SPDX-License-Identifier: GPL-3.0-or-later

import hashlib
import sqlite3
from pathlib import Path
from typing import Dict, List, Tuple, Set, Optional
from dataclasses import dataclass, asdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
from datetime import datetime


@dataclass
class VersionInfo:
    """Database version information"""
    path: str
    version: str
    
    
@dataclass
class ProfileChange:
    """Represents a change in sequence type profiles"""
    st: str
    change_type: str  # 'added', 'removed', 'modified'
    old_alleles: Optional[Dict[str, str]] = None
    new_alleles: Optional[Dict[str, str]] = None
    
    
@dataclass
class AlleleChange:
    """Represents a change in allele sequences"""
    locus: str
    allele_id: str
    change_type: str  # 'added', 'removed', 'modified'
    old_hash: Optional[str] = None
    new_hash: Optional[str] = None
    old_sequence: Optional[str] = None
    new_sequence: Optional[str] = None


@dataclass
class SchemeComparison:
    """Results of comparing a scheme between two databases"""
    scheme_name: str
    old_version: Optional[str]
    new_version: Optional[str]
    version_changed: bool
    profile_changes: List[ProfileChange]
    allele_changes: List[AlleleChange]
    
    @property
    def profiles_added(self) -> int:
        return sum(1 for p in self.profile_changes if p.change_type == 'added')
    
    @property
    def profiles_removed(self) -> int:
        return sum(1 for p in self.profile_changes if p.change_type == 'removed')
    
    @property
    def profiles_modified(self) -> int:
        return sum(1 for p in self.profile_changes if p.change_type == 'modified')
    
    @property
    def alleles_added(self) -> int:
        return sum(1 for a in self.allele_changes if a.change_type == 'added')
    
    @property
    def alleles_removed(self) -> int:
        return sum(1 for a in self.allele_changes if a.change_type == 'removed')
    
    @property
    def alleles_modified(self) -> int:
        return sum(1 for a in self.allele_changes if a.change_type == 'modified')
    
    @property
    def has_changes(self) -> bool:
        return (self.version_changed or 
                len(self.profile_changes) > 0 or 
                len(self.allele_changes) > 0)


class HashDatabase:
    """SQLite database for storing and comparing allele hashes"""
    
    def __init__(self, db_path: Path):
        self.db_path = db_path
        self.conn = None
        
    def __enter__(self):
        self.conn = sqlite3.connect(self.db_path)
        self._create_tables()
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.conn:
            self.conn.close()
            
    def _create_tables(self):
        """Create tables for storing hash information"""
        cursor = self.conn.cursor()
        
        # Table for allele hashes
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS allele_hashes (
                scheme TEXT NOT NULL,
                locus TEXT NOT NULL,
                allele_id TEXT NOT NULL,
                sequence_hash TEXT NOT NULL,
                sequence_length INTEGER NOT NULL,
                timestamp TEXT NOT NULL,
                PRIMARY KEY (scheme, locus, allele_id)
            )
        """)
        
        # Table for comparison metadata
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS comparison_metadata (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                old_db_path TEXT NOT NULL,
                new_db_path TEXT NOT NULL,
                comparison_date TEXT NOT NULL,
                comparison_type TEXT NOT NULL
            )
        """)
        
        self.conn.commit()
        
    def store_allele_hash(self, scheme: str, locus: str, allele_id: str, 
                         sequence_hash: str, sequence_length: int):
        """Store or update an allele hash"""
        cursor = self.conn.cursor()
        timestamp = datetime.now().isoformat()
        
        cursor.execute("""
            INSERT OR REPLACE INTO allele_hashes 
            (scheme, locus, allele_id, sequence_hash, sequence_length, timestamp)
            VALUES (?, ?, ?, ?, ?, ?)
        """, (scheme, locus, allele_id, sequence_hash, sequence_length, timestamp))
        
        self.conn.commit()
        
    def get_allele_hash(self, scheme: str, locus: str, allele_id: str) -> Optional[Tuple[str, int]]:
        """Retrieve hash and length for an allele"""
        cursor = self.conn.cursor()
        cursor.execute("""
            SELECT sequence_hash, sequence_length 
            FROM allele_hashes 
            WHERE scheme = ? AND locus = ? AND allele_id = ?
        """, (scheme, locus, allele_id))
        
        result = cursor.fetchone()
        return result if result else None
        
    def get_all_alleles_for_scheme(self, scheme: str) -> Dict[Tuple[str, str], Tuple[str, int]]:
        """Get all alleles for a scheme as {(locus, allele_id): (hash, length)}"""
        cursor = self.conn.cursor()
        cursor.execute("""
            SELECT locus, allele_id, sequence_hash, sequence_length
            FROM allele_hashes
            WHERE scheme = ?
        """, (scheme,))
        
        return {(row[0], row[1]): (row[2], row[3]) for row in cursor.fetchall()}
        
    def store_comparison_metadata(self, old_db_path: str, new_db_path: str, 
                                 comparison_type: str):
        """Store metadata about the comparison"""
        cursor = self.conn.cursor()
        timestamp = datetime.now().isoformat()
        
        cursor.execute("""
            INSERT INTO comparison_metadata 
            (old_db_path, new_db_path, comparison_date, comparison_type)
            VALUES (?, ?, ?, ?)
        """, (old_db_path, new_db_path, timestamp, comparison_type))
        
        self.conn.commit()


def compute_sequence_hash(sequence: str) -> str:
    """Compute SHA-256 hash of a sequence (case-insensitive, whitespace-stripped)"""
    # Normalize sequence: uppercase and remove whitespace
    normalized = sequence.upper().strip()
    return hashlib.sha256(normalized.encode()).hexdigest()


def parse_fasta(fasta_path: Path) -> Dict[str, str]:
    """Parse a FASTA file and return {allele_id: sequence}"""
    sequences = {}
    current_id = None
    current_seq = []
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                # Start new sequence
                current_id = line[1:].strip()  # Remove '>'
                current_seq = []
            else:
                current_seq.append(line)
        
        # Save last sequence
        if current_id:
            sequences[current_id] = ''.join(current_seq)
    
    return sequences


def parse_profile(profile_path: Path) -> Tuple[List[str], Dict[str, Dict[str, str]]]:
    """
    Parse a profile file and return (loci_order, {st: {locus: allele}})
    """
    profiles = {}
    loci = []
    
    with open(profile_path, 'r') as f:
        # Read header
        header = f.readline().strip().split('\t')
        if len(header) < 2:
            # Invalid header, return empty
            return [], {}
        loci = header[1:]  # Skip 'ST' column
        
        # Read profiles
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < len(header):
                # Skip malformed lines
                continue
            st = parts[0]
            alleles = {loci[i]: parts[i+1] for i in range(len(loci)) if i+1 < len(parts)}
            profiles[st] = alleles
    
    return loci, profiles


def read_database_version(version_path: Path) -> Optional[str]:
    """Read database version from database_version.txt"""
    if not version_path.exists():
        return None
    
    with open(version_path, 'r') as f:
        return f.read().strip()


def hash_alleles_for_scheme(scheme_path: Path, scheme_name: str, 
                           hash_db: HashDatabase) -> Dict[Tuple[str, str], str]:
    """
    Hash all alleles for a scheme and store in database.
    Returns {(locus, allele_id): hash}
    """
    hashes = {}
    
    # Find all .tfa files
    tfa_files = list(scheme_path.glob('*.tfa'))
    
    for tfa_file in tfa_files:
        locus = tfa_file.stem  # Filename without extension
        sequences = parse_fasta(tfa_file)
        
        for allele_id, sequence in sequences.items():
            seq_hash = compute_sequence_hash(sequence)
            hash_db.store_allele_hash(
                scheme_name, locus, allele_id, seq_hash, len(sequence)
            )
            hashes[(locus, allele_id)] = seq_hash
    
    return hashes


def compare_profiles(old_profiles: Dict[str, Dict[str, str]], 
                    new_profiles: Dict[str, Dict[str, str]]) -> List[ProfileChange]:
    """Compare two sets of profiles and identify changes"""
    changes = []
    
    old_sts = set(old_profiles.keys())
    new_sts = set(new_profiles.keys())
    
    # Added STs
    for st in new_sts - old_sts:
        changes.append(ProfileChange(
            st=st,
            change_type='added',
            new_alleles=new_profiles[st]
        ))
    
    # Removed STs
    for st in old_sts - new_sts:
        changes.append(ProfileChange(
            st=st,
            change_type='removed',
            old_alleles=old_profiles[st]
        ))
    
    # Modified STs
    for st in old_sts & new_sts:
        if old_profiles[st] != new_profiles[st]:
            changes.append(ProfileChange(
                st=st,
                change_type='modified',
                old_alleles=old_profiles[st],
                new_alleles=new_profiles[st]
            ))
    
    return changes


def compare_alleles(old_path: Path, new_path: Path, 
                   scheme_name: str) -> List[AlleleChange]:
    """Compare alleles between two scheme directories"""
    changes = []
    
    # Get all .tfa files from both directories
    old_tfa_files = {f.stem: f for f in old_path.glob('*.tfa')}
    new_tfa_files = {f.stem: f for f in new_path.glob('*.tfa')}
    
    all_loci = set(old_tfa_files.keys()) | set(new_tfa_files.keys())
    
    for locus in all_loci:
        old_sequences = parse_fasta(old_tfa_files[locus]) if locus in old_tfa_files else {}
        new_sequences = parse_fasta(new_tfa_files[locus]) if locus in new_tfa_files else {}
        
        old_alleles = set(old_sequences.keys())
        new_alleles = set(new_sequences.keys())
        
        # Added alleles
        for allele_id in new_alleles - old_alleles:
            seq = new_sequences[allele_id]
            changes.append(AlleleChange(
                locus=locus,
                allele_id=allele_id,
                change_type='added',
                new_hash=compute_sequence_hash(seq),
                new_sequence=seq
            ))
        
        # Removed alleles
        for allele_id in old_alleles - new_alleles:
            seq = old_sequences[allele_id]
            changes.append(AlleleChange(
                locus=locus,
                allele_id=allele_id,
                change_type='removed',
                old_hash=compute_sequence_hash(seq),
                old_sequence=seq
            ))
        
        # Check for modified alleles (same ID, different sequence)
        for allele_id in old_alleles & new_alleles:
            old_seq = old_sequences[allele_id]
            new_seq = new_sequences[allele_id]
            old_hash = compute_sequence_hash(old_seq)
            new_hash = compute_sequence_hash(new_seq)
            
            if old_hash != new_hash:
                changes.append(AlleleChange(
                    locus=locus,
                    allele_id=allele_id,
                    change_type='modified',
                    old_hash=old_hash,
                    new_hash=new_hash,
                    old_sequence=old_seq,
                    new_sequence=new_seq
                ))
    
    return changes


def compare_scheme(old_db_path: Path, new_db_path: Path, 
                  scheme_name: str, comparison_type: str = 'basic') -> SchemeComparison:
    """
    Compare a single scheme between old and new databases.
    
    Args:
        old_db_path: Path to old database root
        new_db_path: Path to new database root
        scheme_name: Name of the scheme to compare
        comparison_type: 'basic' or 'full'
    
    Returns:
        SchemeComparison object with all changes
    """
    old_scheme = old_db_path / scheme_name
    new_scheme = new_db_path / scheme_name
    
    # Check if schemes exist
    if not old_scheme.exists() and not new_scheme.exists():
        raise ValueError(f"Scheme '{scheme_name}' not found in either database")
    
    # Read versions
    old_version = read_database_version(old_scheme / 'database_version.txt') if old_scheme.exists() else None
    new_version = read_database_version(new_scheme / 'database_version.txt') if new_scheme.exists() else None
    version_changed = old_version != new_version
    
    # Find profile files
    old_profile_file = old_scheme / f"{scheme_name}.txt" if old_scheme.exists() else None
    new_profile_file = new_scheme / f"{scheme_name}.txt" if new_scheme.exists() else None
    
    # Handle missing profile files
    if old_profile_file and not old_profile_file.exists():
        old_profile_file = None
    if new_profile_file and not new_profile_file.exists():
        new_profile_file = None
    
    # Compare profiles
    profile_changes = []
    if old_profile_file and new_profile_file:
        _, old_profiles = parse_profile(old_profile_file)
        _, new_profiles = parse_profile(new_profile_file)
        profile_changes = compare_profiles(old_profiles, new_profiles)
    elif new_profile_file:
        # All profiles are new
        _, new_profiles = parse_profile(new_profile_file)
        profile_changes = [
            ProfileChange(st=st, change_type='added', new_alleles=alleles)
            for st, alleles in new_profiles.items()
        ]
    elif old_profile_file:
        # All profiles removed
        _, old_profiles = parse_profile(old_profile_file)
        profile_changes = [
            ProfileChange(st=st, change_type='removed', old_alleles=alleles)
            for st, alleles in old_profiles.items()
        ]
    
    # Compare alleles (only for 'full' comparison)
    allele_changes = []
    if comparison_type == 'full':
        if old_scheme.exists() and new_scheme.exists():
            allele_changes = compare_alleles(old_scheme, new_scheme, scheme_name)
    
    return SchemeComparison(
        scheme_name=scheme_name,
        old_version=old_version,
        new_version=new_version,
        version_changed=version_changed,
        profile_changes=profile_changes,
        allele_changes=allele_changes
    )


def _compare_scheme_worker(args: Tuple[Path, Path, str, str]) -> SchemeComparison:
    """Worker function for parallel scheme comparison"""
    old_db_path, new_db_path, scheme_name, comparison_type = args
    return compare_scheme(old_db_path, new_db_path, scheme_name, comparison_type)


def compare_databases(old_db_path: Path, new_db_path: Path, 
                     comparison_type: str = 'basic',
                     parallel: bool = True,
                     max_workers: Optional[int] = None) -> List[SchemeComparison]:
    """
    Compare all schemes between two databases.
    
    Args:
        old_db_path: Path to old database root
        new_db_path: Path to new database root
        comparison_type: 'basic' or 'full'
        parallel: Whether to use parallel processing
        max_workers: Maximum number of parallel workers (None = CPU count)
    
    Returns:
        List of SchemeComparison objects
    """
    # Find all schemes (directories in the database path)
    old_schemes = {p.name for p in old_db_path.iterdir() if p.is_dir()}
    new_schemes = {p.name for p in new_db_path.iterdir() if p.is_dir()}
    all_schemes = sorted(old_schemes | new_schemes)
    
    if not all_schemes:
        return []
    
    results = []
    
    if parallel and len(all_schemes) > 1:
        # Parallel processing
        args_list = [
            (old_db_path, new_db_path, scheme, comparison_type)
            for scheme in all_schemes
        ]
        
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {executor.submit(_compare_scheme_worker, args): args[2] 
                      for args in args_list}
            
            for future in as_completed(futures):
                try:
                    result = future.result()
                    results.append(result)
                except Exception as e:
                    scheme_name = futures[future]
                    print(f"Error comparing scheme {scheme_name}: {e}")
    else:
        # Sequential processing
        for scheme in all_schemes:
            try:
                result = compare_scheme(old_db_path, new_db_path, scheme, comparison_type)
                results.append(result)
            except Exception as e:
                print(f"Error comparing scheme {scheme}: {e}")
    
    return results


def export_comparison_json(comparisons: List[SchemeComparison], 
                          output_path: Path):
    """Export comparison results to JSON"""
    data = {
        'timestamp': datetime.now().isoformat(),
        'schemes': []
    }
    
    for comp in comparisons:
        scheme_data = {
            'scheme_name': comp.scheme_name,
            'old_version': comp.old_version,
            'new_version': comp.new_version,
            'version_changed': comp.version_changed,
            'summary': {
                'profiles_added': comp.profiles_added,
                'profiles_removed': comp.profiles_removed,
                'profiles_modified': comp.profiles_modified,
                'alleles_added': comp.alleles_added,
                'alleles_removed': comp.alleles_removed,
                'alleles_modified': comp.alleles_modified
            },
            'profile_changes': [asdict(pc) for pc in comp.profile_changes],
            'allele_changes': [asdict(ac) for ac in comp.allele_changes]
        }
        data['schemes'].append(scheme_data)
    
    with open(output_path, 'w') as f:
        json.dump(data, f, indent=2)


def export_comparison_text(comparisons: List[SchemeComparison], 
                          output_path: Path):
    """Export comparison results to human-readable text"""
    with open(output_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("MLST Database Comparison Report\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("=" * 80 + "\n\n")
        
        # Summary
        total_schemes = len(comparisons)
        schemes_with_changes = sum(1 for c in comparisons if c.has_changes)
        total_profiles_added = sum(c.profiles_added for c in comparisons)
        total_profiles_removed = sum(c.profiles_removed for c in comparisons)
        total_profiles_modified = sum(c.profiles_modified for c in comparisons)
        total_alleles_added = sum(c.alleles_added for c in comparisons)
        total_alleles_removed = sum(c.alleles_removed for c in comparisons)
        total_alleles_modified = sum(c.alleles_modified for c in comparisons)
        
        f.write("SUMMARY\n")
        f.write("-" * 80 + "\n")
        f.write(f"Total schemes compared: {total_schemes}\n")
        f.write(f"Schemes with changes: {schemes_with_changes}\n")
        f.write(f"\nProfile Changes:\n")
        f.write(f"  Added: {total_profiles_added}\n")
        f.write(f"  Removed: {total_profiles_removed}\n")
        f.write(f"  Modified: {total_profiles_modified}\n")
        
        if total_alleles_added + total_alleles_removed + total_alleles_modified > 0:
            f.write(f"\nAllele Changes:\n")
            f.write(f"  Added: {total_alleles_added}\n")
            f.write(f"  Removed: {total_alleles_removed}\n")
            f.write(f"  Modified: {total_alleles_modified}\n")
        
        f.write("\n" + "=" * 80 + "\n\n")
        
        # Detailed results per scheme
        for comp in comparisons:
            if not comp.has_changes:
                continue
                
            f.write(f"SCHEME: {comp.scheme_name}\n")
            f.write("-" * 80 + "\n")
            
            # Version info
            if comp.version_changed:
                f.write(f"Database Version: {comp.old_version} → {comp.new_version}\n")
            else:
                f.write(f"Database Version: {comp.old_version or 'N/A'}\n")
            
            # Profile changes
            if comp.profile_changes:
                f.write(f"\nProfile Changes:\n")
                f.write(f"  Added: {comp.profiles_added} STs\n")
                f.write(f"  Removed: {comp.profiles_removed} STs\n")
                f.write(f"  Modified: {comp.profiles_modified} STs\n")
                
                # Show details for removed and modified
                removed = [p for p in comp.profile_changes if p.change_type == 'removed']
                if removed:
                    f.write(f"\n  Removed STs: {', '.join(p.st for p in removed[:20])}")
                    if len(removed) > 20:
                        f.write(f" ... and {len(removed) - 20} more")
                    f.write("\n")
                
                modified = [p for p in comp.profile_changes if p.change_type == 'modified']
                if modified:
                    f.write(f"\n  Modified STs:\n")
                    for p in modified[:10]:  # Show first 10
                        f.write(f"    ST {p.st}: ")
                        changes = []
                        for locus in p.old_alleles.keys():
                            if p.old_alleles.get(locus) != p.new_alleles.get(locus):
                                changes.append(f"{locus}:{p.old_alleles[locus]}→{p.new_alleles[locus]}")
                        f.write(", ".join(changes) + "\n")
                    if len(modified) > 10:
                        f.write(f"    ... and {len(modified) - 10} more modified STs\n")
            
            # Allele changes
            if comp.allele_changes:
                f.write(f"\nAllele Changes:\n")
                f.write(f"  Added: {comp.alleles_added} alleles\n")
                f.write(f"  Removed: {comp.alleles_removed} alleles\n")
                f.write(f"  Modified: {comp.alleles_modified} alleles\n")
                
                # Show modified alleles (sequence changes are critical)
                modified_alleles = [a for a in comp.allele_changes if a.change_type == 'modified']
                if modified_alleles:
                    f.write(f"\n  WARNING: Modified alleles (sequence changed):\n")
                    for a in modified_alleles[:10]:
                        f.write(f"    {a.locus}/{a.allele_id}\n")
                        f.write(f"      Old hash: {a.old_hash}\n")
                        f.write(f"      New hash: {a.new_hash}\n")
                    if len(modified_alleles) > 10:
                        f.write(f"    ... and {len(modified_alleles) - 10} more modified alleles\n")
            
            f.write("\n" + "=" * 80 + "\n\n")
