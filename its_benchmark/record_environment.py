import argparse
import hashlib
import importlib.metadata as metadata
import platform
import sys
from pathlib import Path
import pandas as pd

PACKAGES = ['biopython', 'pandas', 'numpy', 'xlsxwriter', 'requests']


def main():
    parser = argparse.ArgumentParser(description='Record software versions and SHA-256 checksums for analysis outputs.')
    parser.add_argument('--results-dir', required=True)
    args = parser.parse_args()

    results = Path(args.results_dir)
    results.mkdir(parents=True, exist_ok=True)
    rows = [
        {'component': 'Python', 'version': sys.version.split()[0]},
        {'component': 'Platform', 'version': platform.platform()},
    ]
    for package in PACKAGES:
        try:
            version = metadata.version(package)
        except metadata.PackageNotFoundError:
            version = 'NOT_INSTALLED'
        rows.append({'component': package, 'version': version})
    pd.DataFrame(rows).to_csv(results / 'software_versions.csv', index=False)

    checksums = []
    for path in sorted(results.rglob('*')):
        if not path.is_file() or path.name == 'sha256_manifest.tsv':
            continue
        checksums.append({
            'file': str(path.relative_to(results)),
            'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
            'bytes': path.stat().st_size,
        })
    pd.DataFrame(checksums).to_csv(results / 'sha256_manifest.tsv', sep='\t', index=False)
    print(f'Wrote environment and checksum records to {results}')


if __name__ == '__main__':
    main()
