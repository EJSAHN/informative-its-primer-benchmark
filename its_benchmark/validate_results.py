import argparse
import json
from pathlib import Path
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description='Validate benchmark summaries against an expected-metrics file.')
    parser.add_argument('--summary', required=True)
    parser.add_argument('--expected', required=True)
    parser.add_argument('--out', default='')
    args = parser.parse_args()

    summary = pd.read_csv(args.summary)
    expected = json.loads(Path(args.expected).read_text(encoding='utf-8'))
    checks = []

    totals = {
        'total_rows': int(summary['N'].sum()),
        'informative_rows': int(summary['informative'].sum()),
        'eligible_rows': int(summary['eligible'].sum()),
        'strict_hits': int(summary['covered'].sum()),
    }
    for name, observed in totals.items():
        exp = int(expected[name])
        checks.append({'check': name, 'expected': exp, 'observed': observed, 'passed': exp == observed})

    for record in expected.get('key_pairs', []):
        row = summary[(summary['species'] == record['species']) & (summary['pair_name'] == record['pair_name'])]
        if row.empty:
            checks.append({'check': f"{record['species']} | {record['pair_name']}", 'expected': 'row present', 'observed': 'missing', 'passed': False})
            continue
        row = row.iloc[0]
        for field in ['N', 'informative', 'eligible', 'covered']:
            exp = int(record[field])
            obs = int(row[field])
            checks.append({'check': f"{record['species']} | {record['pair_name']} | {field}", 'expected': exp, 'observed': obs, 'passed': exp == obs})

    result = pd.DataFrame(checks)
    if args.out:
        output = Path(args.out)
        output.parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(output, index=False)
    print(result.to_string(index=False))
    if not bool(result['passed'].all()):
        raise SystemExit(1)
    print('Study-result validation: PASS')


if __name__ == '__main__':
    main()
