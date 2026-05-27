import argparse, csv
from pathlib import Path
import pandas as pd

def add_csv(writer, sheet_name, path):
    if not path.exists():
        return False
    df = pd.read_csv(path)
    df.to_excel(writer, sheet_name=sheet_name[:31], index=False)
    ws = writer.sheets[sheet_name[:31]]
    for i, col in enumerate(df.columns):
        width = min(max(10, min(max([len(str(col))] + [len(str(x)) for x in df[col].head(300).fillna('').astype(str)]), 35) + 2), 36)
        ws.set_column(i, i, width)
    ws.freeze_panes(1, 0)
    ws.autofilter(0, 0, len(df), max(0, len(df.columns)-1))
    return True

def main():
    ap = argparse.ArgumentParser(description='Export key CSV outputs as a single Excel workbook.')
    ap.add_argument('--tables-dir', required=True, help='Directory containing output CSV files.')
    ap.add_argument('--main-subdir', default='ncbi_primerlength_metadata', help='Subdirectory containing main summary outputs.')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    tables = Path(args.tables_dir)
    main = tables / args.main_subdir
    out = Path(args.out); out.parent.mkdir(parents=True, exist_ok=True)
    workbook_files = [
        ('Main_summary', main / 'Table_all_pairs_by_species_informative.csv'),
        ('Pair_calls', tables / 'calls_ncbi_primerlength_metadata.csv'),
        ('Site_flags', tables / 'site_flags_ncbi_primerlength_metadata.csv'),
        ('Deduplicated_calls', main / 'calls_deduplicated.csv'),
        ('Scenario_totals', tables / 'scenario_totals.csv'),
        ('Pair_catalog', Path('data') / 'pair_catalog_auto.csv'),
    ]
    with pd.ExcelWriter(out, engine='xlsxwriter') as writer:
        pd.DataFrame({
            'Item': ['Description', 'Main analysis', 'Note'],
            'Value': [
                'Source tables and row-level outputs for the ITS primer benchmark.',
                'NCBI-expanded dataset with primer-length terminal filtering and GenBank metadata flagging.',
                'Large row-level sheets are included for auditability; summary tables provide compact views of these data.'
            ]
        }).to_excel(writer, sheet_name='README', index=False)
        ws = writer.sheets['README']; ws.set_column(0, 0, 22); ws.set_column(1, 1, 95)
        for sheet, path in workbook_files:
            add_csv(writer, sheet, path)
    print(f'Wrote {out}')

if __name__ == '__main__':
    main()
