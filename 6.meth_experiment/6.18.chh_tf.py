#!/usr/bin/env python3
"""Summarize direct and <=50 bp nearby FIMO hits for six CHH sites.

BED sites are 0-based half-open; FIMO coordinates are 1-based closed.
The non-coexpressed group is the unique TF set in coex2/total/all.TF.txt
minus the coexpressed TF set. Annotation availability never limits this set.
Both FIMO strands are included. Per-site/TF motif and flanking columns are
binary indicators: any matching hit counts as one regardless of strand.
No additional significance filtering is applied. Python 3; no dependencies.
"""

import argparse
import csv
from pathlib import Path


DEFAULT_WORKDIR = Path(
    '/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/'
    'Arabidopsis_thaliana/list'
)
SITE_IDS = ('CHH_1634', 'CHH_1635', 'CHH_1651', 'CHH_1676', 'CHH_1698', 'CHH_1771')
CATEGORIES = ('直接落入', '邻近50bp')


def load_tf_symbols(path):
    """Collapse repeated annotation rows without multiplying motif hits."""
    symbols = {}
    with path.open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not {'Gene_ID', 'Symbol'}.issubset(reader.fieldnames or []):
            raise ValueError(f'{path}: 缺少 Gene_ID 或 Symbol 列')
        for line, row in enumerate(reader, 2):
            gene = (row.get('Gene_ID') or '').strip()
            symbol = (row.get('Symbol') or '').strip()
            if not gene:
                raise ValueError(f'{path}:{line}: Gene_ID 为空')
            symbols.setdefault(gene, set())
            if symbol and symbol.upper() != 'NA':
                symbols[gene].add(symbol)
    return {gene: '/'.join(sorted(names)) or 'NA' for gene, names in symbols.items()}


def load_sites(path):
    sites = {}
    with path.open(encoding='utf-8-sig') as handle:
        for line, text in enumerate(handle, 1):
            fields = text.split()
            if not fields or fields[0] not in SITE_IDS:
                continue
            site = fields[0]
            if site in sites:
                raise ValueError(f'{path}:{line}: 候选位点重复: {site}')
            if len(fields) < 4:
                raise ValueError(f'{path}:{line}: 位点坐标列不足')
            try:
                start, end = int(fields[2]), int(fields[3])
            except ValueError as exc:
                raise ValueError(f'{path}:{line}: 位点坐标必须为整数') from exc
            if start < 0 or end != start + 1:
                raise ValueError(f'{path}:{line}: 候选位点须为合法的单碱基 BED 区间')
            sites[site] = (fields[1], start, end)
    missing = [site for site in SITE_IDS if site not in sites]
    if missing:
        raise ValueError(f'{path}: 缺少候选位点: {", ".join(missing)}')
    return sites


def load_groups(workdir):
    """Define membership by set difference; optional annotations only add names."""
    all_tf_ids = {hit[0] for hit in load_motifs(workdir / 'coex2/total/all.TF.txt')}
    coexpressed = load_tf_symbols(
        workdir / 'coex_experiment/total.gene.TF_ids.annotated.txt'
    )
    symbols = {}
    for relative in (
        'coex_experiment/total.other.gene.TF_ids.annotated.txt',
        'coex2/total/total.TF_in_module.annotated.txt',
        'coex2/total/total.other.gene.TF_ids.annotated.txt',
    ):
        path = workdir / relative
        if path.is_file():
            for gene, symbol in load_tf_symbols(path).items():
                if symbol != 'NA':
                    symbols[gene] = symbol
    other = {gene: symbols.get(gene, 'NA') for gene in sorted(all_tf_ids - coexpressed.keys())}
    return {'共表达': coexpressed, '非共表达': other}


def load_motifs(path):
    """Read FIMO with or without a header; deduplicate gene/region/strand hits."""
    motifs = set()
    with path.open(encoding='utf-8-sig') as handle:
        for line, text in enumerate(handle, 1):
            fields = text.split()
            if not fields or text.lstrip().startswith('#') or fields[0] == 'motif_id':
                continue
            if len(fields) != 10:
                raise ValueError(f'{path}:{line}: FIMO 记录应有 10 列')
            gene, _, chrom, start, end, strand = fields[:6]
            try:
                start, end = int(start), int(end)
            except ValueError as exc:
                raise ValueError(f'{path}:{line}: FIMO 坐标必须为整数') from exc
            if start < 1 or end < start or strand not in ('+', '-'):
                raise ValueError(f'{path}:{line}: FIMO 坐标或链方向非法')
            motifs.add((gene, chrom, start, end, strand))
    return sorted(motifs)


def classify_hit(site_chrom, position, motif_chrom, start, end):
    """Nearest nucleotide distance, with inclusive boundaries on both sides."""
    if site_chrom != motif_chrom:
        return None
    distance = max(start - position, position - end, 0)
    if distance == 0:
        return CATEGORIES[0], distance
    if distance <= 50:
        return CATEGORIES[1], distance
    return None


def summarize(sites, groups, motifs):
    rows = []
    for site in SITE_IDS:
        chrom, bed_start, bed_end = sites[site]
        position = bed_start + 1
        row = {
            '位点ID': site,
            '染色体': chrom,
            'BED起点_0based': bed_start,
            'BED终点_不含': bed_end,
            '位点坐标_1based': position,
        }
        for group, symbols in groups.items():
            hits = {category: [] for category in CATEGORIES}
            genes = {category: set() for category in CATEGORIES}
            for gene, motif_chrom, start, end, strand in motifs:
                if gene not in symbols:
                    continue
                match = classify_hit(chrom, position, motif_chrom, start, end)
                if match is None:
                    continue
                category, distance = match
                genes[category].add(gene)
                hits[category].append(
                    f'{gene}|{symbols[gene]}|{motif_chrom}:{start}-{end}'
                    f'(1-based闭区间)|链={strand}|距离={distance}bp'
                )
            for category in CATEGORIES:
                prefix = f'{group}{category}'
                row[f'{prefix}_TF数量'] = len(genes[category])
                row[f'{prefix}_匹配详情'] = '; '.join(hits[category]) or '无'
        rows.append(row)
    return rows


def summarize_tf_counts(sites, symbols, motifs):
    """One row per matched site/TF, with binary motif/flanking indicators.

    Omit unmatched sites/TFs. Positive-only,
    negative-only and both-strand matches each count once per category.
    """
    rows = []
    for site in SITE_IDS:
        chrom, start, _ = sites[site]
        position = start + 1
        counts = {}
        for gene, motif_chrom, motif_start, motif_end, strand in motifs:
            if gene not in symbols:
                continue
            match = classify_hit(chrom, position, motif_chrom, motif_start, motif_end)
            if match is None:
                continue
            values = counts.setdefault(gene, [0, 0])
            values[0 if match[0] == CATEGORIES[0] else 1] = 1
        for gene in sorted(counts):
            symbol = symbols.get(gene, '')
            direct, nearby = counts.get(gene, (0, 0))
            rows.append({
                'site_id': site,
                'location': position,
                'tf_id': gene,
                'symble': '' if symbol.upper() == 'NA' else symbol,
                'motif': direct,
                'flanking': nearby,
            })
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workdir', type=Path, default=DEFAULT_WORKDIR)
    parser.add_argument('--output-dir', type=Path,
                        help='默认使用工作目录下 meth_experiment/CHH_TF_loci')
    args = parser.parse_args()
    try:
        groups = load_groups(args.workdir)
        sites = load_sites(args.workdir / 'meth_experiment/CHH.CV.before.numbered.tsv')
        motifs = load_motifs(args.workdir / 'coex2/all.fimo.out.filt')
        rows = summarize(sites, groups, motifs)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    output_dir = args.output_dir or args.workdir / 'meth_experiment/CHH_TF_loci'
    output_dir.mkdir(parents=True, exist_ok=True)
    output = output_dir / 'CHH_TF_loci.summary.tsv'
    with output.open('w', encoding='utf-8-sig', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    print(f'已生成 {len(rows)} 行汇总表: {output}')
    print(f'纳入正负链结合记录: {len(motifs)} 条；每个位点/TF/类别只计 0 或 1')
    for group, filename in (
        ('共表达', 'CHH_TF_loci.coexpressed.tsv'),
        ('非共表达', 'CHH_TF_loci.non_coexpressed.tsv'),
    ):
        tf_rows = summarize_tf_counts(sites, groups[group], motifs)
        tf_output = output_dir / filename
        with tf_output.open('w', encoding='utf-8-sig', newline='') as handle:
            writer = csv.DictWriter(
                handle, fieldnames=['site_id', 'location', 'tf_id', 'symble', 'motif', 'flanking'],
                delimiter='\t',
            )
            writer.writeheader()
            writer.writerows(tf_rows)
        print(f'已生成{group} TF 表 ({len(tf_rows)} 行): {tf_output}')
    print(f'TF 分组: 共表达 {len(groups["共表达"])}; '
          f'非共表达（all.TF.txt 去除共表达）{len(groups["非共表达"])}')
    for row in rows:
        counts = ', '.join(f'{key}={value}' for key, value in row.items() if key.endswith('_TF数量'))
        print(f'{row["位点ID"]}: {counts}')


if __name__ == '__main__':
    main()
