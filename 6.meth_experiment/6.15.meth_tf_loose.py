#!/usr/bin/env python3
"""Select TF motifs near CHH sites and plot temperature-mean TF/RCA associations.

BED sites are 0-based half-open; FIMO hits are 1-based inclusive.
Expression is ln1p-transformed before aggregation; replicates are not paired.
SPL8 is excluded by default. Spearman uses eight means and all 8! pairings.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import importlib.util
import itertools
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir()) / 'tf-rca-mpl'))
os.environ.setdefault('XDG_CACHE_HOME', str(Path(tempfile.gettempdir()) / 'tf-rca-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.transforms import ScaledTranslation
import numpy as np
import pandas as pd
from scipy.stats import rankdata

PROJECT = Path(__file__).resolve().parents[2]
QA_TOOLS = Path('/datapool/home/2023102768/.codex/skills/nature-figure/scripts')
TEMPERATURES = tuple(range(4, 33, 4))
METRICS = ('total', 'α', 'β', 'α/%', 'β/%')
COLORS = ('#4D4D4D', '#2F6DAE', '#D9822B', '#2F6DAE', '#D9822B')
PREFIX = 'CHH1634_CHH1635_TF_RCA'
PERMUTATIONS = np.array(list(itertools.permutations(range(8))), dtype=np.int16)


def read_tsv(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as h:
        rows = list(csv.reader(h, delimiter='\t'))
    if not rows or not rows[0] or any(not c for c in rows[0]):
        raise ValueError(f'{path}: missing/blank header')
    if len(set(rows[0])) != len(rows[0]):
        raise ValueError(f'{path}: duplicate column names')
    if any(len(r) != len(rows[0]) for r in rows[1:]):
        raise ValueError(f'{path}: inconsistent field count')
    return pd.DataFrame(rows[1:], columns=rows[0])


def check_design(frame, label):
    counts = frame.groupby('temperature').size().to_dict()
    if counts != {t: 3 for t in TEMPERATURES}:
        raise ValueError(f'{label}: expected eight temperatures 4–32 °C, three replicates each; got {counts}')


def numeric(values, label):
    try:
        a = pd.to_numeric(values, errors='raise').astype(float)
    except (ValueError, TypeError) as exc:
        raise ValueError(f'{label}: invalid numeric value') from exc
    if not np.isfinite(a).all() or (a < 0).any():
        raise ValueError(f'{label}: values must be finite and nonnegative')
    return a


def load_inputs(expression_file, rca_file, gene, symbol):
    expr = read_tsv(expression_file)
    if expr.columns[0] != 'Name':
        raise ValueError('expression: first column must be Name')
    target = expr.loc[expr.Name == gene]
    if len(target) != 1:
        raise ValueError(f'expected exactly one {gene} row; found {len(target)}')
    samples = list(expr.columns[1:])
    if any(re.fullmatch(r'\d+_[1-9]\d*', s) is None for s in samples):
        raise ValueError('expression sample IDs must be <temperature>_<positive replicate>')
    spl = pd.DataFrame({'sample': samples, 'temperature': [int(s.split('_')[0]) for s in samples],
                        'raw': numeric(target.iloc[0, 1:], 'TF').to_numpy()})
    if len({tuple(map(int, s.split('_'))) for s in samples}) != len(samples):
        raise ValueError('expression: duplicate normalized sample IDs')
    check_design(spl, 'TF')
    spl['value'] = np.log1p(spl.raw)
    rca = read_tsv(rca_file)
    if any(c not in rca for c in ('Temperature', *METRICS)):
        raise ValueError('RCA: missing required columns')
    rca = rca[['Temperature', *METRICS]].rename(columns={'Temperature': 'temperature'})
    rca['source_row'] = np.arange(2, len(rca) + 2)
    rca['temperature'] = numeric(rca.temperature, 'RCA temperature')
    check_design(rca, 'RCA')
    rca['temperature'] = rca.temperature.astype(int)
    for m in METRICS:
        rca[m] = numeric(rca[m], m)
    if (rca[list(METRICS[3:])] > 100).any().any():
        raise ValueError('RCA: proportions must be within 0–100')
    if not np.allclose(rca.total, rca['α'] + rca['β'], atol=0.02, rtol=1e-6):
        raise ValueError('RCA: total differs from α + β')
    if (abs(rca['α/%'] + rca['β/%'] - 100) > 0.02 + np.spacing(100.)).any():
        raise ValueError('RCA: percentage sum differs from 100')
    return spl.sort_values('sample').reset_index(drop=True), rca.sort_values(['temperature', 'total']).reset_index(drop=True)


def exact_spearman(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if x.shape != (8,) or y.shape != (8,) or not np.isfinite([x, y]).all():
        raise ValueError('Spearman requires eight finite pairs')
    a, b = rankdata(x), rankdata(y)
    a -= a.mean(); b -= b.mean()
    denom = np.linalg.norm(a) * np.linalg.norm(b)
    if denom == 0:
        return np.nan, np.nan, 0, 'constant_variable'
    observed = a @ b
    extreme = int(np.count_nonzero(abs(a[PERMUTATIONS] @ b) >= abs(observed) - 1e-12))
    return float(np.clip(observed / denom, -1, 1)), extreme / 40320, extreme, 'ok'


def bh_adjust(p):
    p = np.asarray(p, float)
    q = np.full_like(p, np.nan)
    finite = np.flatnonzero(np.isfinite(p))
    order = finite[np.argsort(p[finite], kind='stable')]
    if len(order):
        q[order] = np.minimum(1, np.minimum.accumulate((p[order] * len(order) / np.arange(1, len(order)+1))[::-1])[::-1])
    return q


def summarize(spl, rca, gene, symbol):
    groups, points, tests = [], [], []
    for metric in METRICS:
        rows = []
        for t in TEMPERATURES:
            x = spl.loc[spl.temperature == t, 'value']
            block = rca.loc[rca.temperature == t]
            raw = block[metric].to_numpy()
            y = np.log1p(raw) if metric in METRICS[:3] else raw
            row = dict(metric=metric, temperature=t, tf_mean=float(x.mean()), tf_sd=float(x.std(ddof=1)),
                       rca_mean=float(y.mean()), rca_sd=float(y.std(ddof=1)), n_tf=len(x), n_rca=len(y),
                       x_unit='ln(TPM+1)', y_unit='ln(TPM+1)' if metric in METRICS[:3] else '%')
            rows.append(row); groups.append(row)
            for source, rawvalue, value in zip(block.source_row, raw, y):
                points.append(dict(metric=metric, temperature=t, rca_source_row=int(source),
                                   rca_raw=rawvalue, tf_mean=row['tf_mean'], y_plot=value))
        rho, p, extreme, status = exact_spearman([r['tf_mean'] for r in rows], [r['rca_mean'] for r in rows])
        tests.append(dict(gene=gene, symbol=symbol, metric=metric, n=8, rho=rho, p=p,
                          extreme_permutations=extreme, total_permutations=40320, status=status,
                          method='two-sided exhaustive pairing permutation Spearman'))
    stats = pd.DataFrame(tests)
    stats['q_bh'] = bh_adjust(stats.p)
    group_frame, point_frame = pd.DataFrame(groups), pd.DataFrame(points)
    order = group_frame[group_frame.metric == 'total'].sort_values(
        ['tf_mean', 'temperature'], kind='stable')
    positions = {int(t): i for i, t in enumerate(order.temperature)}
    for frame in (group_frame, point_frame):
        frame['x_plot'] = frame.temperature.map(positions).astype(int)
    return group_frame, point_frame, stats


def load_auditor():
    spec = importlib.util.spec_from_file_location('meth_tf_alignment', QA_TOOLS / 'audit_panel_alignment.py')
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def make_figure(groups, points, stats, qa, symbol):
    plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
                         'font.size': 7, 'axes.labelsize': 7, 'axes.titlesize': 8,
                         'xtick.labelsize': 6, 'ytick.labelsize': 6, 'axes.spines.top': False,
                         'axes.spines.right': False, 'axes.linewidth': 0.6,
                         'mathtext.fontset': 'custom', 'mathtext.bf': 'DejaVu Sans:bold:italic',
                         'mathtext.rm': 'DejaVu Sans', 'mathtext.it': 'DejaVu Sans:italic',
                         'mathtext.cal': 'DejaVu Sans',
                         'pdf.fonttype': 42, 'svg.fonttype': 'none', 'savefig.facecolor': 'white'})
    width_mm = 183
    height_mm = 150
    fig = plt.figure(figsize=(width_mm / 25.4, height_mm / 25.4))
    grid = fig.add_gridspec(2, 6, left=0.09, right=0.98, bottom=0.25, top=0.85,
                           wspace=1.15, hspace=0.90)
    audit_rows = []
    panel_labels = []
    for i, (metric, (row, col)) in enumerate(zip(METRICS, [(0,0), (0,2), (0,4), (1,1), (1,3)])):
        ax = fig.add_subplot(grid[row, col:col+2])
        ax.set_gid(chr(97+i))
        sg = groups[groups.metric == metric].sort_values('x_plot')
        sp = points[points.metric == metric]
        st = stats.iloc[i]
        ax.scatter(sp.x_plot, sp.y_plot, s=9, color=COLORS[i], alpha=.4, edgecolors='none', zorder=3)
        ax.plot(sg.x_plot, sg.rca_mean, color=COLORS[i], lw=.7, alpha=.8, zorder=2)
        ax.errorbar(sg.x_plot, sg.rca_mean, yerr=sg.rca_sd,
                    fmt='o', ms=3, color=COLORS[i], elinewidth=.55, capsize=2,
                    capthick=.55, markeredgecolor='white', markeredgewidth=.35, zorder=4)
        bottom = min(sp.y_plot.min(), (sg.rca_mean-sg.rca_sd).min())
        top = max(sp.y_plot.max(), (sg.rca_mean+sg.rca_sd).max())
        margin = max(top-bottom, .1) * .12
        ax.set_ylim(bottom-margin, top+margin)
        ax.set_xlim(-.5, 7.5)
        ax.set_xticks(sg.x_plot, labels=[f'{value:.3f}' for value in sg.tf_mean])
        for label in ax.get_xticklabels():
            label.set_rotation(60)
            label.set_rotation_mode('anchor')
            label.set_ha('right')
            label.set_va('top')
        ax.set_xlabel(f'Mean {symbol} expression, ln(TPM+1)', labelpad=8)
        ax.set_ylabel('ln(TPM+1)' if i<3 else 'Proportion (%)', labelpad=3)
        ax.set_title((r'Total $\mathbf{RCA}$', r'$\mathbf{RCA}$ α', r'$\mathbf{RCA}$ β',
                      r'$\mathbf{RCA}$ α proportion', r'$\mathbf{RCA}$ β proportion')[i],
                     pad=25, fontweight='bold')
        panel_labels.append(ax.text(0, 1, f'({chr(97+i)})', fontweight='bold', fontsize=8,
                ha='left', va='bottom',
                transform=ax.transAxes + ScaledTranslation(-20/72, 24/72, fig.dpi_scale_trans)))
        stat = f'ρ = {st.rho:.3f}; P = {st.p:.3f}' if st.status == 'ok' else 'ρ = NA; P = NA'
        ax.text(.5, 1, stat, ha='center', va='bottom', fontsize=6,
                transform=ax.transAxes + ScaledTranslation(0, 8/72, fig.dpi_scale_trans))
        ax.tick_params(length=2.5, width=.5)
        audit_rows.append(dict(panel=chr(97+i), metric=metric, rca_replicates=len(sp), means=len(sg),
                               horizontal_sd_bars=0, vertical_sd_bars=len(sg), distinct_x_ticks=len(ax.get_xticks()),
                               temperature_order=','.join(map(str, sg.temperature))))
    captions = [
        r'Faint points: all $\mathbf{RCA}$ replicates; solid points: mean ± sample SD (n = 3 per temperature). Lines are descriptive only.',
        f'X: eight equally spaced temperature groups ordered by mean {symbol} expression; ticks show measured means.',
        'Spearman: n = 8 temperature-mean pairs; exact two-sided permutation P. BH q values are in the statistics table.',
        'Expression: ln(TPM+1) before aggregation; proportions: %. Cross-temperature associations do not imply regulation.'
    ]
    for y, caption in zip((.120, .093, .066, .039), captions):
        fig.text(.5, y, caption, ha='center', fontsize=6)
    fig.canvas.draw()
    auditor = load_auditor()
    auditor.require_matplotlib_panel_alignment(fig, tolerance_pt=1.5,
        gutter_tolerance_pt=1.5, strict=True)
    # The automatic label detector only recognizes bare letters. Measure the
    # displayed parenthesized labels explicitly for the same alignment gate.
    manifest = auditor.matplotlib_layout_manifest(fig)
    for panel, label in zip(manifest['panels'], panel_labels):
        xy = label.get_transform().transform(label.get_position())
        anchor = fig.dpi_scale_trans.inverted().transform(xy) * 72
        panel['panel_label'] = label.get_text()
        panel['panel_label_anchor_pt'] = list(map(float, anchor))
    report = auditor.audit_layout_manifest(manifest,
        tolerance_pt=1.5, gutter_tolerance_pt=1.5, require_panel_labels=True)
    auditor.write_json_report(report, qa / 'alignment.json')
    if auditor.exit_code(report, strict=True):
        raise RuntimeError(auditor.render_text(report, strict=True))
    return fig, audit_rows


def run_audit(name, command, qa):
    result = subprocess.run([sys.executable, *map(str, command)], capture_output=True, text=True)
    (qa / name).write_text(result.stdout + result.stderr, encoding='utf-8')
    if result.returncode:
        raise RuntimeError(f'{name} failed: {result.stdout[-2000:]} {result.stderr[-500:]}')


def site_distance(chrom, start, end, motif_chrom, motif_start, motif_end):
    """Nearest nucleotide distance: BED site versus 1-based closed motif.

    A single-base BED [s,s+1) denotes the 1-based nucleotide s+1.
    Adjacent nucleotides have distance one, overlapping nucleotides zero.
    """
    if chrom != motif_chrom:
        return None
    return max(motif_start - end, start + 1 - motif_end, 0)


def select_candidates(candidate_file, fimo_file, sites_file, site_ids, window, excluded):
    candidates = read_tsv(candidate_file)
    required = {'motif_id', 'symbol', 'set_type', 'region'}
    if not required.issubset(candidates.columns) or candidates.motif_id.duplicated().any():
        raise ValueError('candidate table: missing fields or duplicate motif_id')
    if window < 0 or not site_ids or len(set(site_ids)) != len(site_ids):
        raise ValueError('window must be nonnegative and site IDs nonempty and unique')
    positions = {}
    with Path(sites_file).open() as handle:
        for line, row in enumerate(csv.reader(handle, delimiter='\t'), 1):
            if len(row) != 6:
                raise ValueError(f'sites line {line}: expected six fields')
            identifier, chrom, start, end = row[:4]
            start, end = int(start), int(end)
            if start < 0 or end <= start or identifier in positions:
                raise ValueError(f'sites: invalid interval or duplicate ID {identifier}')
            positions[identifier] = (chrom, start, end)
    if set(site_ids) - positions.keys():
        raise ValueError(f'missing site IDs: {sorted(set(site_ids) - positions.keys())}')
    motifs = []
    with Path(fimo_file).open() as handle:
        for line, row in enumerate(csv.reader(handle, delimiter='\t'), 1):
            if not row or row[0].startswith('#'):
                continue
            if len(row) != 10:
                raise ValueError(f'FIMO line {line}: expected ten fields')
            gene, alt, chrom, start, end, strand, score, p, q, sequence = row
            start, end = int(start), int(end)
            if start < 1 or end < start or strand not in ('+', '-'):
                raise ValueError(f'FIMO line {line}: invalid interval/strand')
            values = np.array([score, p, q], dtype=float)
            if not np.isfinite(values).all() or not (0 <= values[1] <= 1 and 0 <= values[2] <= 1):
                raise ValueError(f'FIMO line {line}: invalid score/probability')
            motifs.append(dict(fimo_line=line, motif_id=gene, motif_alt_id=alt,
                               chromosome=chrom, motif_start_1based=start, motif_end_1based=end,
                               motif_start_0based=start-1, motif_end_exclusive=end,
                               strand=strand, score=values[0], fimo_p=values[1], fimo_q=values[2],
                               matched_sequence=sequence))
    evidence, selection = [], []
    for candidate in candidates.to_dict('records'):
        gene = candidate['motif_id']
        records = [m for m in motifs if m['motif_id'] == gene]
        matched_lines, matched_sites, distances = set(), set(), []
        for motif in records:
            for site in site_ids:
                chrom, start, end = positions[site]
                distance = site_distance(chrom, start, end, motif['chromosome'],
                                         motif['motif_start_1based'], motif['motif_end_1based'])
                hit = distance is not None and distance <= window
                if distance is not None:
                    distances.append(distance)
                if hit:
                    matched_lines.add(motif['fimo_line'])
                    matched_sites.add(site)
                evidence.append(dict(**motif, symbol=candidate['symbol'], site=site,
                                     site_chromosome=chrom, site_start_0based=start,
                                     site_end_exclusive=end, distance_bp=distance, within_window=hit))
        spatial = bool(matched_lines)
        analyze = spatial and gene not in excluded
        reason = ('selected' if analyze else 'already_plotted_excluded' if spatial
                  else 'no_fimo_motif' if not records else 'outside_window_or_chromosome')
        selection.append(dict(**candidate, spatial_match=spatial, analyze=analyze,
                              reason=reason, exclusion_note='已绘制，本次排除' if spatial and gene in excluded else '',
                              fimo_records=len(records), matched_fimo_records=len(matched_lines),
                              matched_sites=','.join(s for s in site_ids if s in matched_sites),
                              minimum_distance_bp=min(distances) if distances else None))
    return pd.DataFrame(selection), pd.DataFrame(evidence)


def write_table(table, path):
    table.to_csv(path, sep='\t', index=False, na_rep='NA', float_format='%.12g')


def main(argv=None):
    base = PROJECT / 'list'
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate-file', type=Path, default=base / 'RCA/TF_CHH/AT2G39730_CHH_combined_historical.panel-b-points.tsv')
    parser.add_argument('--fimo-file', type=Path, default=base / 'coex2/all.fimo.out.filt')
    parser.add_argument('--sites-file', type=Path, default=base / 'meth_experiment/CHH.CV.before.numbered.tsv')
    parser.add_argument('--expression-file', type=Path, default=base / 'coex_experiment/RCA.tran.txt')
    parser.add_argument('--rca-file', type=Path, default=base / 'RCA/RCA.rnadata.tsv')
    parser.add_argument('--output-dir', type=Path, default=base / 'meth_experiment')
    parser.add_argument('--sites', nargs='+', default=['CHH_1634', 'CHH_1635'])
    parser.add_argument('--window', type=int, default=100)
    parser.add_argument('--exclude', nargs='*', default=['AT1G02065'])
    parser.add_argument('--validate-only', action='store_true')
    args = parser.parse_args(argv)
    selection, evidence = select_candidates(args.candidate_file, args.fimo_file, args.sites_file,
                                            args.sites, args.window, set(args.exclude))
    chosen = selection.loc[selection.analyze]
    if chosen.empty:
        raise ValueError('no TFs remain after spatial selection and exclusions')
    all_groups, all_points, all_stats, tf_sources = [], [], [], []
    for row in chosen.itertuples():
        tf, rca = load_inputs(args.expression_file, args.rca_file, row.motif_id, row.symbol)
        groups, points, stats = summarize(tf, rca, row.motif_id, row.symbol)
        for frame in (tf, groups, points):
            frame.insert(0, 'symbol', row.symbol)
            frame.insert(0, 'gene', row.motif_id)
        all_groups.append(groups)
        all_points.append(points)
        all_stats.append(stats)
        tf_sources.append(tf)
    groups, points, stats, sources = [pd.concat(frames, ignore_index=True) for frames in
                                     (all_groups, all_points, all_stats, tf_sources)]
    stats['q_bh'] = bh_adjust(stats.p)
    summary = dict(candidate_tfs=len(selection), spatially_matched_tfs=int(selection.spatial_match.sum()),
                   analyzed_tfs=len(chosen), tests=len(stats), sites=args.sites, window_bp=args.window,
                   excluded=args.exclude, tf_order=list(chosen.symbol),
                   motif_hits=int(evidence.loc[evidence.within_window, 'fimo_line'].nunique()),
                   gene_site_motif_hits=int(evidence.within_window.sum()))
    if args.validate_only:
        print(json.dumps(summary, ensure_ascii=False, indent=2))
        return 0
    if importlib.util.find_spec('fitz') is None:
        raise RuntimeError('PyMuPDF is required for PDF collision QA; NOT AUDITABLE')
    args.output_dir.mkdir(parents=True, exist_ok=True)
    qa_root = args.output_dir / f'{PREFIX}_QA'
    qa_root.mkdir(exist_ok=True)
    run_audit('source-preflight.txt', [QA_TOOLS / 'validate_figure.py', __file__], qa_root)
    rca_sources = rca.melt(id_vars=['temperature', 'source_row'], value_vars=list(METRICS),
                           var_name='metric', value_name='raw')
    rca_sources['value'] = rca_sources.raw
    expression_mask = rca_sources.metric.isin(METRICS[:3])
    rca_sources.loc[expression_mask, 'value'] = np.log1p(rca_sources.loc[expression_mask, 'raw'])
    for name, table in [('selection', selection), ('site-motif-evidence', evidence),
                        ('site-motif-hits', evidence.loc[evidence.within_window]),
                        ('tf-source', sources), ('rca-source', rca_sources),
                        ('groups', groups), ('points', points), ('statistics', stats)]:
        write_table(table, args.output_dir / f'{PREFIX}.{name}.tsv')
    for row in chosen.itertuples():
        stem = args.output_dir / f'{PREFIX}_{row.motif_id}_{row.symbol}'
        qa = qa_root / row.motif_id
        qa.mkdir(exist_ok=True)
        sg, sp, st = [frame.loc[frame.gene == row.motif_id].reset_index(drop=True)
                      for frame in (groups, points, stats)]
        fig, panels = make_figure(sg, sp, st, qa, row.symbol)
        try:
            fig.savefig(stem.with_suffix('.pdf'))
            run_audit('pdf-text.json', [QA_TOOLS / 'audit_pdf_text.py', stem.with_suffix('.pdf'), '--min-pt', '5', '--json'], qa)
            run_audit('collision-audit.txt', [QA_TOOLS / 'audit_figure_collisions.py', stem.with_suffix('.pdf'),
                      '--strict', '--json-out', qa / 'collision-audit.json'], qa)
            fig.savefig(stem.with_suffix('.svg'))
            fig.savefig(stem.with_suffix('.png'), dpi=600)
            fig.savefig(stem.with_suffix('.tiff'), dpi=600, pil_kwargs={'compression': 'tiff_lzw'})
        finally:
            plt.close(fig)
        write_table(pd.DataFrame(panels), qa / 'panels.tsv')
        print(f'Exported {row.symbol} ({row.motif_id}): five panels, automated QA PASS', flush=True)
    summary.update(status='automated_QA_PASS', figure_size_mm=[183, 150], raster_dpi=600,
                   question='Are TFs with motifs near the selected CHH sites associated with RCA across temperatures?',
                   archetype='quantitative grid; total expression, isoforms and complementary proportions',
                   coordinates='BED 0-based half-open sites; FIMO 1-based closed; nearest nucleotide distance <= window',
                   transforms='ln1p each expression replicate before mean and sample SD (ddof=1); percentages unchanged',
                   sample_pairing='temperature only; no individual pairing inferred',
                   spearman='8 means; average tied ranks; exhaustive 40320 two-sided pairing permutations',
                   multiplicity=f'BH across all {len(stats)} tests; undefined tests retained as NA and omitted from BH',
                   x_axis='equally spaced groups ordered by measured TF means, ties by temperature',
                   uncertainty='RCA sample SD drawn; TF sample SD retained in groups table',
                   interpretation='descriptive across-temperature association; not evidence of regulation; proportions complementary',
                   input_sha256={name: hashlib.sha256(getattr(args, name).read_bytes()).hexdigest()
                                 for name in ('candidate_file', 'fimo_file', 'sites_file', 'expression_file', 'rca_file')},
                   source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   versions={'python': sys.version.split()[0], 'numpy': np.__version__,
                             'pandas': pd.__version__, 'matplotlib': matplotlib.__version__})
    (qa_root / 'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    (qa_root / 'README.md').write_text(
        '# TF–RCA temperature associations\n\n'
        'SPL8 is excluded by default because it has already been plotted. Selection evidence retains it.\n'
        'Panels a–c show total RCA, alpha and beta abundance; d–e show complementary percentages.\n'
        'Each panel retains all 24 RCA observations; filled points show eight means with sample SD (n=3).\n'
        'Replicate biological/technical status is not encoded in the source tables.\n'
        'TF and RCA are joined only by temperature. Expression is ln(TPM+1) transformed before averaging.\n'
        'X positions are categorical, ordered by measured TF mean; tick labels show those means.\n'
        'Equal means remain separate groups, ordered by temperature. Lines are descriptive, not fitted.\n'
        'Spearman uses measured means, not categorical positions, with exact two-sided permutation P.\n'
        'BH q values cover all analyzed TF–metric tests and are reported in the statistics table.\n'
        'No imputation or expression-based exclusions. Constant-variable tests are reported as NA.\n'
        'These across-temperature associations do not establish regulation or methylation-mediated effects.\n'
        'The summary records input hashes, coordinate rules, exclusions, and numerical methods.\n', encoding='utf-8')
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (ValueError, OSError, RuntimeError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        raise SystemExit(1)
