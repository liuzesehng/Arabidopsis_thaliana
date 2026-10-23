#!/usr/bin/env python3
"""SPL8 versus RCA: temperature-mean Spearman associations with all RCA replicates.

No individual pairing is inferred. Expression is ln1p-transformed before
aggregation. Eight equally spaced temperature groups are ordered by SPL8 mean;
vertical error bars are RCA sample SD and connecting lines are descriptive only.
Spearman uses measured group means, never categorical positions, with all 8! pairings.
An additional two-panel figure fits quadratic models to the two RCA proportions.
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

os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir()) / 'spl8-rca-mpl'))
os.environ.setdefault('XDG_CACHE_HOME', str(Path(tempfile.gettempdir()) / 'spl8-rca-cache'))
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
PREFIX = 'SPL8_RCA_correlation'
GENE = 'AT1G02065'
PERMUTATIONS = np.array(list(itertools.permutations(range(8))), dtype=np.int16)


# Optional external QA: keep each plotting script usable on its own.
_OPTIONAL_QA_SKIPPED = []


def _skip_optional_qa(reason, output=None):
    """Record an unavailable check without claiming that it passed."""
    if reason not in _OPTIONAL_QA_SKIPPED:
        _OPTIONAL_QA_SKIPPED.append(reason)
        print(f"[QA SKIPPED] {reason}", file=sys.stderr)
    report = {"verdict": "SKIPPED", "status": "SKIPPED", "auditable": False,
              "reason": reason}
    if output is not None:
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return report


def _skip_missing_qa_command(command, output=None):
    """Only missing optional scripts are skipped; executed checks stay strict."""
    script = Path(command[1])
    if script.is_file():
        return False
    reason = f"Optional QA tool not found: {script}"
    _skip_optional_qa(reason, output)
    if "--json-out" in command:
        _skip_optional_qa(reason, command[command.index("--json-out") + 1])
    return True


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


def load_inputs(expression_file, rca_file):
    expr = read_tsv(expression_file)
    if expr.columns[0] != 'Name':
        raise ValueError('expression: first column must be Name')
    target = expr.loc[expr.Name == GENE]
    if len(target) != 1:
        raise ValueError(f'expected exactly one {GENE} row; found {len(target)}')
    samples = list(expr.columns[1:])
    if any(re.fullmatch(r'\d+_[1-9]\d*', s) is None for s in samples):
        raise ValueError('expression sample IDs must be <temperature>_<positive replicate>')
    spl = pd.DataFrame({'sample': samples, 'temperature': [int(s.split('_')[0]) for s in samples],
                        'raw': numeric(target.iloc[0, 1:], 'SPL8').to_numpy()})
    if len({tuple(map(int, s.split('_'))) for s in samples}) != len(samples):
        raise ValueError('expression: duplicate normalized sample IDs')
    check_design(spl, 'SPL8')
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


def summarize(spl, rca):
    groups, points, tests = [], [], []
    for metric in METRICS:
        rows = []
        for t in TEMPERATURES:
            x = spl.loc[spl.temperature == t, 'value']
            block = rca.loc[rca.temperature == t]
            raw = block[metric].to_numpy()
            y = np.log1p(raw) if metric in METRICS[:3] else raw
            row = dict(metric=metric, temperature=t, spl8_mean=float(x.mean()), spl8_sd=float(x.std(ddof=1)),
                       rca_mean=float(y.mean()), rca_sd=float(y.std(ddof=1)), n_spl8=len(x), n_rca=len(y),
                       x_unit='ln(TPM+1)', y_unit='ln(TPM+1)' if metric in METRICS[:3] else '%')
            rows.append(row); groups.append(row)
            for source, rawvalue, value in zip(block.source_row, raw, y):
                points.append(dict(metric=metric, temperature=t, rca_source_row=int(source),
                                   rca_raw=rawvalue, spl8_mean=row['spl8_mean'], y_plot=value))
        rho, p, extreme, status = exact_spearman([r['spl8_mean'] for r in rows], [r['rca_mean'] for r in rows])
        tests.append(dict(gene=GENE, symbol='SPL8', metric=metric, n=8, rho=rho, p=p,
                          extreme_permutations=extreme, total_permutations=40320, status=status,
                          method='two-sided exhaustive pairing permutation Spearman'))
    stats = pd.DataFrame(tests)
    stats['q_bh'] = bh_adjust(stats.p)
    group_frame, point_frame = pd.DataFrame(groups), pd.DataFrame(points)
    order = group_frame[group_frame.metric == 'total'].sort_values(
        ['spl8_mean', 'temperature'], kind='stable')
    positions = {int(t): i for i, t in enumerate(order.temperature)}
    for frame in (group_frame, point_frame):
        frame['x_plot'] = frame.temperature.map(positions).astype(int)
    return group_frame, point_frame, stats


def load_auditor():
    if not Path(QA_TOOLS / 'audit_panel_alignment.py').is_file():
        _skip_optional_qa(f"Optional QA tool not found: {QA_TOOLS / 'audit_panel_alignment.py'}")
        return None
    spec = importlib.util.spec_from_file_location('spl8_alignment', QA_TOOLS / 'audit_panel_alignment.py')
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    try:
        spec.loader.exec_module(module)
    except ImportError as exc:
        _skip_optional_qa(f"Optional alignment dependency unavailable: {exc}")
        return None
    return module


def make_figure(groups, points, stats, qa):
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
        ax.set_xticks(sg.x_plot, labels=[f'{value:.3f}' for value in sg.spl8_mean])
        for label in ax.get_xticklabels():
            label.set_rotation(60)
            label.set_rotation_mode('anchor')
            label.set_ha('right')
            label.set_va('top')
        ax.set_xlabel('Mean SPL8 expression, ln(TPM+1)', labelpad=8)
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
        'X: eight equally spaced temperature groups ordered by mean SPL8 expression; ticks show measured means.',
        'Spearman: n = 8 temperature-mean pairs; exact two-sided permutation P. BH q values are in the statistics table.',
        'Expression: ln(TPM+1) before aggregation; proportions: %. Cross-temperature associations do not imply regulation.'
    ]
    for y, caption in zip((.120, .093, .066, .039), captions):
        fig.text(.5, y, caption, ha='center', fontsize=6)
    fig.canvas.draw()
    auditor = load_auditor()
    if auditor is None:
        _skip_optional_qa("audit_panel_alignment.py is unavailable", qa / 'alignment.json')
    else:
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
    if _skip_missing_qa_command([sys.executable, *map(str, command)], qa / name):
        return
    result = subprocess.run([sys.executable, *map(str, command)], capture_output=True, text=True)
    (qa / name).write_text(result.stdout + result.stderr, encoding='utf-8')
    if result.returncode:
        raise RuntimeError(f'{name} failed: {result.stdout[-2000:]} {result.stderr[-500:]}')


def fit_quadratic(x, y):
    """Unweighted OLS on eight group means; conditional, pointwise mean CI.

    Center/scale x for stability. Test the added quadratic term against a
    linear model with a nested F(1, 5) test, separately from overall F(2, 5).
    """
    from scipy.stats import f as f_distribution, t as t_distribution
    x, y = np.asarray(x, float), np.asarray(y, float)
    if x.shape != (8,) or y.shape != (8,) or not np.isfinite([x, y]).all():
        raise ValueError('quadratic fit requires eight finite group-mean pairs')
    center, scale = float(x.mean()), float(x.std(ddof=0))
    if scale == 0 or len(np.unique(x)) < 3:
        raise ValueError('quadratic fit requires at least three distinct SPL8 means')
    z = (x-center) / scale
    design = np.column_stack([np.ones(8), z, z*z])
    coefficient, _, rank, _ = np.linalg.lstsq(design, y, rcond=None)
    if rank != 3:
        raise ValueError('quadratic fit design is rank deficient')
    fitted = design @ coefficient
    residual = y-fitted
    sse, sst = float(residual @ residual), float(np.sum((y-y.mean())**2))
    if sst <= np.finfo(float).eps * max(float(y@y), 1):
        raise ValueError('quadratic fit cannot infer curvature for a constant response')
    mse = sse / 5
    linear = design[:, :2] @ np.linalg.lstsq(design[:, :2], y, rcond=None)[0]
    linear_sse = float(np.sum((y-linear)**2))
    p_curvature = float(f_distribution.sf(max(0., linear_sse-sse)/mse, 1, 5)) if mse else 0.
    p_overall = float(f_distribution.sf(max(0., sst-sse)/2/mse, 2, 5)) if mse else 0.
    covariance = np.linalg.inv(design.T @ design) * mse
    critical = float(t_distribution.ppf(.975, 5))
    # Convert back to y = intercept + linear*x + quadratic*x² for source data.
    transform = np.array([[1, -center/scale, center**2/scale**2],
                          [0, 1/scale, -2*center/scale**2], [0, 0, 1/scale**2]])
    raw_coef = transform @ coefficient
    raw_se = np.sqrt(np.maximum(0., np.diag(transform @ covariance @ transform.T)))
    curved = abs(coefficient[2]) > 1e-12 * max(1., float(np.max(abs(coefficient))))
    vertex = float(center-scale*coefficient[1]/(2*coefficient[2])) if curved else np.nan
    inside = bool(curved and x.min() < vertex < x.max())
    shape = ('U-shaped' if coefficient[2] > 0 else 'inverted-U-shaped') if inside else 'no interior turning point'
    grid = np.linspace(x.min(), x.max(), 301)
    zg = (grid-center)/scale
    grid_design = np.column_stack([np.ones(len(grid)), zg, zg*zg])
    prediction = grid_design @ coefficient
    se_mean = np.sqrt(np.maximum(0., np.einsum('ij,jk,ik->i', grid_design, covariance, grid_design)))
    summary = dict(n=8, residual_df=5, model='unweighted quadratic OLS on temperature means',
                   intercept=raw_coef[0], linear=raw_coef[1], quadratic=raw_coef[2],
                   quadratic_ci_low=raw_coef[2]-critical*raw_se[2],
                   quadratic_ci_high=raw_coef[2]+critical*raw_se[2],
                   r_squared=1-sse/sst, adjusted_r_squared=1-(sse/5)/(sst/7),
                   rmse=np.sqrt(mse), p_overall=p_overall, p_curvature=p_curvature,
                   vertex_x=vertex, vertex_inside_observed_range=inside, fitted_shape=shape,
                   xmin=float(x.min()), xmax=float(x.max()), x_center=center, x_scale=scale)
    curve = pd.DataFrame(dict(x=grid, fitted=prediction,
                              ci_low=prediction-critical*se_mean, ci_high=prediction+critical*se_mean))
    return summary, curve, fitted, residual


def format_quadratic_ticks(ax):
    """Eight equal-spaced labels, all at 45 degrees on one baseline."""
    for label in ax.get_xticklabels():
        label.set_rotation(45)
        label.set_rotation_mode('anchor')
        label.set_ha('right')
        label.set_va('top')
    ax.figure.canvas.draw()
    renderer = ax.figure.canvas.get_renderer()
    lowest = min(label.get_window_extent(renderer).y0 for label in ax.get_xticklabels())
    xlabel_y = ax.transAxes.inverted().transform((ax.bbox.x0, lowest-9*ax.figure.dpi/72))[1]
    ax.xaxis.set_label_coords(.5,xlabel_y)


def export_quadratic(groups, points, output_dir):
    """Additional two-panel proportion figure using the first figure's equal-spaced group axis."""
    prefix = 'SPL8_RCA_proportion_quadratic'
    qa = output_dir / f'{prefix}_QA'
    qa.mkdir(exist_ok=True)
    fits, curves, observations = [], [], []
    for metric in METRICS[3:]:
        group = groups[groups.metric == metric].sort_values('spl8_mean').copy()
        summary, curve, fitted, residual = fit_quadratic(group.spl8_mean, group.rca_mean)
        fits.append(dict(metric=metric, **summary))
        # Fit on measured x, then map the curve to the same categorical
        # display coordinates as the first figure; do not refit on ranks.
        measured_x = group.spl8_mean.to_numpy()
        if not np.all(np.diff(measured_x) > 0):
            raise ValueError('equal-spaced quadratic display requires distinct group means')
        curve['x_plot'] = np.interp(curve.x, measured_x, group.x_plot)
        curve.insert(0, 'metric', metric)
        curves.append(curve)
        group['fitted'] = fitted
        group['residual'] = residual
        observations.append(group)
    stats = pd.DataFrame(fits)
    stats['q_curvature_bh'] = bh_adjust(stats.p_curvature)
    stats.to_csv(output_dir / f'{prefix}.statistics.tsv', sep='\t', index=False, na_rep='NA', float_format='%.12g')
    pd.concat(curves, ignore_index=True).to_csv(output_dir / f'{prefix}.curves.tsv', sep='\t', index=False, float_format='%.12g')
    pd.concat(observations, ignore_index=True).to_csv(output_dir / f'{prefix}.groups.tsv', sep='\t', index=False, float_format='%.12g')
    raw_points = points[points.metric.isin(METRICS[3:])].copy()
    raw_points.to_csv(output_dir / f'{prefix}.points.tsv', sep='\t', index=False, float_format='%.12g')
    width_mm = 183
    height_mm = 125
    fig = plt.figure(figsize=(width_mm / 25.4, height_mm / 25.4))
    grid = fig.add_gridspec(1, 2, left=.09, right=.98, bottom=.44, top=.82, wspace=.30)
    labels = []
    for i, metric in enumerate(METRICS[3:]):
        ax = fig.add_subplot(grid[0,i]); ax.set_gid(chr(97+i))
        group, curve, fit = observations[i], curves[i], fits[i]
        raw = raw_points[raw_points.metric == metric]
        color = COLORS[i+3]
        ax.fill_between(curve.x_plot.to_numpy(), curve.ci_low.to_numpy(), curve.ci_high.to_numpy(),
                        color=color, alpha=.15, linewidth=0, zorder=1)
        ax.plot(curve.x_plot, curve.fitted, color=color, lw=1.1, zorder=2)
        ax.scatter(raw.x_plot, raw.y_plot, s=10, color=color, alpha=.35, edgecolors='none', zorder=3)
        ax.errorbar(group.x_plot, group.rca_mean, yerr=group.rca_sd,
                    fmt='o', ms=3.8, color=color, elinewidth=.65, capsize=2,
                    capthick=.65, markeredgecolor='white', markeredgewidth=.4, zorder=4)
        low = min(float(curve.ci_low.min()), float((group.rca_mean-group.rca_sd).min()), float(raw.y_plot.min()))
        high = max(float(curve.ci_high.max()), float((group.rca_mean+group.rca_sd).max()), float(raw.y_plot.max()))
        margin = max(high-low, .1)*.08
        ax.set_ylim(low-margin, high+margin)
        ax.set_xlim(-.5, 7.5)
        ax.set_xticks(group.x_plot, labels=[
            f'{row.spl8_mean:.3f}' for row in group.itertuples()])
        ax.set_xlabel('Mean SPL8 expression, ln(TPM+1)', labelpad=5)
        ax.set_ylabel('Proportion (%)', labelpad=4)
        ax.set_title((r'$\mathbf{RCA}$ α proportion', r'$\mathbf{RCA}$ β proportion')[i], pad=26, fontweight='bold')
        labels.append(ax.text(0, 1, f'({chr(97+i)})', fontweight='bold', fontsize=8, ha='left', va='bottom',
                      transform=ax.transAxes+ScaledTranslation(-20/72,25/72,fig.dpi_scale_trans)))
        pv=fit['p_curvature']; ptext=f'{pv:.3f}' if pv>=.001 else f'{pv:.2e}'
        ax.text(.5,1,f'R² = {fit["r_squared"]:.3f}; P (curvature) = {ptext}', ha='center', va='bottom',
                fontsize=6.5, transform=ax.transAxes+ScaledTranslation(0,8/72,fig.dpi_scale_trans))
        ax.tick_params(length=2.5,width=.5)
        format_quadratic_ticks(ax)
    captions = [
        'X: eight equally spaced groups, as in the first figure; ticks show mean ln(TPM+1). Fits use n = 8 measured pairs.',
        r'Faint points: all 3 $\mathbf{RCA}$ replicates per temperature at the SPL8 group mean; solid points: $\mathbf{RCA}$ mean ± sample SD.',
        'Shading: pointwise 95% confidence interval for the fitted mean, conditional on x (OLS residual df = 5).',
        'P (curvature): added quadratic term versus a linear model. α and β proportions are complementary.',
        'Curves and bands are mapped from measured x to group positions for display; fits are exploratory, not causal.'
    ]
    for y,text in zip((.210,.176,.142,.108,.074),captions):
        fig.text(.5,y,text,ha='center',fontsize=6)
    fig.canvas.draw()
    auditor = load_auditor()
    if auditor is None:
        _skip_optional_qa("audit_panel_alignment.py is unavailable", qa / 'alignment.json')
    else:
        auditor.require_matplotlib_panel_alignment(fig,tolerance_pt=1.5,gutter_tolerance_pt=1.5,strict=True)
        manifest=auditor.matplotlib_layout_manifest(fig)
        for panel,label in zip(manifest['panels'],labels):
            xy=label.get_transform().transform(label.get_position())
            panel['panel_label']=label.get_text()
            panel['panel_label_anchor_pt']=list(map(float,fig.dpi_scale_trans.inverted().transform(xy)*72))
        report=auditor.audit_layout_manifest(manifest,tolerance_pt=1.5,gutter_tolerance_pt=1.5,require_panel_labels=True)
        auditor.write_json_report(report,qa/'alignment.json')
        if auditor.exit_code(report,strict=True):
            raise RuntimeError(auditor.render_text(report,strict=True))
    stem=output_dir/prefix
    try:
        fig.savefig(stem.with_suffix('.pdf'))
        run_audit('pdf-text.json',[QA_TOOLS/'audit_pdf_text.py',stem.with_suffix('.pdf'),'--min-pt','5','--json'],qa)
        run_audit('collision-audit.txt',[QA_TOOLS/'audit_figure_collisions.py',stem.with_suffix('.pdf'),
                  '--strict','--json-out',qa/'collision-audit.json'],qa)
        fig.savefig(stem.with_suffix('.svg'))
        fig.savefig(stem.with_suffix('.png'),dpi=600)
        fig.savefig(stem.with_suffix('.tiff'),dpi=600,pil_kwargs={'compression':'tiff_lzw'})
    finally:
        plt.close(fig)
    (qa/'README.md').write_text(
        '# Quadratic SPL8–RCA proportion fits\n\n'
        'Question: do the eight temperature means show curvature in RCA isoform proportions versus SPL8 expression?\n'
        'All eight means receive equal OLS weight; the 24 replicate points are displayed but are not 24 fitted pairs.\n'
        'The model is y = intercept + linear*x + quadratic*x², fit after centering/scaling measured x for stability.\n'
        'x is the mean of replicate ln(TPM+1); y is the mean of untransformed percentages. No categorical positions enter the fit.\n'
        'All eight temperature groups have ticks labelled only with measured SPL8 means; labels rotate 45 degrees on a common baseline at equal-spaced group positions.\n'
        'The second figure uses the first figure\'s eight equal-spaced group positions and identical tick values.\n'
        'Curves and CI are fitted on measured expression, then mapped piecewise linearly to these display positions.\n'
        'Bands are pointwise 95% Student-t mean confidence intervals with residual df=5, not prediction intervals or SD.\n'
        'OLS inference assumes independent, normal, equal-variance residuals and treats measured x as fixed; it does not propagate SPL8 measurement error.\n'
        'The nested F(1,5) curvature test compares quadratic versus linear; overall F(2,5) and BH across the two curvature tests are in the table.\n'
        'The proportions are complementary and do not constitute independent validation. This small-sample, exploratory model was chosen after visual inspection.\n'
        'Turning points are reported on the transformed expression scale; no curve is extrapolated beyond observed x.\n'
        'Positive/negative curvature with an interior vertex describes a fitted U/inverted U; it does not establish a regulatory mechanism.\n',encoding='utf-8')
    summary=dict(status=("automated_QA_INCOMPLETE" if _OPTIONAL_QA_SKIPPED else "automated_QA_PASS"),figure_size_mm=[183,125],raster_dpi=600,
                 panels=2,temperature_means_per_panel=8,rca_replicates_per_panel=24,
                 x_ticks_per_panel=8, tick_labels='SPL8 mean only; 45-degree labels at the same eight categorical positions as the first figure',
                 model='quadratic OLS; measured SPL8 means; equal weight per temperature',
                 ci='pointwise 95% fitted-mean t interval, df=5; conditional on measured x',
                 source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (qa/'summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    print(stats[['metric','r_squared','p_curvature','vertex_x','fitted_shape']].to_string(index=False))
    return stats


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--expression-file', type=Path, default=PROJECT / 'list/coex_experiment/RCA.tran.txt')
    parser.add_argument('--rca-file', type=Path, default=PROJECT / 'list/RCA/RCA.rnadata.tsv')
    parser.add_argument('--output-dir', type=Path, default=PROJECT / 'list/meth_experiment')
    parser.add_argument('--validate-only', action='store_true')
    args = parser.parse_args(argv)
    spl, rca = load_inputs(args.expression_file, args.rca_file)
    if args.validate_only:
        print('PASS: SPL8 and RCA each have 24 observations, eight temperatures, three replicates per temperature')
        return 0
    groups, points, stats = summarize(spl, rca)
    # PDF QA dependencies are checked only when the optional tools run.
    args.output_dir.mkdir(parents=True, exist_ok=True)
    qa = args.output_dir / f'{PREFIX}_QA'
    qa.mkdir(exist_ok=True)
    run_audit('source-preflight.txt', [QA_TOOLS / 'validate_figure.py', __file__], qa)
    for name, table in [('groups', groups), ('points', points), ('statistics', stats), ('spl8-source', spl)]:
        table.to_csv(args.output_dir / f'{PREFIX}.{name}.tsv', sep='\t', index=False, na_rep='NA', float_format='%.12g')
    fig, panels = make_figure(groups, points, stats, qa)
    stem = args.output_dir / PREFIX
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
    pd.DataFrame(panels).to_csv(qa / 'panels.tsv', sep='\t', index=False)
    report = dict(status=("automated_QA_INCOMPLETE" if _OPTIONAL_QA_SKIPPED else "automated_QA_PASS"), panels=panels, figure_size_mm=[183,150], raster_dpi=600,
                  excluded_observations=0, sample_pairing='temperature only; no individual pairing inferred',
                  transforms='ln1p each expression replicate before mean and sample SD; percentages unchanged',
                  spearman='8 means; average tied ranks; exhaustive 40320 two-sided pairing permutations',
                  multiplicity='BH across finite tests in the five-metric family',
                  x_axis='eight equally spaced groups ordered by SPL8 mean, ties by temperature; x_plot is categorical',
                  uncertainty='RCA sample SD only; SPL8 SD retained in groups table, not drawn',
                  trend='means connected in ascending SPL8 mean order; descriptive only',
                  interpretation='descriptive across temperatures; no causal or regulatory inference; proportions complementary',
                  input_sha256={k:hashlib.sha256(p.read_bytes()).hexdigest() for k,p in
                                [('expression',args.expression_file),('rca',args.rca_file)]},
                  versions={'python':sys.version.split()[0], 'numpy':np.__version__, 'pandas':pd.__version__, 'matplotlib':matplotlib.__version__})
    (qa / 'summary.json').write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    (qa / 'README.md').write_text(
        '# SPL8 and RCA associations\n\n'
        'Panels a–c show abundance and its components; d–e show complementary proportions.\n'
        'All 24 RCA observations per panel appear in eight equally spaced temperature groups ordered by SPL8 mean.\n'
        'Tick labels show measured SPL8 means (3 decimals); x_plot is the categorical group position, not a measured value.\n'
        'Equal means retain separate groups, with temperature as the sorting tie-breaker.\n'
        'They are not individually matched SPL8/RCA measurements. Temperature is the only joining key.\n'
        'Vertical bars show RCA sample SD (ddof=1, three observations); SPL8 SD is retained in the groups table.\n'
        'Expression is transformed before aggregation. The lines connect means descriptively; there is no regression.\n'
        'All correlations use eight means and exact two-sided permutation probabilities. BH q is reported separately.\n'
        'Spearman uses measured group means, never categorical positions. No exclusions or imputation.\n'
        'These are descriptive associations, not evidence of regulation.\n'
        'Source, alignment, PDF text and collision audit results accompany this report.\n', encoding='utf-8')
    export_quadratic(groups, points, args.output_dir)
    print(stats[['metric','rho','p','q_bh']].to_string(index=False))
    print(f'Exported PDF/SVG/PNG/TIFF to {args.output_dir}; automated QA {"INCOMPLETE" if _OPTIONAL_QA_SKIPPED else "PASS"}')
    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (ValueError, OSError, RuntimeError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        raise SystemExit(1)
