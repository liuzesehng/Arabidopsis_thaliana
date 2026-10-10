#!/usr/bin/env python3
"""按前三列匹配 meth 的实际行号，为 meth_experiment 记录添加编号。

默认输出 meth_experiment/{CG,CHG,CHH}.CV.before.numbered.tsv（无表头）。
保留目标文件原始列及顺序，匹配位点使用参考文件行号。
开头 n 个未匹配位点编号为 1-n 到 0；其余未匹配位点在前一个匹配
编号后添加 .1、.2 等后缀（包括末尾位点）。后缀是字符串，不是小数。
相同位点重复出现时复用编号；非空目标文件完全无匹配时报告错误。
新增编号列后不再是标准 bedGraph 格式，因此使用 .tsv 后缀。
"""

import argparse
from pathlib import Path
import sys


DEFAULT_LIST_DIR = Path(
    "/datapool/home/2023102768/lico_share_dir/life-gongl/"
    "zesheng/Arabidopsis_thaliana/list"
)
SUFFIX = ".CV.before.bedGraph"


def read_records(path):
    """忽略空行和 bedGraph 说明行，但行号始终按物理行计数。"""
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            raw = line.rstrip("\r\n")
            fields = raw.split()
            if not fields or raw.lstrip().startswith("#"):
                continue
            if fields[0] in {"track", "browser"}:
                continue
            if len(fields) < 4:
                raise ValueError(f"{path}:{line_number}: 至少需要四列")
            yield line_number, tuple(fields[:3]), raw


def number_records(reference, target):
    prefix = reference.name.removesuffix(SUFFIX)
    index = {}
    for line_number, key, _ in read_records(reference):
        if key in index:
            raise ValueError(
                f"{reference}:{line_number}: 前三列重复 {key}，"
                f"与第 {index[key]} 行冲突，无法唯一编号"
            )
        index[key] = line_number

    records = list(read_records(target))
    first_match = next(
        (i for i, (_, key, _) in enumerate(records) if key in index), None
    )
    if records and first_match is None:
        raise ValueError(f"{target}: 没有任何匹配位点，无法按相邻匹配编号插入")

    # 开头未匹配位点以 0 结尾；重复坐标仅占一个编号。
    leading_keys = dict.fromkeys(key for _, key, _ in records[:first_match])
    assigned = {
        key: f"{prefix}_{number}"
        for number, key in enumerate(leading_keys, 1 - len(leading_keys))
    }
    last_match = None
    suffix_counts = {}
    result = []
    missing = 0
    for line_number, key, raw in records:
        reference_line = index.get(key)
        if reference_line is None:
            if key not in assigned:
                suffix_counts[last_match] = suffix_counts.get(last_match, 0) + 1
                assigned[key] = f"{prefix}_{last_match}.{suffix_counts[last_match]}"
            identifier = assigned[key]
            missing += 1
            print(
                f"提示：{target}:{line_number}: 未匹配 {key}，插入编号 {identifier}",
                file=sys.stderr,
            )
        else:
            last_match = reference_line
            identifier = f"{prefix}_{reference_line}"
        result.append(f"{identifier}\t{raw}\n")
    return result, missing


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--list-dir", type=Path, default=DEFAULT_LIST_DIR,
        help="包含 meth 和 meth_experiment 的目录",
    )
    parser.add_argument(
        "--output-dir", type=Path,
        help="输出目录；默认使用 list-dir/meth_experiment",
    )
    args = parser.parse_args()
    output_dir = args.output_dir or args.list_dir / "meth_experiment"

    # 先检查全部输入，避免格式错误时仅生成部分结果。
    prepared = []
    try:
        for context in ("CG", "CHG", "CHH"):
            name = context + SUFFIX
            rows, missing = number_records(
                args.list_dir / "meth" / name,
                args.list_dir / "meth_experiment" / name,
            )
            prepared.append((context, rows, missing))
        output_dir.mkdir(parents=True, exist_ok=True)
        for context, rows, missing in prepared:
            output = output_dir / f"{context}.CV.before.numbered.tsv"
            output.write_text("".join(rows), encoding="utf-8")
            print(
                f"{context}: 总行数={len(rows)}, "
                f"匹配={len(rows) - missing}, 未匹配={missing} -> {output}"
            )
    except (OSError, ValueError) as exc:
        parser.exit(1, f"错误：{exc}\n")


if __name__ == "__main__":
    main()
