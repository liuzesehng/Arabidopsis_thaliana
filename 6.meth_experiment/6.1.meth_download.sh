#!/usr/bin/env bash
set -Eeuo pipefail

WORK_DIR="${WORK_DIR:-/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/Arabidopsis_thaliana/list}"
export WORK_DIR
command -v python3 >/dev/null 2>&1 || { printf '错误：需要 python3\n' >&2; exit 1; }
command -v curl >/dev/null 2>&1 || { printf '错误：需要 curl\n' >&2; exit 1; }
command -v md5sum >/dev/null 2>&1 || { printf '错误：需要 md5sum\n' >&2; exit 1; }

python3 - "$@" <<'PY'
import base64
import datetime
import fcntl
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import subprocess
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET


WORK_DIR = Path(os.environ['WORK_DIR']).expanduser().resolve()
OUT_DIR = WORK_DIR / 'methdata'
BUCKET = 'https://ss-berry-biotech.tos-cn-beijing.volces.com'
SHARES = (
    ('08', '20260907T062509Z', 'a6a2823389418f1b2785e1786ca4636b7879f8b7870e1f9122feb72776c8a17f'),
    ('09', '20260907T062611Z', '11bdd724b431852d7a557569400532bea41eded51c01c76f9be3903bd62b5ee8'),
)


def log(message):
    print(message, flush=True)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def safe_name(name):
    require(bool(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.+-]*', name)),
            f'不安全或不支持的文件/样本名：{name!r}')
    return name


def fastq(name):
    return name.lower().endswith(('.fq.gz', '.fastq.gz', '.fq', '.fastq'))


def share_auth(number, timestamp, signature):
    prefix = f'FJABFC20250784/FJABFC20250784-{number}/'
    policy = {'conditions': [['starts-with', '$key', prefix],
                             ['eq', '$bucket', 'ss-berry-biotech']]}
    query = urllib.parse.urlencode({
        'X-Tos-Algorithm': 'TOS4-HMAC-SHA256',
        'X-Tos-Credential': 'AKLTYzQzMzY0YTRmYWJiNDZkMjhhMzExMjg1MGI1YTZhY2Q/20260907/tos-cn-beijing.volces.com/tos/request',
        'X-Tos-Date': timestamp,
        'X-Tos-Expires': '2592000',
        'X-Tos-Policy': base64.b64encode(json.dumps(policy, separators=(',', ':')).encode()).decode(),
        'X-Tos-Signature': signature,
    })
    expires = datetime.datetime.strptime(timestamp, '%Y%m%dT%H%M%SZ').replace(
        tzinfo=datetime.timezone.utc) + datetime.timedelta(seconds=2592000)
    require(datetime.datetime.now(datetime.timezone.utc) < expires,
            f'分享 {number} 已过期（{expires.isoformat()}），请更新分享凭据')
    return prefix, query


def parse_listing(payload):
    if payload.lstrip().startswith(b'{'):
        return json.loads(payload)
    root = ET.fromstring(payload)
    for element in root.iter():
        element.tag = element.tag.rsplit('}', 1)[-1]
    require(root.tag == 'ListBucketResult', '对象清单 XML 格式不正确')
    return {
        'Contents': [{child.tag: child.text for child in entry}
                     for entry in root.findall('Contents')],
        'IsTruncated': root.findtext('IsTruncated'),
        'NextMarker': root.findtext('NextMarker'),
    }


def list_objects(prefix, auth):
    marker = ''
    seen = set()
    objects = []
    while True:
        parameters = urllib.parse.urlencode({'prefix': prefix, 'max-keys': 1000, 'marker': marker})
        url = f'{BUCKET}/?{auth}&{parameters}'
        for attempt in range(3):
            try:
                with urllib.request.urlopen(url, timeout=60) as response:
                    page = parse_listing(response.read())
                break
            except Exception:
                if attempt == 2:
                    raise RuntimeError('读取远程清单失败，请检查网络、分享权限或有效期') from None
                time.sleep(5)
        require('IsTruncated' in page, '清单缺少分页状态，拒绝使用不完整清单')
        entries = page.get('Contents', []) or []
        for entry in entries:
            key = entry['Key']
            require(key.startswith(prefix) and key not in seen, '清单路径越界或重复')
            seen.add(key)
            if not key.endswith('/'):
                objects.append((entry, auth))
        truncated = page['IsTruncated']
        require(truncated in (True, False, 'true', 'false'), '无效分页状态')
        if truncated in (False, 'false'):
            return objects
        next_marker = page.get('NextMarker') or (entries[-1]['Key'] if entries else '')
        require(next_marker and next_marker > marker, '分页未推进，拒绝遗漏对象')
        marker = next_marker


def md5(path):
    digest = hashlib.md5()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def download(entry, auth, destination, expected):
    require(not destination.is_symlink(), f'拒绝符号链接：{destination}')
    size = int(entry['Size'])
    if destination.is_file() and destination.stat().st_size == size and md5(destination) == expected:
        log(f'校验通过，跳过：{destination.name}')
        return
    temporary = destination.with_name(destination.name + '.part')
    require(not temporary.is_symlink(), f'拒绝符号链接：{temporary}')
    url = f"{BUCKET}/{urllib.parse.quote(entry['Key'], safe='/')}?{auth}"
    for attempt in range(1, 4):
        if temporary.is_file() and temporary.stat().st_size >= size:
            if temporary.stat().st_size == size and md5(temporary) == expected:
                temporary.replace(destination)
                return
            temporary.unlink()
        log(f'下载 {destination.parent.name}/{destination.name}（{attempt}/3）')
        result = subprocess.run([
            'curl', '--fail', '--silent', '--show-error', '--retry', '3',
            '--retry-delay', '5', '--connect-timeout', '60', '--speed-limit', '1024',
            '--speed-time', '120', '--continue-at', '-', '--output', str(temporary), url,
        ], stderr=subprocess.PIPE)
        if result.returncode:
            if result.returncode == 33:
                temporary.unlink(missing_ok=True)
            continue
        if temporary.stat().st_size == size and md5(temporary) == expected:
            temporary.replace(destination)
            return
        log(f'MD5 或大小不匹配，清除临时文件重试：{destination.name}')
        temporary.unlink(missing_ok=True)
    raise RuntimeError(f'下载或完整性校验失败：{destination.name}；网络中断的 .part 文件保留供续传')


def parse_manifest(path, expected_names):
    checksums = {}
    for line_number, line in enumerate(path.read_text(encoding='utf-8').splitlines(), 1):
        require(bool(line), f'{path.name}:{line_number} 存在空行')
        match = re.fullmatch(r'([0-9a-fA-F]{32}) [ *](.+)', line)
        require(match is not None, f'{path.name}:{line_number} 不是标准 MD5 记录')
        checksum, name = match.groups()
        if name.startswith('./'):
            name = name[2:]
        safe_name(name)
        require(name not in checksums, f'MD5 重复记录：{name}')
        checksums[name] = checksum.lower()
    require(set(checksums) == set(expected_names),
            f'MD5 与远程 FASTQ 清单不一致：缺少 {sorted(set(expected_names) - set(checksums))}；'
            f'多出 {sorted(set(checksums) - set(expected_names))}')
    return checksums


def build_inventory(objects):
    groups = {}
    for entry, auth in objects:
        remote = PurePosixPath(entry['Key'])
        if fastq(remote.name) or remote.name.lower() == 'md5.txt':
            safe_name(remote.name)
            groups.setdefault(str(remote.parent), []).append((entry, auth))
    batches = {}
    samples = {}
    for remote_dir, entries in sorted(groups.items()):
        batch = safe_name('_'.join(PurePosixPath(remote_dir).parts[-2:]))
        require(batch not in batches, f'批次标识冲突：{batch}')
        manifests = [(entry, auth) for entry, auth in entries
                     if PurePosixPath(entry['Key']).name.lower() == 'md5.txt']
        reads = {PurePosixPath(entry['Key']).name: (entry, auth)
                 for entry, auth in entries if fastq(entry['Key'])}
        require(len(reads) == sum(fastq(entry['Key']) for entry, _ in entries),
                f'{batch} 存在重复 FASTQ')
        require(len(manifests) == 1 and bool(reads), f'{batch} 缺少唯一 MD5.txt 或 FASTQ')
        entry, auth = manifests[0]
        checksum = entry.get('ETag', '').strip('"').lower()
        require(bool(re.fullmatch('[0-9a-f]{32}', checksum)),
                f'{batch} 的 MD5 文件 ETag 不是单段 MD5，无法验证校验文件')
        batches[batch] = (entry, auth, checksum, reads)
        pairs = {}
        for filename in sorted(reads):
            match = re.fullmatch(r'(.+)_L\d+_[A-Za-z0-9]+\.R([12])\.fastq\.gz', filename)
            require(match is not None, f'无法从 FASTQ 文件名识别样本：{filename}')
            sample, read_number = match.groups()
            safe_name(sample)
            pair = pairs.setdefault(sample, {})
            require(read_number not in pair, f'{sample} 存在重复 R{read_number}，停止以避免遗漏')
            pair[read_number] = filename
        for sample, pair in pairs.items():
            require(set(pair) == {'1', '2'}, f'{sample} 缺少配对 FASTQ')
            require(pair['1'].replace('.R1.', '.R2.') == pair['2'], f'{sample} 双端文件不配对')
            require(sample not in samples, f'多个批次存在重名样本：{sample}')
            samples[sample] = (batch, [pair['1'], pair['2']])
    require(len(samples) == 22, f'应有 22 个样本，实际为 {len(samples)}（{len(groups)} 个批次目录）；停止下载')
    return batches, samples


def main():
    if sys.argv[1:] == ['--help']:
        print('用法：bash 6.1.meth_download.sh [--list-only]\n'
              '依赖：python3、curl、md5sum。按 FASTQ 文件名识别 22 个样本，\n'
              '下载到 $WORK_DIR/methdata/<样本>/，并校验 MD5 文件及全部 FASTQ。\n'
              '原始批次 MD5 保存在 methdata/.manifests/；每个样本生成对应的 MD5.txt。\n'
              '--list-only 只读取远程对象清单，不下载文件、不创建目录。\n'
              '支持 .part 断点续传；已完整文件校验后跳过。分享凭据请勿公开。')
        return
    require(sys.argv[1:] in ([], ['--list-only']), '不支持的参数；请使用 --help')
    objects = []
    for number, timestamp, signature in SHARES:
        prefix, auth = share_auth(number, timestamp, signature)
        objects.extend(list_objects(prefix, auth))
    batches, samples = build_inventory(objects)
    log(f'完整清单：{len(batches)} 个批次目录，22 个样本，44 个 FASTQ')
    if sys.argv[1:] == ['--list-only']:
        for sample, (batch, filenames) in sorted(samples.items()):
            log(f'{sample}\t{batch}\t' + '\t'.join(filenames))
        return
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    os.chdir(WORK_DIR)
    require(not OUT_DIR.is_symlink(), f'拒绝符号链接：{OUT_DIR}')
    OUT_DIR.mkdir(exist_ok=True)
    lock_path = OUT_DIR / '.download.lock'
    require(not lock_path.is_symlink(), f'拒绝符号链接：{lock_path}')
    with lock_path.open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        manifest_dir = OUT_DIR / '.manifests'
        require(not manifest_dir.is_symlink(), f'拒绝符号链接：{manifest_dir}')
        manifest_dir.mkdir(exist_ok=True)
        batch_checksums = {}
        for batch, (entry, auth, checksum, reads) in sorted(batches.items()):
            manifest = manifest_dir / f'{batch}.md5.txt'
            download(entry, auth, manifest, checksum)
            batch_checksums[batch] = parse_manifest(manifest, reads)
        for sample, (batch, filenames) in sorted(samples.items()):
            directory = OUT_DIR / sample
            require(not directory.is_symlink(), f'拒绝符号链接：{directory}')
            directory.mkdir(exist_ok=True)
            manifest = directory / 'MD5.txt'
            require(not manifest.is_symlink(), f'拒绝符号链接：{manifest}')
            manifest_text = ''.join(f'{batch_checksums[batch][filename]}  {filename}\n'
                                    for filename in filenames)
            temporary_manifest = manifest.with_name('MD5.txt.part')
            require(not temporary_manifest.is_symlink(), f'拒绝符号链接：{temporary_manifest}')
            temporary_manifest.write_text(manifest_text, encoding='utf-8')
            temporary_manifest.replace(manifest)
            for filename in filenames:
                expected = batch_checksums[batch][filename]
                read_entry, read_auth = batches[batch][3][filename]
                download(read_entry, read_auth, directory / filename, expected)
            subprocess.run(['md5sum', '--check', '--strict', manifest.name],
                           cwd=directory, check=True)
            log(f'样本全部校验通过：{sample}')
    log(f'全部 22 个样本下载及 MD5 完整性检验通过：{OUT_DIR}')


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit('已中断；保留 .part 文件供下次续传')
    except Exception as error:
        sys.exit(f'错误：{error}')
PY
