#!/usr/bin/env bash

# 下载 FJABFC20250784-07 RNA-seq 双端测序数据及 MD5 文件，并校验完整性。
# 分享链接有效期至 2026-10-01 10:35:41（北京时间）；过期后需替换 AUTH_QUERY。

set -Eeuo pipefail

WORK_DIR="/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/Arabidopsis_thaliana/list"
OUT_DIR="${WORK_DIR}/rnadata"
BUCKET_URL="https://ss-berry-biotech.tos-cn-beijing.volces.com"
OBJECT_PREFIX="FJABFC20250784/FJABFC20250784-07/FJABFC20250784-07"
AUTH_QUERY='X-Tos-Algorithm=TOS4-HMAC-SHA256&X-Tos-Credential=AKLTYzQzMzY0YTRmYWJiNDZkMjhhMzExMjg1MGI1YTZhY2Q%2F20260901%2Ftos-cn-beijing.volces.com%2Ftos%2Frequest&X-Tos-Date=20260901T023541Z&X-Tos-Expires=2592000&X-Tos-Policy=eyJjb25kaXRpb25zIjpbWyJzdGFydHMtd2l0aCIsIiRrZXkiLCJGSkFCRkMyMDI1MDc4NC9GSkFCRkMyMDI1MDc4NC0wNy8iXSxbImVxIiwiJGJ1Y2tldCIsInNzLWJlcnJ5LWJpb3RlY2giXV19&X-Tos-Signature=09bfd343053ebda362c4e64ad90cff59de050254ab0976bb998d32abe91ba7a4'
LINK_EXPIRES_EPOCH=1790822141

samples=(
    T4_2 T4_4 T4_5
    T8_1 T8_2 T8_4
    T12_2 T12_3 T12_6
    T16_3 T16_5 T16_6
    T20_1 T20_2 T20_4
    T24_1 T24_2 T24_4
    T28_1 T28_4 T28_6
    T32_1 T32_3 T32_4
)

# 对象清单中各 MD5.txt 的 ETag。文件均为单段对象，因此 ETag 即文件 MD5。
declare -A md5_list_hash=(
    [T4_2]="953718fe4924c7db384b93e8a5d78081"
    [T4_4]="0e1b5e441c4eaa16cc41a9b960ed9a53"
    [T4_5]="059a361219461cdc617879da34d415ff"
    [T8_1]="aaf2d70365daa820ccdd010f936e2b6a"
    [T8_2]="0a8cd2c9f7e674470d7d03b285322616"
    [T8_4]="c3da20599ecfac6b1092bbda0b7a8ef5"
    [T12_2]="a43750b5101e4014a99cccd941276356"
    [T12_3]="7349294df0ca7b0e8366424154b594cd"
    [T12_6]="50593603e89249f1e42b2334bdeb73b9"
    [T16_3]="6dcb6abc8665de066008d124c7a045b2"
    [T16_5]="d362e1d474b14c2455647cc61c82bd05"
    [T16_6]="b3fe10bd2c8e82097d3b1be0894749b5"
    [T20_1]="fc7567242d4f6efd8497056cb2f263fc"
    [T20_2]="507e04c6eb73e92beeb07cc5d8f2285a"
    [T20_4]="2eae32027fbcf3c2eb133e6015a91ac5"
    [T24_1]="ded8693c7b6c2d0ae90860237763ab35"
    [T24_2]="14c846e49235c8e6af2948beb1c964ed"
    [T24_4]="e5f7101cac42abf306b9695c5f0ace45"
    [T28_1]="d705eb451c9461d9c06a72e3d55dff02"
    [T28_4]="6cde2091d850ab600e72a386ac450edf"
    [T28_6]="99be7e4d22fcd37744e5803a1e897786"
    [T32_1]="f13209ab5cf92d730d8db26dc535ffde"
    [T32_3]="08dd1bd203ac6e4a1e5293da4b972325"
    [T32_4]="d65304fb1d06901426f2ef64fa11c3b7"
)

log() {
    printf '[%(%F %T)T] %s\n' -1 "$*"
}

die() {
    printf '错误：%s\n' "$*" >&2
    exit 1
}

file_matches_md5() {
    local file=$1
    local expected=$2

    [[ -f "$file" ]] || return 1
    printf '%s  %s\n' "$expected" "$file" | md5sum --check --status -
}

download_raw() {
    local url=$1
    local output=$2

    case "$DOWNLOADER" in
        aria2c)
            aria2c \
                --continue=true \
                --max-connection-per-server=16 \
                --split=16 \
                --min-split-size=1M \
                --file-allocation=none \
                --retry-wait=30 \
                --max-tries=10 \
                --connect-timeout=60 \
                --timeout=600 \
                --allow-overwrite=true \
                --auto-file-renaming=false \
                --console-log-level=notice \
                --summary-interval=10 \
                --dir="$(dirname "$output")" \
                --out="$(basename "$output")" \
                "$url"
            ;;
        wget)
            wget \
                --continue \
                --tries=10 \
                --timeout=600 \
                --retry-connrefused \
                --output-document="$output" \
                "$url"
            ;;
        curl)
            curl \
                --fail \
                --location \
                --continue-at - \
                --retry 10 \
                --retry-delay 30 \
                --retry-all-errors \
                --connect-timeout 60 \
                --max-time 0 \
                --output "$output" \
                "$url"
            ;;
    esac
}

download_and_verify() {
    local url=$1
    local output=$2
    local expected=$3
    local temporary="${output}.part"
    local attempt

    if file_matches_md5 "$output" "$expected"; then
        log "完整性校验通过，跳过：$output"
        return 0
    fi

    for attempt in 1 2; do
        log "下载（第 ${attempt} 次）：$output"
        download_raw "$url" "$temporary"

        if file_matches_md5 "$temporary" "$expected"; then
            mv -f -- "$temporary" "$output"
            log "完整性校验通过：$output"
            return 0
        fi

        log "MD5 不匹配，删除临时文件后重试：$temporary"
        rm -f -- "$temporary" "${temporary}.aria2"
    done

    die "两次下载后 MD5 仍不匹配：$output"
}

read_expected_md5() {
    local md5_file=$1
    local filename=$2

    awk -v target="$filename" '
        NF >= 2 {
            name = $2
            sub(/^\*/, "", name)
            sub(/^\.\//, "", name)
            if (name == target) {
                print tolower($1)
                found++
            }
        }
        END {
            if (found != 1) {
                exit 1
            }
        }
    ' "$md5_file"
}

make_url() {
    local sample=$1
    local filename=$2

    printf '%s/%s/%s/%s?%s&response-content-type=application%%2Foctet-stream' \
        "$BUCKET_URL" "$OBJECT_PREFIX" "$sample" "$filename" "$AUTH_QUERY"
}

command -v md5sum >/dev/null 2>&1 || die "未找到 md5sum"

if command -v aria2c >/dev/null 2>&1; then
    DOWNLOADER=aria2c
elif command -v wget >/dev/null 2>&1; then
    DOWNLOADER=wget
elif command -v curl >/dev/null 2>&1; then
    DOWNLOADER=curl
else
    die "未找到下载工具；请安装 aria2c、wget 或 curl"
fi

if (( $(date +%s) >= LINK_EXPIRES_EPOCH )); then
    die "分享链接已于 2026-10-01 10:35:41（北京时间）过期，请更新链接和 AUTH_QUERY"
fi

mkdir -p -- "$OUT_DIR"
log "下载工具：$DOWNLOADER"
log "输出目录：$OUT_DIR"

for sample in "${samples[@]}"; do
    sample_dir="${OUT_DIR}/${sample}"
    md5_file="${sample_dir}/MD5.txt"
    mkdir -p -- "$sample_dir"

    log "开始处理样本：$sample"
    download_and_verify \
        "$(make_url "$sample" MD5.txt)" \
        "$md5_file" \
        "${md5_list_hash[$sample]}"

    for read_number in 1 2; do
        filename="${sample}_${read_number}.fq.gz"
        if ! expected_md5=$(read_expected_md5 "$md5_file" "$filename"); then
            die "${md5_file} 中缺少或重复记录：$filename"
        fi
        [[ "$expected_md5" =~ ^[0-9a-f]{32}$ ]] || \
            die "${md5_file} 中的 MD5 格式错误：$filename"

        download_and_verify \
            "$(make_url "$sample" "$filename")" \
            "${sample_dir}/${filename}" \
            "$expected_md5"
    done

    # 严格检查 MD5.txt 中的全部条目，避免遗漏额外记录或格式错误。
    (
        cd "$sample_dir"
        md5sum --check --strict MD5.txt
    ) || die "样本完整性校验失败：$sample"
done

log "全部 24 个样本、48 个 FASTQ 文件均已通过 MD5 完整性校验。"
