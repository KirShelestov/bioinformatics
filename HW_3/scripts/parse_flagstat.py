#!/usr/bin/env python3
"""
parse_flagstat.py — разбор результатов samtools flagstat
Домашнее задание 3, НГУ 2026

Использование:
    python3 parse_flagstat.py results/sample.flagstat.txt
    python3 parse_flagstat.py results/sample.flagstat.txt --threshold 90
"""

import re
import sys
import argparse
from pathlib import Path


def parse_flagstat(filepath: str) -> dict:
    """Парсит файл samtools flagstat и возвращает словарь с метриками."""
    metrics = {}
    path = Path(filepath)

    if not path.exists():
        print(f"[ERROR] Файл не найден: {filepath}", file=sys.stderr)
        sys.exit(1)

    with open(path) as f:
        for line in f:
            line = line.strip()

            # Пример строки: "1099585 + 0 in total (QC-passed reads + QC-failed reads)"
            m = re.match(r"(\d+) \+ (\d+) (.+)", line)
            if not m:
                continue

            qc_pass = int(m.group(1))
            qc_fail = int(m.group(2))
            label   = m.group(3)

            # Извлекаем процент, если есть: (99.83% : N/A)
            pct_match = re.search(r"\((\d+\.?\d*)%", label)
            pct = float(pct_match.group(1)) if pct_match else None

            # Чистый ключ без скобок
            key = re.sub(r"\s*\(.*\)", "", label).strip()

            metrics[key] = {
                "qc_pass": qc_pass,
                "qc_fail": qc_fail,
                "percent": pct,
            }

    return metrics


def get_mapped_percent(metrics: dict) -> float | None:
    """Возвращает процент картированных ридов."""
    for key, val in metrics.items():
        if key.startswith("mapped"):
            return val["percent"]
    return None


def print_report(metrics: dict, threshold: float) -> None:
    """Печатает отчёт и оценку качества."""
    print("=" * 55)
    print("  samtools flagstat — отчёт")
    print("=" * 55)

    for key, val in metrics.items():
        pct_str = f"  ({val['percent']}%)" if val["percent"] is not None else ""
        print(f"  {val['qc_pass']:>10}  {key}{pct_str}")

    print("=" * 55)

    mapped_pct = get_mapped_percent(metrics)
    if mapped_pct is None:
        print("[WARN] Не удалось найти строку 'mapped' в файле.")
        return

    print(f"\n  % картированных ридов: {mapped_pct}%")
    print(f"  Порог:                 {threshold}%")

    if mapped_pct > threshold:
        print(f"\n  ✓ Mapping quality: OK")
    else:
        print(f"\n  ✗ Mapping quality: not OK")

    print()


def main():
    parser = argparse.ArgumentParser(
        description="Парсер samtools flagstat"
    )
    parser.add_argument("flagstat_file", help="Путь к файлу flagstat")
    parser.add_argument(
        "--threshold", type=float, default=90.0,
        help="Порог %% картирования (по умолчанию: 90)"
    )
    parser.add_argument(
        "--json", action="store_true",
        help="Вывод в формате JSON"
    )
    args = parser.parse_args()

    metrics = parse_flagstat(args.flagstat_file)

    if args.json:
        import json
        mapped_pct = get_mapped_percent(metrics)
        result = {
            "file": args.flagstat_file,
            "mapped_percent": mapped_pct,
            "threshold": args.threshold,
            "status": "OK" if (mapped_pct or 0) > args.threshold else "not OK",
            "metrics": metrics,
        }
        print(json.dumps(result, ensure_ascii=False, indent=2))
    else:
        print_report(metrics, args.threshold)


if __name__ == "__main__":
    main()
