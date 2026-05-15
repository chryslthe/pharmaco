"""
PharmacoVEP full local pipeline for PyCharm.

Research prototype for a pharmacogenetic hackathon.
Not a medical device and not a prescription.

Expected input files in the same folder as this script by default:
    dataset_full.csv
    PharmaVEP_final.csv
    1kG_Full_dataset.vcf

Run:
    python pharmavep_pycharm.py

Outputs:
    pharmavep_outputs/
"""

from __future__ import annotations

import gzip
import json
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from jinja2 import Template
from pydantic import BaseModel
pip install -r requirements_pycharm.txt
python pharmavep_pycharm.py

# =============================================================================
# 1. Paths
# =============================================================================
# Если входные файлы лежат рядом с этим скриптом, ничего менять не нужно.
# Если файлы лежат в другой папке, поменяй BASE_DIR, например:
# BASE_DIR = Path("/Users/arinabolsuhina/Downloads")

BASE_DIR = Path(pharmavep).resolve().parent

RULES_PATH = BASE_DIR / "dataset_full.csv"
VCF_PATH = BASE_DIR / "1kG_Full_dataset.vcf"
ANNOTATION_CSV_PATH = BASE_DIR / "PharmaVEP_final.csv"

OUT_DIR = BASE_DIR / "pharmavep_outputs"
OUT_DIR.mkdir(parents=True, exist_ok=True)

ASSEMBLY = "GRCh38"


# =============================================================================
# 2. Basic helpers
# =============================================================================

RSID_RE = re.compile(r"rs\d+")


def normalize_chrom(chrom: str) -> str:
    chrom = str(chrom).strip()
    return chrom if chrom.startswith("chr") else "chr" + chrom


def make_variant_key(chrom: str, pos: int, ref: str, alt: str) -> str:
    return f"{normalize_chrom(chrom)}:{int(pos)}:{str(ref).upper()}:{str(alt).upper()}"


def extract_rule_rsids(value) -> list[str]:
    if pd.isna(value):
        return []
    return RSID_RE.findall(str(value))


def normalize_allele_name(value) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    match = re.search(r"\*[\w]+", text)
    return match.group(0) if match else text


def check_input_files() -> None:
    print("Input paths")
    print("RULES:", RULES_PATH, RULES_PATH.exists())
    print("VCF:", VCF_PATH, VCF_PATH.exists())
    print("ANNOTATION:", ANNOTATION_CSV_PATH, ANNOTATION_CSV_PATH.exists())
    print("OUT_DIR:", OUT_DIR, OUT_DIR.exists())

    missing = [
        path
        for path in [RULES_PATH, VCF_PATH, ANNOTATION_CSV_PATH]
        if not path.exists()
    ]
    if missing:
        raise FileNotFoundError(
            "Missing input files:\n"
            + "\n".join(str(path) for path in missing)
            + "\n\nPut the files next to pharmavep_pycharm.py or edit BASE_DIR."
        )


# =============================================================================
# 3. Rule validation and loading
# =============================================================================


class RuleModel(BaseModel):
    Drug: str
    Gene: str
    Marker_type: str
    rsID: str | None = None
    Allele: str
    Function: str | None = None
    Phenotype: str | None = None
    Recommendation: str
    CPIC_Level: str | None = None
    Source: str | None = None
    Disease: str | None = None


def load_rules(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path).fillna("")

    required = [
        "Drug",
        "Gene",
        "Marker_type",
        "rsID",
        "Allele",
        "Function",
        "Phenotype",
        "Recommendation",
        "CPIC_Level",
        "Source",
        "Disease",
    ]
    missing = set(required) - set(df.columns)
    if missing:
        raise ValueError(f"dataset_full.csv is missing columns: {sorted(missing)}")

    errors = []
    for idx, row in df.iterrows():
        try:
            RuleModel(**row[required].to_dict())
        except Exception as exc:
            errors.append((idx, exc))

    if errors:
        print("Rule validation errors:")
        for idx, exc in errors[:10]:
            print(idx, exc)
        raise ValueError("dataset_full.csv has validation errors")

    df["rule_id"] = range(1, len(df) + 1)
    df["rsids_norm"] = df["rsID"].apply(extract_rule_rsids)
    df["allele_norm"] = df["Allele"].apply(normalize_allele_name)
    df["marker_type_norm"] = df["Marker_type"].astype(str).str.lower()
    df["is_hla"] = df["Gene"].astype(str).str.upper().eq("HLA-B")

    unique_rsids = sorted({x for xs in df["rsids_norm"] for x in xs})
    print("\nRules loaded")
    print("Rules:", len(df))
    print("Drugs:", df["Drug"].nunique())
    print("Genes:", df["Gene"].nunique())
    print("Unique rsIDs:", len(unique_rsids))
    return df


# =============================================================================
# 4. Build rsID -> coordinate cache
# =============================================================================


def build_rsid_coordinate_cache_from_annotation(
    annotation_path: Path, rules: pd.DataFrame, cache_path: Path
) -> pd.DataFrame:
    ann = pd.read_csv(annotation_path).fillna("")

    required = {"CHROM", "POS", "REF", "ALT"}
    missing = required - set(ann.columns)
    if missing:
        raise ValueError(f"PharmaVEP_final.csv is missing columns: {sorted(missing)}")

    rsid_columns = [col for col in ["rsID", "ID"] if col in ann.columns]
    if not rsid_columns:
        raise ValueError("PharmaVEP_final.csv must have an rsID or ID column")

    wanted_rsids = sorted({rsid for ids in rules["rsids_norm"] for rsid in ids})
    rows = []

    for _, row in ann.iterrows():
        row_rsids = set()
        for col in rsid_columns:
            row_rsids.update(RSID_RE.findall(str(row.get(col, ""))))

        matched = row_rsids.intersection(wanted_rsids)
        if not matched:
            continue

        chrom = normalize_chrom(row["CHROM"])
        pos = int(row["POS"])
        ref = str(row["REF"]).upper()

        for alt in str(row["ALT"]).split(","):
            alt = alt.strip().upper()
            if not alt:
                continue
            key = make_variant_key(chrom, pos, ref, alt)

            for rsid in matched:
                rows.append(
                    {
                        "rsid": rsid,
                        "chrom": chrom,
                        "pos": pos,
                        "ref": ref,
                        "alt": alt,
                        "key": key,
                        "status": "resolved",
                        "source": annotation_path.name,
                    }
                )

    found = {row["rsid"] for row in rows}
    for rsid in wanted_rsids:
        if rsid not in found:
            rows.append(
                {
                    "rsid": rsid,
                    "chrom": "",
                    "pos": "",
                    "ref": "",
                    "alt": "",
                    "key": "",
                    "status": "not_resolved",
                    "source": annotation_path.name,
                }
            )

    cache = pd.DataFrame(rows).drop_duplicates()
    cache.to_csv(cache_path, index=False)

    print("\nCoordinate cache")
    print(cache["status"].value_counts())
    unresolved = cache[cache["status"].eq("not_resolved")]
    if len(unresolved):
        print("Not resolved rsIDs:", sorted(unresolved["rsid"].tolist()))

    return cache


# =============================================================================
# 5. HLA-B proxy SNP
# =============================================================================


HLA_PROXY_RULES = {
    "HLA-B*58:01": {
        "proxy_snp": "rs9263726",
        "risk_allele": "A",
        "drug": "Allopurinol",
        "gene": "HLA-B",
        "phenotype": "High risk of severe cutaneous adverse reactions (SCAR)",
        "recommendation": (
            "Обнаружен proxy SNP, связанный с HLA-B*58:01: возможен повышенный риск "
            "тяжелых кожных нежелательных реакций на аллопуринол. "
            "Аллопуринол следует избегать до подтверждения прямым HLA-типированием."
        ),
        "confidence": "LOW",
        "method_note": (
            "HLA-B*58:01 оценен через proxy SNP rs9263726, а не прямым HLA typing. "
            "Связь proxy SNP с HLA-B*58:01 зависит от популяции."
        ),
    }
}


def add_hla_proxy_coordinate(rsid_coords: pd.DataFrame, coord_cache_path: Path) -> pd.DataFrame:
    hla_proxy_json_path = OUT_DIR / "hla_proxy_rules.json"
    with open(hla_proxy_json_path, "w", encoding="utf-8") as handle:
        json.dump(HLA_PROXY_RULES, handle, ensure_ascii=False, indent=2)

    # Manual GRCh38 coordinate for rs9263726.
    manual_hla_proxy_coords = pd.DataFrame(
        [
            {
                "rsid": "rs9263726",
                "chrom": "chr6",
                "pos": 31272120,
                "ref": "G",
                "alt": "A",
                "key": make_variant_key("chr6", 31272120, "G", "A"),
                "status": "resolved",
                "source": "manual GRCh38 coordinate for HLA proxy SNP",
            }
        ]
    )

    rsid_coords = rsid_coords[rsid_coords["rsid"].ne("rs9263726")].copy()
    rsid_coords = pd.concat([rsid_coords, manual_hla_proxy_coords], ignore_index=True)
    rsid_coords = rsid_coords.drop_duplicates()
    rsid_coords.to_csv(coord_cache_path, index=False)
    print("\nHLA proxy SNP added: rs9263726, confidence LOW")
    return rsid_coords


# =============================================================================
# 6. Stream-read target variants from VCF
# =============================================================================


@dataclass
class VcfCall:
    sample: str
    chrom: str
    pos: int
    ref: str
    alt: str
    gt: str
    alt_count: int | None
    rsid: str | None = None
    key: str | None = None


def is_gzip_file(path: Path) -> bool:
    with open(path, "rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def test_gzip(path: Path) -> bool:
    if not is_gzip_file(path):
        print("\nVCF is plain text, gzip test skipped")
        return True

    result = subprocess.run(
        ["gzip", "-t", str(path)],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode == 0:
        print("\nVCF gzip integrity: OK")
        return True

    print("\nWARNING: VCF gzip integrity check failed")
    print(result.stderr.strip())
    print("The file may be truncated. Results can be incomplete.")
    return False


def iter_vcf_lines(path: Path):
    if is_gzip_file(path):
        proc = subprocess.Popen(
            ["gzip", "-dc", str(path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
        for line in proc.stdout:
            yield line

        stderr = proc.stderr.read()
        return_code = proc.wait()
        if return_code != 0:
            print("WARNING: gzip exited with code", return_code)
            print(stderr[:1000])
    else:
        with open(path, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                yield line


def gt_alt_count_for_alt(gt: str | None, alt_index: int) -> int | None:
    if gt is None:
        return None
    alleles = str(gt).replace("|", "/").split("/")
    if "." in alleles:
        return None

    count = 0
    for allele in alleles:
        try:
            if int(allele) == alt_index:
                count += 1
        except ValueError:
            return None
    return count


def parse_sample_gt(sample_field: str, format_keys: list[str]) -> str:
    values = str(sample_field).split(":")
    data = dict(zip(format_keys, values))
    return data.get("GT", "./.")


def read_vcf_calls_filtered(
    vcf_path: Path, rsid_coords_df: pd.DataFrame
) -> tuple[dict, dict, pd.DataFrame, list[str]]:
    resolved = rsid_coords_df[rsid_coords_df["status"].eq("resolved")].copy()
    key_to_rsids = (
        resolved.groupby("key")["rsid"].apply(lambda x: sorted(set(x))).to_dict()
    )
    target_keys = set(key_to_rsids.keys())

    samples = []
    by_rsid = {}
    by_coord = {}
    variant_rows = []
    total_seen = 0
    total_kept = 0
    malformed_lines = 0

    for line in iter_vcf_lines(vcf_path):
        line = line.rstrip("\n")
        if not line or line.startswith("##"):
            continue

        if line.startswith("#CHROM"):
            header = line.split("\t")
            samples = header[9:]
            by_rsid = {sample: {} for sample in samples}
            by_coord = {sample: {} for sample in samples}
            continue

        if line.startswith("#"):
            continue

        total_seen += 1
        parts = line.split("\t")
        if len(parts) < 10:
            malformed_lines += 1
            continue

        chrom, pos, variant_id, ref, alt_string = (
            parts[0],
            parts[1],
            parts[2],
            parts[3],
            parts[4],
        )
        chrom = normalize_chrom(chrom)
        try:
            pos = int(pos)
        except ValueError:
            malformed_lines += 1
            continue

        ref = str(ref).upper()
        alts = [item.upper() for item in str(alt_string).split(",")]
        format_keys = parts[8].split(":")
        sample_fields = parts[9:]

        for alt_index, alt in enumerate(alts, start=1):
            key = make_variant_key(chrom, pos, ref, alt)
            if key not in target_keys:
                continue

            total_kept += 1
            rsids = key_to_rsids.get(key, [])
            variant_rows.append(
                {
                    "chrom": chrom,
                    "pos": pos,
                    "id": variant_id,
                    "ref": ref,
                    "alt": alt,
                    "key": key,
                    "mapped_rsids": ";".join(rsids),
                }
            )

            for sample, sample_field in zip(samples, sample_fields):
                gt = parse_sample_gt(sample_field, format_keys)
                alt_count = gt_alt_count_for_alt(gt, alt_index)
                call = VcfCall(
                    sample=sample,
                    chrom=chrom,
                    pos=pos,
                    ref=ref,
                    alt=alt,
                    gt=gt,
                    alt_count=alt_count,
                    key=key,
                )
                by_coord[sample][key] = call
                for rsid in rsids:
                    by_rsid[sample][rsid] = call

    variant_index = pd.DataFrame(variant_rows).drop_duplicates()

    print("\nVCF scan")
    print("Samples in VCF:", len(samples))
    print("Variants scanned:", total_seen)
    print("Target variants found:", total_kept)
    print(
        "Unique target coordinates found:",
        variant_index["key"].nunique() if len(variant_index) else 0,
    )
    print("Malformed lines skipped:", malformed_lines)

    variant_index.to_csv(OUT_DIR / "found_target_variants.csv", index=False)

    missing_keys = (
        sorted(target_keys - set(variant_index["key"])) if len(variant_index) else sorted(target_keys)
    )
    if missing_keys:
        print("Missing target coordinates:", len(missing_keys))
        print(missing_keys[:50])

    return by_rsid, by_coord, variant_index, samples


# =============================================================================
# 7. Marker interpretation
# =============================================================================


STANDARD_NO_PGX_ACTION_RU = (
    "Фармакогенетических оснований для изменения терапии по оцененным маркерам не выявлено; "
    "используйте стандартный режим с обычным клиническим мониторингом."
)

INCOMPLETE_ASSESSMENT_RU = (
    "Оценка неполная: часть клинически значимых маркеров не найдена в доступных данных; "
    "отрицательный результат следует интерпретировать осторожно."
)

REFERENCE_WORDS = {
    "negative",
    "absent",
    "wild-type",
    "wild type",
    "reference",
    "normal",
    "no alternate",
    "no risk marker detected",
}


def is_reference_or_absent_rule(row) -> bool:
    allele = str(row.get("Allele", "")).lower()
    phenotype = str(row.get("Phenotype", "")).lower()
    function = str(row.get("Function", "")).lower()
    text = " ".join([allele, phenotype, function])
    if row.get("allele_norm") == "*1":
        return True
    return any(word in text for word in REFERENCE_WORDS)


def deduplicate_markers(markers: list[dict]) -> list[dict]:
    seen = set()
    output = []
    for marker in markers:
        key = (
            marker.get("rsid"),
            marker.get("allele"),
            marker.get("gt"),
            marker.get("key"),
        )
        if key in seen:
            continue
        seen.add(key)
        output.append(marker)
    return output


def genotype_alleles_from_call(call: VcfCall) -> list[str]:
    alleles = []
    for part in re.split(r"[\/|]", str(call.gt)):
        if part == ".":
            alleles.append(".")
        elif part == "0":
            alleles.append(call.ref.upper())
        elif part == "1":
            alleles.append(call.alt.upper())
        else:
            alleles.append("?")
    return alleles


def infer_marker_status(sample: str, group: pd.DataFrame, vcf_by_rsid: dict):
    drug = group["Drug"].iloc[0]
    gene = group["Gene"].iloc[0]
    gene_upper = str(gene).upper()
    detected = []
    missing = []
    warnings = []

    if gene_upper == "HLA-B":
        for hla_allele, proxy_rule in HLA_PROXY_RULES.items():
            if str(proxy_rule["drug"]).lower() != str(drug).lower():
                continue

            proxy_snp = proxy_rule["proxy_snp"]
            risk_allele = proxy_rule["risk_allele"].upper()
            call = vcf_by_rsid.get(sample, {}).get(proxy_snp)
            if call is None:
                return (
                    "not_assessed",
                    [],
                    [proxy_snp],
                    [
                        "HLA-B не удалось оценить: proxy SNP не найден в VCF.",
                        "Для клинического вывода нужен прямой HLA typing.",
                    ],
                )

            genotype_alleles = genotype_alleles_from_call(call)
            risk_detected = risk_allele in genotype_alleles
            warnings.append("HLA оценен через proxy SNP, а не прямым HLA typing.")
            warnings.append(proxy_rule["method_note"])

            if risk_detected:
                detected.append(
                    {
                        "rsid": proxy_snp,
                        "allele": hla_allele + " proxy",
                        "allele_norm": hla_allele,
                        "function": "Proxy marker for HLA risk allele",
                        "phenotype": proxy_rule["phenotype"],
                        "recommendation": proxy_rule["recommendation"],
                        "gt": call.gt,
                        "copies": call.alt_count,
                        "key": call.key,
                        "marker_type": "HLA proxy SNP",
                    }
                )
                return (
                    f"{hla_allele} proxy-positive via {proxy_snp} ({'/'.join(genotype_alleles)})",
                    detected,
                    [],
                    warnings,
                )

            return (
                f"{hla_allele} proxy-negative via {proxy_snp} ({'/'.join(genotype_alleles)})",
                [],
                [],
                warnings,
            )

    if gene_upper == "CYP2D6" and group["marker_type_norm"].str.contains(
        "copy number|variation"
    ).any():
        warnings.append(
            "CYP2D6 CNV/duplication оценивается неполно: нужен отдельный CYP2D6 CNV/star-allele caller."
        )
    if gene_upper == "CFTR":
        warnings.append(
            "CFTR gating mutations оцениваются панельно; одиночный rsID не покрывает все клинически значимые варианты."
        )
    if gene_upper == "NAT2":
        warnings.append("NAT2 фенотип обычно гаплотипный; текущая оценка упрощена по доступным SNP.")
    if gene_upper == "G6PD":
        warnings.append("G6PD X-сцепленный ген; без информации о поле пациента интерпретация упрощена.")

    for _, row in group.iterrows():
        if not row["rsids_norm"]:
            if not row["is_hla"]:
                missing.append(str(row["rsID"]))
            continue

        for rsid in row["rsids_norm"]:
            call = vcf_by_rsid.get(sample, {}).get(rsid)
            if call is None or call.alt_count is None:
                missing.append(rsid)
                continue

            if call.alt_count > 0 and not is_reference_or_absent_rule(row):
                detected.append(
                    {
                        "rsid": rsid,
                        "allele": row["Allele"],
                        "allele_norm": row["allele_norm"],
                        "function": row["Function"],
                        "phenotype": row["Phenotype"],
                        "recommendation": row["Recommendation"],
                        "gt": call.gt,
                        "copies": call.alt_count,
                        "key": call.key,
                        "marker_type": row["Marker_type"],
                    }
                )

    detected = deduplicate_markers(detected)
    missing = sorted(set(missing))
    if detected:
        status = "; ".join(
            [f"{marker['allele']} ({marker['rsid']}, GT={marker['gt']})" for marker in detected]
        )
    else:
        status = "no alternate/risk marker detected"
    return status, detected, missing, warnings


# =============================================================================
# 8. Phenotype mapping
# =============================================================================


STAR_ALLELE_RE = re.compile(r"\*[\w]+")


def build_simple_diplotype(gene: str, detected_markers: list[dict]) -> str:
    star_hits = []
    for marker in detected_markers:
        allele = marker.get("allele_norm", "")
        copies = marker.get("copies", 0)
        if not str(allele).startswith("*"):
            continue
        for _ in range(int(copies)):
            star_hits.append(allele)

    if not star_hits:
        return "*1/*1"
    if len(star_hits) == 1:
        return f"*1/{star_hits[0]}"
    return "/".join(star_hits[:2])


def fallback_star_phenotype(gene: str, diplotype: str) -> str:
    gene = str(gene).upper()
    alleles = STAR_ALLELE_RE.findall(str(diplotype))

    if gene == "CYP2C19":
        if "*17" in alleles and not any(a in {"*2", "*3"} for a in alleles):
            return "Rapid Metabolizer" if alleles.count("*17") == 1 else "Ultrarapid Metabolizer"
        if alleles.count("*2") + alleles.count("*3") == 2:
            return "Poor Metabolizer"
        if any(a in {"*2", "*3"} for a in alleles):
            return "Intermediate Metabolizer"
        return "Normal Metabolizer"

    if gene == "CYP2C9":
        bad = sum(a in {"*2", "*3", "*5", "*6", "*8", "*11"} for a in alleles)
        return ["Normal Metabolizer", "Intermediate Metabolizer", "Poor Metabolizer"][min(bad, 2)]

    if gene == "CYP3A5":
        if alleles == ["*3", "*3"]:
            return "Poor Metabolizer"
        if "*3" in alleles:
            return "Intermediate Metabolizer"
        return "Normal Metabolizer"

    if gene == "CYP2B6":
        bad = sum(a in {"*6", "*18"} for a in alleles)
        return ["Normal Metabolizer", "Intermediate Metabolizer", "Poor Metabolizer"][min(bad, 2)]

    if gene == "TPMT":
        bad = sum(a in {"*2", "*3A", "*3B", "*3C"} for a in alleles)
        return ["Normal Metabolizer", "Intermediate Metabolizer", "Poor Metabolizer"][min(bad, 2)]

    if gene == "CYP2D6":
        bad = sum(a in {"*3", "*4", "*5", "*6"} for a in alleles)
        reduced = sum(a in {"*10", "*17", "*41"} for a in alleles)
        if bad == 2:
            return "Poor Metabolizer"
        if bad == 1 or reduced > 0:
            return "Intermediate Metabolizer"
        return "Normal Metabolizer"

    return "Normal Metabolizer" if diplotype == "*1/*1" else "unknown"


def infer_phenotype(gene: str, status: str, detected_markers: list[dict]) -> str:
    if status == "not_assessed":
        return "not assessed"

    if detected_markers:
        if any(str(marker.get("allele_norm", "")).startswith("*") for marker in detected_markers):
            diplotype = build_simple_diplotype(gene, detected_markers)
            return fallback_star_phenotype(gene, diplotype)

        phenotypes = sorted(
            {
                str(marker.get("phenotype"))
                for marker in detected_markers
                if marker.get("phenotype") and str(marker.get("phenotype")) != "nan"
            }
        )
        if phenotypes:
            return "; ".join(phenotypes)

    return "normal/reference-like or no risk marker detected"


# =============================================================================
# 9. Recommendations, outputs, and plots
# =============================================================================


def select_recommendation(group, detected_markers, phenotype, missing_markers, status) -> str:
    if status == "not_assessed":
        return (
            "Рекомендация не сформирована: генотип/фенотип не удалось надежно оценить "
            "по доступным VCF-данным. "
            + INCOMPLETE_ASSESSMENT_RU
        )

    if detected_markers:
        recommendations = []
        for marker in detected_markers:
            rec = marker.get("recommendation")
            if rec and str(rec) != "nan":
                recommendations.append(str(rec))
        if recommendations:
            return " ".join(dict.fromkeys(recommendations))

    rec = STANDARD_NO_PGX_ACTION_RU
    if missing_markers:
        rec += " " + INCOMPLETE_ASSESSMENT_RU
    return rec


def classify_actionability(status, detected_markers, missing_markers, warnings) -> str:
    if status == "not_assessed":
        return "not_assessed"
    if missing_markers or warnings:
        if detected_markers:
            return "actionable_incomplete"
        return "incomplete"
    if detected_markers:
        return "actionable"
    return "standard"


def interpret_all_samples(rules_df: pd.DataFrame, samples: list[str], vcf_by_rsid: dict) -> pd.DataFrame:
    results = []
    grouped = rules_df.groupby(["Drug", "Gene"], sort=True)

    for sample in samples:
        for (drug, gene), group in grouped:
            status, detected, missing, warnings = infer_marker_status(sample, group, vcf_by_rsid)
            phenotype = infer_phenotype(gene, status, detected)
            recommendation = select_recommendation(group, detected, phenotype, missing, status)
            actionability = classify_actionability(status, detected, missing, warnings)

            if detected and any(str(marker.get("allele_norm", "")).startswith("*") for marker in detected):
                diplotype_or_status = build_simple_diplotype(gene, detected)
            else:
                diplotype_or_status = status

            results.append(
                {
                    "sample": sample,
                    "drug": drug,
                    "gene": gene,
                    "diplotype_or_status": diplotype_or_status,
                    "phenotype": phenotype,
                    "recommendation": recommendation,
                    "actionability": actionability,
                    "detected_markers": "\n".join(
                        [
                            f"{marker['rsid']} {marker['allele']} GT={marker['gt']} copies={marker['copies']}"
                            for marker in detected
                        ]
                    )
                    if detected
                    else "—",
                    "coverage": f"covered: {len(group) - len(missing)}; missing: {len(missing)}",
                    "warnings": "\n".join(warnings) if warnings else "—",
                }
            )

    return pd.DataFrame(results)


def save_html_report(result_df: pd.DataFrame) -> Path:
    html_template = Template(
        """
<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>Personal pharmacogenetic report</title>
<style>
body {
    background: #111827;
    color: #f9fafb;
    font-family: Arial, sans-serif;
    margin: 32px;
}
h1, h2 { color: #ffffff; }
.sample { margin-bottom: 40px; }
table {
    border-collapse: collapse;
    width: 100%;
    margin-bottom: 24px;
    font-size: 13px;
}
th, td {
    border: 1px solid #374151;
    padding: 8px;
    vertical-align: top;
}
th { background: #1f2937; }
td { background: #111827; }
.actionable { color: #fecaca; font-weight: bold; }
.standard { color: #bbf7d0; }
.incomplete, .actionable_incomplete, .not_assessed { color: #fde68a; }
small { color: #d1d5db; }
</style>
</head>
<body>
<h1>Personal pharmacogenetic report</h1>
<p><small>Research prototype for hackathon. Not a medical prescription.</small></p>

{% for sample, rows in grouped %}
<div class="sample">
<h2>Sample: {{ sample }}</h2>
<table>
<thead>
<tr>
<th>Drug</th>
<th>Gene</th>
<th>Diplotype/status</th>
<th>Phenotype</th>
<th>Actionability</th>
<th>Recommendation</th>
<th>Detected markers</th>
<th>Coverage / warnings</th>
</tr>
</thead>
<tbody>
{% for r in rows %}
<tr>
<td>{{ r.drug }}</td>
<td>{{ r.gene }}</td>
<td>{{ r.diplotype_or_status }}</td>
<td>{{ r.phenotype }}</td>
<td class="{{ r.actionability }}">{{ r.actionability }}</td>
<td>{{ r.recommendation }}</td>
<td>{{ r.detected_markers | replace('\\n', '<br>') }}</td>
<td>{{ r.coverage }}<br>{{ r.warnings | replace('\\n', '<br>') }}</td>
</tr>
{% endfor %}
</tbody>
</table>
</div>
{% endfor %}
</body>
</html>
"""
    )

    grouped = [
        (sample, rows.to_dict(orient="records"))
        for sample, rows in result_df.groupby("sample", sort=True)
    ]
    html = html_template.render(grouped=grouped)
    html_path = OUT_DIR / "personal_pharmacogenetic_report_full.html"
    html_path.write_text(html, encoding="utf-8")
    return html_path


def save_plots(result_df: pd.DataFrame) -> Path:
    plot_dir = OUT_DIR / "plots"
    plot_dir.mkdir(exist_ok=True)

    plt.figure(figsize=(10, 5))
    sns.countplot(
        data=result_df,
        x="actionability",
        order=result_df["actionability"].value_counts().index,
    )
    plt.title("Actionability summary")
    plt.xticks(rotation=30, ha="right")
    plt.tight_layout()
    plt.savefig(plot_dir / "actionability_summary.png", dpi=200)
    plt.close()

    plt.figure(figsize=(14, 6))
    sns.countplot(data=result_df, x="drug", hue="actionability")
    plt.title("Actionability by drug")
    plt.xticks(rotation=45, ha="right")
    plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(plot_dir / "actionability_by_drug.png", dpi=200)
    plt.close()

    plt.figure(figsize=(14, 6))
    sns.countplot(data=result_df, x="gene", hue="actionability")
    plt.title("Actionability by gene")
    plt.xticks(rotation=45, ha="right")
    plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(plot_dir / "actionability_by_gene.png", dpi=200)
    plt.close()

    return plot_dir


def main() -> None:
    check_input_files()
    test_gzip(VCF_PATH)

    rules = load_rules(RULES_PATH)

    coord_cache_path = OUT_DIR / "rsid_coordinate_cache.csv"
    rsid_coords = build_rsid_coordinate_cache_from_annotation(
        ANNOTATION_CSV_PATH, rules, coord_cache_path
    )
    rsid_coords = add_hla_proxy_coordinate(rsid_coords, coord_cache_path)

    vcf_by_rsid, _vcf_by_coord, variant_index, samples = read_vcf_calls_filtered(
        VCF_PATH, rsid_coords
    )

    result_df = interpret_all_samples(rules, samples, vcf_by_rsid)

    result_csv = OUT_DIR / "personal_pharmacogenetic_report_full.csv"
    coverage_csv = OUT_DIR / "variant_coverage_full.csv"

    result_df.to_csv(result_csv, index=False)
    coverage_df = (
        result_df.groupby(["drug", "gene", "actionability"])
        .size()
        .reset_index(name="n")
    )
    coverage_df.to_csv(coverage_csv, index=False)

    html_path = save_html_report(result_df)
    plot_dir = save_plots(result_df)

    print("\nFinal summary")
    print("Interpretations:", len(result_df))
    print("Patients:", result_df["sample"].nunique())
    print("Drugs:", result_df["drug"].nunique())
    print("Drug-gene pairs:", result_df[["drug", "gene"]].drop_duplicates().shape[0])
    print("\nActionability:")
    print(result_df["actionability"].value_counts())

    print("\nSaved outputs")
    print(result_csv)
    print(coverage_csv)
    print(html_path)
    print(plot_dir)
    print(variant_index.head())


if __name__ == "__main__":
    main()
