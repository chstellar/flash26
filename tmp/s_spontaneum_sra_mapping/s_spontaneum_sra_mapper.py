#!/usr/bin/env python3
"""
Extract Supplementary Table 12 accessions from the local PDF and map each
accession to SRA SRX/SRR IDs using NCBI E-utilities.

Designed to avoid HTTP 414 errors:
- esearch is done one accession at a time
- esummary is done with POST in chunks
"""

from __future__ import annotations

import csv
import html
import json
import os
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

from pypdf import PdfReader


DESKTOP = Path(r"D:\OneDrive - Stanford\Desktop")
PDF = DESKTOP / "41588_2026_2641_MOESM1_ESM.pdf"
OUTDIR = DESKTOP / "s_spontaneum_sra_mapping"
OUTDIR.mkdir(exist_ok=True)

BIOPROJECT = "PRJNA1303125"
EMAIL = os.environ.get("NCBI_EMAIL", "")
API_KEY = os.environ.get("NCBI_API_KEY", "")
SLEEP_SECONDS = 0.12 if API_KEY else 0.45


def fetch_text(url: str, data: bytes | None = None, retries: int = 4) -> str:
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, data=data)
            if data is not None:
                req.add_header("Content-Type", "application/x-www-form-urlencoded")
            with urllib.request.urlopen(req, timeout=90) as response:
                return response.read().decode("utf-8")
        except Exception as e:
            if attempt == retries - 1:
                raise
            # NCBI can return 429 when testing repeatedly. Back off gently.
            code = getattr(e, "code", None)
            if code == 429:
                time.sleep(15 * (attempt + 1))
            else:
                time.sleep(2 * (attempt + 1))
    raise RuntimeError("unreachable")


def eutils_params(extra: dict[str, str]) -> dict[str, str]:
    params = {"tool": "s_spontaneum_sra_mapper"}
    if EMAIL:
        params["email"] = EMAIL
    if API_KEY:
        params["api_key"] = API_KEY
    params.update(extra)
    return params


def normalize_query_text(x: str) -> str:
    roman = {
        "Ⅰ": "I",
        "Ⅱ": "II",
        "Ⅲ": "III",
        "Ⅳ": "IV",
        "Ⅴ": "V",
        "Ⅵ": "VI",
    }
    for k, v in roman.items():
        x = x.replace(k, v)
    return x


def extract_table12(pdf_path: Path) -> list[dict[str, str]]:
    reader = PdfReader(str(pdf_path))
    text = "\n".join(page.extract_text() or "" for page in reader.pages)

    # Rows look like:
    # ss_001 yunnan82-17 11.75 39.41 NA
    row_re = re.compile(
        r"\b(ss_\d{3})\s+(.+?)\s+"
        r"((?:NA|[-+]?\d+(?:\.\d+)?))\s+"
        r"((?:NA|[-+]?\d+(?:\.\d+)?))\s+"
        r"((?:NA|[-+]?\d+(?:\.\d+)?))"
        r"(?=\s+ss_\d{3}|\s*$)"
    )
    rows = []
    for sid, accession_name, brix, purity, tiller in row_re.findall(text):
        rows.append(
            {
                "id": sid,
                "accession_name": accession_name.strip(),
                "brix_percent": brix,
                "apparent_purity_percent": purity,
                "tiller_number": tiller,
            }
        )

    # The final row can sit at the very end of the table/page and may be missed
    # depending on PDF extraction whitespace. Add any missing ss rows with a
    # simpler local scan.
    seen = {r["id"] for r in rows}
    loose_re = re.compile(
        r"\b(ss_\d{3})\s+([^\s]+)\s+"
        r"((?:NA|[-+]?\d+(?:\.\d+)?))\s+"
        r"((?:NA|[-+]?\d+(?:\.\d+)?))\s+"
        r"((?:NA|[-+]?\d+(?:\.\d+)?))"
    )
    for sid, accession_name, brix, purity, tiller in loose_re.findall(text):
        if sid not in seen:
            rows.append(
                {
                    "id": sid,
                    "accession_name": accession_name.strip(),
                    "brix_percent": brix,
                    "apparent_purity_percent": purity,
                    "tiller_number": tiller,
                }
            )
            seen.add(sid)

    rows.sort(key=lambda r: int(r["id"].split("_")[1]))
    return rows


def esearch_sra(term: str) -> dict:
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
    params = eutils_params({"db": "sra", "term": term, "retmode": "json", "retmax": "10"})
    text = fetch_text(url + urllib.parse.urlencode(params))
    time.sleep(SLEEP_SECONDS)
    return json.loads(text)


def esearch_project_uids(bioproject: str) -> list[str]:
    """Return SRA UIDs for a BioProject, in one small URL."""
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
    params = eutils_params(
        {
            "db": "sra",
            "term": f"{bioproject}[BioProject]",
            "retmode": "json",
            "retmax": "5000",
        }
    )
    text = fetch_text(url + urllib.parse.urlencode(params), retries=6)
    time.sleep(SLEEP_SECONDS)
    return json.loads(text).get("esearchresult", {}).get("idlist", [])


def choose_uid_for_row(row: dict[str, str]) -> tuple[list[str], str]:
    accession = row["accession_name"]
    sid = row["id"]
    ss_num = sid.split("_")[1]
    queries = [
        accession,
        f'"{accession}"',
        normalize_query_text(accession),
        f"Ss_{ss_num}",
        sid,
    ]

    tried = []
    for query in dict.fromkeys(queries):
        tried.append(query)
        js = esearch_sra(query)
        ids = js.get("esearchresult", {}).get("idlist", [])
        if len(ids) == 1:
            return ids, query
        if len(ids) > 1:
            return ids, query
    return [], " | ".join(tried)


def esummary_sra(uids: list[str], chunk_size: int = 100) -> dict[str, dict]:
    out = {}
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    for i in range(0, len(uids), chunk_size):
        chunk = uids[i : i + chunk_size]
        params = eutils_params({"db": "sra", "id": ",".join(chunk), "retmode": "json"})
        data = urllib.parse.urlencode(params).encode("utf-8")
        js = json.loads(fetch_text(url, data=data))
        result = js.get("result", {})
        for uid in chunk:
            if uid in result:
                out[uid] = result[uid]
        time.sleep(SLEEP_SECONDS)
    return out


def parse_sra_summary(summary: dict) -> dict[str, str]:
    expxml = "<root>" + html.unescape(summary.get("expxml", "")) + "</root>"
    runsxml = "<root>" + html.unescape(summary.get("runs", "")) + "</root>"

    parsed = {
        "title": "",
        "srx": "",
        "srp": "",
        "srs": "",
        "srr": "",
        "prjna": "",
        "samn": "",
        "library_name": "",
        "organism": "",
        "instrument": "",
    }

    try:
        root = ET.fromstring(expxml)
        exp = root.find(".//Experiment")
        study = root.find(".//Study")
        sample = root.find(".//Sample")
        org = root.find(".//Organism")
        lib = root.find(".//LIBRARY_NAME")
        title = root.find(".//Title")
        instrument = root.find(".//Instrument")
        parsed["srx"] = exp.attrib.get("acc", "") if exp is not None else ""
        parsed["srp"] = study.attrib.get("acc", "") if study is not None else ""
        parsed["srs"] = sample.attrib.get("acc", "") if sample is not None else ""
        parsed["organism"] = org.attrib.get("ScientificName", "") if org is not None else ""
        parsed["library_name"] = lib.text or "" if lib is not None else ""
        parsed["title"] = title.text or "" if title is not None else ""
        parsed["instrument"] = next(iter(instrument.attrib.values())) if instrument is not None and instrument.attrib else ""
        bioproject = root.find(".//Bioproject")
        biosample = root.find(".//Biosample")
        parsed["prjna"] = bioproject.text or "" if bioproject is not None else ""
        parsed["samn"] = biosample.text or "" if biosample is not None else ""
    except ET.ParseError:
        pass

    try:
        root = ET.fromstring(runsxml)
        parsed["srr"] = ";".join(run.attrib.get("acc", "") for run in root.findall(".//Run"))
    except ET.ParseError:
        pass
    return parsed


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    rows = extract_table12(PDF)
    write_csv(
        OUTDIR / "s_spontaneum_table12_accessions.csv",
        rows,
        ["id", "accession_name", "brix_percent", "apparent_purity_percent", "tiller_number"],
    )
    print(f"Extracted Table 12 rows: {len(rows)}")

    print(f"Fetching SRA summaries for {BIOPROJECT} ...")
    project_uids = esearch_project_uids(BIOPROJECT)
    print(f"{BIOPROJECT} SRA UID count: {len(project_uids)}")
    summaries = esummary_sra(project_uids)

    by_library = {}
    parsed_by_uid = {}
    for uid, summary in summaries.items():
        parsed = parse_sra_summary(summary)
        parsed_by_uid[uid] = parsed
        lib = parsed.get("library_name", "")
        if lib:
            by_library.setdefault(lib.lower(), []).append(uid)

    search_rows = []
    for i, row in enumerate(rows, start=1):
        ss_num = row["id"].split("_")[1]
        lib_query = f"Ss_{ss_num}"
        ids = by_library.get(lib_query.lower(), [])
        query = f"{BIOPROJECT} library_name={lib_query}"

        # Fallback only if the project/library-name join fails.
        if not ids:
            ids, query = choose_uid_for_row(row)

        chosen = ids[0] if ids else ""
        search_rows.append(
            {
                **row,
                "query_used": query,
                "sra_uid_count": str(len(ids)),
                "sra_uids": ";".join(ids),
                "chosen_uid": chosen,
            }
        )
        print(f"{i:3d}/{len(rows)} {row['id']} {row['accession_name']} -> {len(ids)} project/library uid(s)")

    mapped_rows = []
    for row in search_rows:
        parsed = parsed_by_uid.get(row["chosen_uid"], {})
        if not parsed and row["chosen_uid"]:
            parsed = parse_sra_summary(summaries.get(row["chosen_uid"], {}))
        mapped_rows.append({**row, **parsed})

    fields = [
        "id",
        "accession_name",
        "brix_percent",
        "apparent_purity_percent",
        "tiller_number",
        "query_used",
        "sra_uid_count",
        "sra_uids",
        "chosen_uid",
        "srx",
        "srr",
        "srp",
        "prjna",
        "srs",
        "samn",
        "library_name",
        "organism",
        "instrument",
        "title",
    ]
    write_csv(OUTDIR / "s_spontaneum_290_sra_mapping.csv", mapped_rows, fields)

    with (OUTDIR / "s_spontaneum_290_srr_list.txt").open("w", encoding="utf-8") as fh:
        for r in mapped_rows:
            for srr in r.get("srr", "").split(";"):
                if srr:
                    fh.write(srr + "\n")
    with (OUTDIR / "s_spontaneum_290_srx_list.txt").open("w", encoding="utf-8") as fh:
        for r in mapped_rows:
            if r.get("srx"):
                fh.write(r["srx"] + "\n")

    n_mapped = sum(bool(r.get("srr")) for r in mapped_rows)
    n_ambig = sum(r.get("sra_uid_count") not in ("1", 1) for r in mapped_rows)
    print(f"Mapped rows with SRR: {n_mapped}/{len(mapped_rows)}")
    print(f"Rows with non-unique/no SRA search result: {n_ambig}")
    print(f"Output folder: {OUTDIR}")


if __name__ == "__main__":
    main()
