from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, Tuple


def strip_version(accession: str) -> str:
    return re.sub(r"\.\d+$", "", accession.strip())

def load_nm_cds_positions(path: Path) -> Dict[str, Tuple[str, int, int]]:
    """Carrega posicoes CDS do TSV, retorna nm -> cds_start, cds_end 1-based"""
    result: Dict[str, Tuple[str, int, int]] = {}
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("NM\t"):  # pula cabecalho
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            nm = strip_version(parts[0].strip())
            np = parts[1].strip()
            try:
                cds_start = int(parts[2])
                cds_end = int(parts[3])
            except ValueError:
                continue
            result[nm] = (np, cds_start, cds_end)
    return result

def load_nm_np(path: Path) -> Dict[str, str]:
    """Carrega mapeamento NM para NP do TSV"""
    result: Dict[str, str] = {}
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            nm = parts[0].strip()
            np = parts[1].strip()
            result[nm] = np
    return result

def load_nm_np_with_base(path: Path) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Carrega mapeamento NM para NP com versao exata e sem versao"""
    exact: Dict[str, str] = {}
    base: Dict[str, str] = {}
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            nm = parts[0].strip()
            np = parts[1].strip()
            exact[nm] = np
            base[strip_version(nm)] = np
    return exact, base