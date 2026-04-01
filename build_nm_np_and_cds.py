import os
import re
import sys
import gzip
import logging
import urllib.request
import urllib.parse
from collections import defaultdict
from datetime import datetime
from typing import List, Tuple

try:
    from Bio import SeqIO  # type: ignore
except ImportError:
    print(
        "Erro biblioteca BioPython nao encontrada instale com pip install biopython",
        file=sys.stderr,
    )
    raise

BASE_URL = "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/"
FILE_PATTERN = re.compile(r"human\.\d+\.rna\.gbff\.gz")
DOWNLOAD_DIR = "tmp/data_cds"
RESULT_FILE = "data/nm_cds_positions.tsv"
NM_NP_FILE = "data/nm.np.txt"
LOG_FILE = "data/cds_extraction.log"

def setup_logging() -> None:
    os.makedirs("data", exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(LOG_FILE, mode="w", encoding="utf-8"),
            logging.StreamHandler(sys.stdout),
        ],
    )

def fetch_directory_listing() -> List[str]:
    logging.info("Consultando indice FTP em %s", BASE_URL)
    response = urllib.request.urlopen(BASE_URL)
    html = response.read().decode("utf-8", errors="replace")
    files = FILE_PATTERN.findall(html)

    seen = set()
    unique_files = []

    for fname in files:
        if fname not in seen:
            seen.add(fname)
            unique_files.append(fname)

    logging.info("%d arquivo(s) .gbff.gz encontrado(s) no indice", len(unique_files))
    return unique_files

def get_remote_last_modified(url: str) -> datetime:
    req = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(req) as response:
        last_mod = response.headers.get("Last-Modified")
        if not last_mod:
            return datetime.fromtimestamp(0)
        return datetime.strptime(last_mod, "%a, %d %b %Y %H:%M:%S %Z")

def download_file(filename: str) -> None:
    os.makedirs(DOWNLOAD_DIR, exist_ok=True)
    local_path = os.path.join(DOWNLOAD_DIR, filename)
    remote_url = urllib.parse.urljoin(BASE_URL, filename)

    try:
        remote_date = get_remote_last_modified(remote_url)
    except Exception as exc:
        logging.warning("Nao foi possivel obter Last-Modified de %s: %s", remote_url, exc)
        remote_date = datetime.max

    if os.path.exists(local_path):
        local_mtime = datetime.fromtimestamp(os.path.getmtime(local_path))
        if local_mtime >= remote_date:
            logging.info("Arquivo %s ja esta atualizado ignorando download", filename)
            return

    logging.info("Baixando %s", filename)

    try:
        with urllib.request.urlopen(remote_url) as response, open(local_path, "wb") as out_f:
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                out_f.write(chunk)

        mtime = remote_date.timestamp() if remote_date != datetime.max else None
        if mtime:
            os.utime(local_path, (mtime, mtime))

        logging.info("Download de %s concluido", filename)
    except Exception as exc:
        logging.error("Falha ao baixar %s: %s", filename, exc)

def parse_gbff(file_path: str) -> List[Tuple[str, str, int, int, str]]:
    results: List[Tuple[str, str, int, int, str]] = []
    file_name = os.path.basename(file_path)
    logging.info("Processando arquivo %s", file_name)
    count = 0
    missing_cds = 0
    missing_np = 0

    try:
        with gzip.open(file_path, "rt", encoding="utf-8", errors="replace") as handle:
            for record in SeqIO.parse(handle, "genbank"):
                count += 1
                acc = record.id

                if not acc.startswith("NM_"):
                    continue

                cds_features = [feat for feat in record.features if feat.type == "CDS"]
                if not cds_features:
                    missing_cds += 1
                    logging.warning(
                        "Transcrito %s no arquivo %s nao possui CDS anotado", acc, file_name
                    )
                    continue

                for cds in cds_features:
                    try:
                        start = int(cds.location.start) + 1
                        end = int(cds.location.end)

                        if start >= end:
                            logging.warning(
                                "CDS invalido em %s no arquivo %s start=%s end=%s",
                                acc,
                                file_name,
                                start,
                                end,
                            )
                            continue

                        if end > len(record.seq):
                            logging.warning(
                                "CDS fora do tamanho do transcrito em %s no arquivo %s end=%s tamanho=%s",
                                acc,
                                file_name,
                                end,
                                len(record.seq),
                            )
                            continue

                        np_id = ""
                        if "protein_id" in cds.qualifiers:
                            np_id = cds.qualifiers["protein_id"][0]
                        else:
                            missing_np += 1
                            logging.warning(
                                "CDS de %s no arquivo %s nao possui protein_id",
                                acc,
                                file_name,
                            )

                        results.append((acc, np_id, start, end, file_name))
                    except Exception as exc:
                        logging.error(
                            "Erro ao extrair CDS de %s no arquivo %s: %s", acc, file_name, exc
                        )
                        continue

        logging.info(
            "Arquivo %s: %d registros processados %d transcritos sem CDS %d CDS sem NP",
            file_name,
            count,
            missing_cds,
            missing_np,
        )
    except Exception as exc:
        logging.error("Erro ao processar %s: %s", file_name, exc)

    return results

def escrever_nm_cds(all_results: List[Tuple[str, str, int, int, str]]) -> None:
    logging.info("Escrevendo resultados para %s", RESULT_FILE)

    try:
        with open(RESULT_FILE, "w", encoding="utf-8") as out_f:
            out_f.write("NM\tNP\tCDS_start\tCDS_end\torigem\n")
            for nm, np_id, start, end, src in all_results:
                out_f.write(f"{nm}\t{np_id}\t{start}\t{end}\t{src}\n")
        logging.info("Foram gravadas %d linhas em %s", len(all_results), RESULT_FILE)
    except Exception as exc:
        logging.error("Erro ao escrever arquivo de resultados: %s", exc)

def escrever_nm_np(all_results: List[Tuple[str, str, int, int, str]]) -> None:
    logging.info("Gerando arquivo de relacao NM NP em %s", NM_NP_FILE)

    pares_nm_np = sorted(
        set((nm, np_id) for nm, np_id, start, end, src in all_results if np_id)
    )

    mapa_nm_para_np = defaultdict(set)
    for nm, np_id in pares_nm_np:
        mapa_nm_para_np[nm].add(np_id)

    nm_com_multiplos_np = {
        nm: sorted(nps)
        for nm, nps in mapa_nm_para_np.items()
        if len(nps) > 1
    }

    if nm_com_multiplos_np:
        logging.warning(
            "Foram encontrados %d NM com mais de um NP associado",
            len(nm_com_multiplos_np),
        )
        for nm, nps in sorted(nm_com_multiplos_np.items()):
            logging.warning("NM %s possui multiplos NP: %s", nm, ", ".join(nps))

    try:
        with open(NM_NP_FILE, "w", encoding="utf-8") as out_f:
            for nm, np_id in pares_nm_np:
                out_f.write(f"{nm}\t{np_id}\n")
        logging.info("Foram gravadas %d linhas em %s", len(pares_nm_np), NM_NP_FILE)
    except Exception as exc:
        logging.error("Erro ao escrever arquivo nm_np: %s", exc)

def main() -> None:
    setup_logging()
    logging.info("Iniciando extracao de CDS de transcritos NM")

    try:
        files_to_download = fetch_directory_listing()
    except Exception as exc:
        logging.error("Falha ao obter lista de arquivos: %s", exc)
        return

    for fname in files_to_download:
        download_file(fname)

    all_results: List[Tuple[str, str, int, int, str]] = []
    gbff_files = [
        os.path.join(DOWNLOAD_DIR, fname)
        for fname in files_to_download
        if os.path.exists(os.path.join(DOWNLOAD_DIR, fname))
    ]

    for gbff in gbff_files:
        all_results.extend(parse_gbff(gbff))

    escrever_nm_cds(all_results)
    escrever_nm_np(all_results)

    logging.info("Processo concluido")

if __name__ == "__main__":
    main()
