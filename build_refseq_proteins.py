import os
import re
import sys
import logging
import urllib.request
import urllib.parse
import gzip
from datetime import datetime
from typing import List

BASE_URL = "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/"
PROTEIN_PATTERN = re.compile(r"human\.\d+\.protein\.faa\.gz")
DATA_DIR = "tmp/data_proteins"
OUTPUT_FILE = "data/refseqHumanFullNP.fasta"
LOG_FILE = "data/aggregate_proteins.log"

def setup_logging() -> None:
    """Configura logging para console e arquivo."""
    os.makedirs("data", exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(LOG_FILE, mode="w", encoding="utf-8"),
            logging.StreamHandler(sys.stdout),
        ],
    )


def fetch_protein_files() -> List[str]:
    logging.info("Consultando índice FTP em %s para arquivos de proteínas", BASE_URL)
    response = urllib.request.urlopen(BASE_URL)
    html = response.read().decode("utf-8", errors="replace")
    files = PROTEIN_PATTERN.findall(html)
    # Remover duplicatas preservando a ordem
    unique = []
    seen = set()
  
    for f in files:
        if f not in seen:
            seen.add(f)
            unique.append(f)
    logging.info("Foram encontrados %d arquivos de proteínas.", len(unique))
    return unique

def get_remote_last_modified(url: str) -> datetime:
    """Retorna a data de modificação do arquivo remoto via requisição HEAD."""
    req = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(req) as response:
        last_mod = response.headers.get("Last-Modified")
        if not last_mod:
            return datetime.fromtimestamp(0)
        return datetime.strptime(last_mod, "%a, %d %b %Y %H:%M:%S %Z")

def download_if_needed(filename: str) -> str:
    os.makedirs(DATA_DIR, exist_ok=True)
    local_path = os.path.join(DATA_DIR, filename)
    remote_url = urllib.parse.urljoin(BASE_URL, filename)
    
    try:
        remote_date = get_remote_last_modified(remote_url)
    except Exception as exc:
        logging.warning("Não foi possível obter Last-Modified de %s: %s", remote_url, exc)
        remote_date = datetime.max
    if os.path.exists(local_path):
        local_mtime = datetime.fromtimestamp(os.path.getmtime(local_path))
        if local_mtime >= remote_date:
            logging.info("Arquivo %s já está atualizado.", filename)
            return local_path
    
    logging.info("Baixando %s...", filename)
    try:
        with urllib.request.urlopen(remote_url) as resp, open(local_path, "wb") as out:
            while True:
                chunk = resp.read(1024 * 1024)
                if not chunk:
                    break
                out.write(chunk)
       
        if remote_date != datetime.max:
            ts = remote_date.timestamp()
            os.utime(local_path, (ts, ts))
        logging.info("Download concluído de %s.", filename)
    except Exception as exc:
        logging.error("Erro ao baixar %s: %s", filename, exc)
    return local_path

def concatenate_proteins(file_paths: List[str]) -> None:
    logging.info(
        "Concatenando sequências de %d arquivo(s) em %s", len(file_paths), OUTPUT_FILE
    )
  
    count_files = 0
    count_sequences = 0
    with open(OUTPUT_FILE, "w", encoding="utf-8") as out_f:
        for gz_path in file_paths:
            count_files += 1
            try:
                with gzip.open(gz_path, "rt", encoding="utf-8", errors="replace") as gz_f:
                    for line in gz_f:
                        out_f.write(line)
                        if line.startswith(">"):
                            count_sequences += 1
                logging.info("Arquivo %s concatenado.", os.path.basename(gz_path))
            except Exception as exc:
                logging.error(
                    "Erro ao concatenar arquivo %s: %s", os.path.basename(gz_path), exc
                )
    logging.info(
        "Concatenação concluída. %d arquivos processados, %d sequências escritas.",
        count_files,
        count_sequences,
    )

def main() -> None:
    setup_logging()
    logging.info("Iniciando agregação de proteínas RefSeq humanas.")
    try:
        protein_files = fetch_protein_files()
    except Exception as exc:
        logging.error("Falha ao obter lista de arquivos: %s", exc)
        return
    # Download de cada arquivo
    local_files = []
    for fname in protein_files:
        local_path = download_if_needed(fname)
        local_files.append(local_path)

    concatenate_proteins(local_files)
    logging.info("Processo concluído com sucesso.")

if __name__ == "__main__":
    main()
