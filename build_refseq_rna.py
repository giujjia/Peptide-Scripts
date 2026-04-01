import gzip
import logging
import os
import re
import sys
import urllib.parse
import urllib.request
from datetime import datetime
from typing import List

BASE_URL = "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/"
PADRAO_RNA = re.compile(r"human\.\d+\.rna\.fna\.gz")
PASTA_DADOS = "tmp/data_rna"
ARQUIVO_SAIDA = "data/refseqHumanFullNM.fasta"
ARQUIVO_LOG = "data/aggregate_rna.log"

def configurar_logging() -> None:
    os.makedirs("data", exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(ARQUIVO_LOG, mode="w", encoding="utf-8"),
            logging.StreamHandler(sys.stdout),
        ],
    )

def listar_arquivos_rna() -> List[str]:
    logging.info("Consultando índice FTP em %s", BASE_URL)

    with urllib.request.urlopen(BASE_URL) as resposta:
        html = resposta.read().decode("utf-8", errors="replace")

    arquivos = PADRAO_RNA.findall(html)

    arquivos_unicos = []
    vistos = set()

    for arquivo in arquivos:
        if arquivo not in vistos:
            vistos.add(arquivo)
            arquivos_unicos.append(arquivo)

    logging.info("%d arquivo(s) RNA encontrados no índice", len(arquivos_unicos))
    return arquivos_unicos

def obter_last_modified_remoto(url: str) -> datetime:
    requisicao = urllib.request.Request(url, method="HEAD")

    with urllib.request.urlopen(requisicao) as resposta:
        last_modified = resposta.headers.get("Last-Modified")

    if not last_modified:
        return datetime.fromtimestamp(0)

    return datetime.strptime(last_modified, "%a, %d %b %Y %H:%M:%S %Z")

def baixar_se_necessario(nome_arquivo: str) -> str:
    os.makedirs(PASTA_DADOS, exist_ok=True)

    caminho_local = os.path.join(PASTA_DADOS, nome_arquivo)
    url_remota = urllib.parse.urljoin(BASE_URL, nome_arquivo)

    try:
        data_remota = obter_last_modified_remoto(url_remota)
    except Exception as erro:
        logging.warning(
            "Não foi possível obter Last-Modified de %s: %s",
            url_remota,
            erro,
        )
        data_remota = datetime.max

    if os.path.exists(caminho_local):
        data_local = datetime.fromtimestamp(os.path.getmtime(caminho_local))

        if data_local >= data_remota:
            logging.info("Arquivo %s já está atualizado", nome_arquivo)
            return caminho_local

    logging.info("Baixando %s...", nome_arquivo)

    try:
        with urllib.request.urlopen(url_remota) as resposta, open(caminho_local, "wb") as arquivo_saida:
            while True:
                bloco = resposta.read(1024 * 1024)
                if not bloco:
                    break
                arquivo_saida.write(bloco)

        if data_remota != datetime.max:
            timestamp = data_remota.timestamp()
            os.utime(caminho_local, (timestamp, timestamp))

        logging.info("Download concluído de %s", nome_arquivo)

    except Exception as erro:
        logging.error("Erro ao baixar %s: %s", nome_arquivo, erro)

    return caminho_local

def concatenar_fastas(caminhos_arquivos: List[str]) -> None:
    logging.info(
        "Concatenando %d arquivo(s) em %s",
        len(caminhos_arquivos),
        ARQUIVO_SAIDA,
    )

    total_arquivos = 0
    total_sequencias = 0

    with open(ARQUIVO_SAIDA, "w", encoding="utf-8") as saida:
        for caminho in caminhos_arquivos:
            total_arquivos += 1

            try:
                with gzip.open(caminho, "rt", encoding="utf-8", errors="replace") as entrada:
                    for linha in entrada:
                        saida.write(linha)
                        if linha.startswith(">"):
                            total_sequencias += 1

                logging.info("Arquivo %s concatenado", os.path.basename(caminho))

            except Exception as erro:
                logging.error(
                    "Erro ao concatenar arquivo %s: %s",
                    os.path.basename(caminho),
                    erro,
                )

    logging.info(
        "Concatenação concluída. %d arquivo(s), %d sequência(s) FASTA",
        total_arquivos,
        total_sequencias,
    )

def main() -> None:
    configurar_logging()
    logging.info("Iniciando agregação dos arquivos RefSeq RNA humanos")

    try:
        arquivos_rna = listar_arquivos_rna()
    except Exception as erro:
        logging.error("Falha ao obter lista de arquivos: %s", erro)
        return

    caminhos_locais = []

    for nome_arquivo in arquivos_rna:
        caminho = baixar_se_necessario(nome_arquivo)
        caminhos_locais.append(caminho)

    concatenar_fastas(caminhos_locais)

    logging.info("Processo finalizado com sucesso")

if __name__ == "__main__":
    main()
