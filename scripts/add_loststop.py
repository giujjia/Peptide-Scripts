import re
import sys
import os
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from utils import load_nm_cds_positions, strip_version

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L",
    "CTA": "L", "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S",
    "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A",
    "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "CAT": "H", "CAC": "H",
    "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "TGT": "C", "TGC": "C",
    "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGT": "S",
    "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G",
    "GGG": "G", "TAA": "STOP", "TAG": "STOP", "TGA": "STOP"
}

STOP_CODONS = {"TAA", "TAG", "TGA"}

RE_CDS_SNP = re.compile(r"^c\.(\d+)([ACGT])>([ACGT])$", re.IGNORECASE)

def parse_snp_hgvs(hgvs: str):
    """Retorna cds_pos, ref, alt ou None"""
    m = RE_CDS_SNP.match(hgvs)
    if not m:
        return None
    return int(m.group(1)), m.group(2).upper(), m.group(3).upper()

def apply_snp(nucleotides: str, abs_pos: int, ref: str, alt: str):
    """Aplica SNP em abs_pos 0-based, retorna sequencia mutada ou None se base nao bate"""
    if abs_pos >= len(nucleotides):
        return None
    if nucleotides[abs_pos].upper() != ref:
        return None
    return nucleotides[:abs_pos] + alt + nucleotides[abs_pos + 1:]

def translate_readthrough(nucleotides: str, orf: int, stop_pos: int, new_aa: str) -> str:
    """Traduz do orf substituindo stop_pos pelo new_aa e continua ate o proximo stop"""
    protein = ""
    for i in range(orf, len(nucleotides) - 2, 3):
        aa_pos = (i - orf) // 3 + 1
        codon = nucleotides[i:i + 3].upper()
        if codon not in CODON_TABLE:
            break
        aa = CODON_TABLE[codon]
        if aa == "STOP":
            if aa_pos == stop_pos:
                protein += new_aa
            else:
                break
        else:
            protein += aa
    return protein

def translate_reference(nucleotides: str, orf: int) -> str:
    """Traduz do orf ate o primeiro stop"""
    protein = ""
    for i in range(orf, len(nucleotides) - 2, 3):
        codon = nucleotides[i:i + 3].upper()
        if codon in CODON_TABLE:
            if CODON_TABLE[codon] == "STOP":
                break
            protein += CODON_TABLE[codon]
        else:
            break
    return protein

def find_tryptic(protein: str, stop_pos: int):
    """Encontra ultimo R/K antes de stop_pos e retorna sitio, peptideo triptico e posicao absoluta"""
    site_pos = None
    limit = min(stop_pos - 1, len(protein))
    for i in range(limit):
        if protein[i] in ("R", "K"):
            site_pos = i

    if site_pos is None:
        return None, None, None

    seq_sub = protein[site_pos:]
    m = re.search(r"([^RK]+(R|K|$))", seq_sub)
    if not m:
        return None, None, None

    pep = m.group(1)
    if not (7 <= len(pep) <= 35):
        return None, None, None

    pep_abs = protein.find(pep)
    return site_pos, pep, pep_abs

def resolve_output_path(category, path_value):
    if os.path.dirname(path_value):
        out_path = path_value
    else:
        out_path = os.path.join("output", category, path_value)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    return out_path

def main():
    if len(sys.argv) != 5:
        print("Uso: python add_loststop.py <stop_loss_nm.txt> <transcritos> <saida> <relacao>")
        sys.exit(1)

    stop_loss_nm, transcritos, saida, relacao = sys.argv[1:5]
    saida = resolve_output_path("stoploss", saida)
    relacao = resolve_output_path("stoploss", relacao)

    nm_cds = load_nm_cds_positions(Path("data/nm_cds_positions.tsv"))

    # carrega transcritos do FASTA canonico
    transcripts_dict = {}
    current_nm = None
    with open(transcritos, "r") as f:
        for line in f:
            line = line.strip().replace("\r", "")
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:].split()
                if header:
                    current_nm = strip_version(header[0])
                    if current_nm not in transcripts_dict:
                        transcripts_dict[current_nm] = ""
            elif current_nm is not None:
                transcripts_dict[current_nm] += line

    try:
        with open(stop_loss_nm, "r") as src, \
             open(saida, "w") as out_f, \
             open(relacao, "w") as rel_f:
            for raw in src:
                line = raw.strip().replace("\r", "")
                if not line:
                    continue

                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                left = parts[0]
                var_id = parts[1].strip()
                if ":" not in left or not left.startswith("NM_"):
                    continue

                nm_versioned, hgvs = left.split(":", 1)
                nm_id = strip_version(nm_versioned)
                rs_id = f"rs{var_id}"

                parsed = parse_snp_hgvs(hgvs)
                if parsed is None:
                    continue
                cds_pos, ref, alt = parsed

                if nm_id not in nm_cds or nm_id not in transcripts_dict:
                    continue

                np_id, cds_start, cds_end = nm_cds[nm_id]
                orf = cds_start - 1

                # posicao 0-based do SNP no transcript
                abs_pos = cds_start + cds_pos - 2

                nucleotides = transcripts_dict[nm_id]

                # aplica SNP, descarta se base referencia nao bate
                mutated = apply_snp(nucleotides, abs_pos, ref, alt)
                if mutated is None:
                    continue

                # verifica se stop codon foi destruido
                stop_codon = mutated[cds_end - 3:cds_end].upper()
                if stop_codon in STOP_CODONS:
                    continue

                new_aa = CODON_TABLE.get(stop_codon, "X")
                if new_aa == "STOP":
                    continue

                # posicao do stop na proteina 1-based
                stop_pos = (cds_end - cds_start + 1) // 3

                protein = translate_readthrough(mutated, orf, stop_pos, new_aa)
                if not protein:
                    continue

                site_pos, pep, pep_abs = find_tryptic(protein, stop_pos)
                if pep is None:
                    continue

                ref_protein = translate_reference(nucleotides, orf)
                ref_sub = ref_protein[site_pos:site_pos + (stop_pos - site_pos) - 1]
                m = re.search(r"([^RK]+(R|K|$))", ref_sub)
                pep_ref = m.group(1) if m else ref_sub

                out_f.write(f"{np_id}\n")
                out_f.write(f"{pep}\n")
                rel_f.write(f"{np_id}\t{rs_id}\t{pep_abs}\n{pep_ref}\t{pep}\n")

    except FileNotFoundError as e:
        print(f"Arquivo nao encontrado: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Erro: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
