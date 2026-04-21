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

RE_DELINS    = re.compile(r"^c\.(\d+)_(\d+)del([ACGTN]+)ins([ACGTN]+)$", re.IGNORECASE)
RE_INS       = re.compile(r"^c\.(\d+)_(\d+)ins([ACGTN]+)$", re.IGNORECASE)
RE_DEL_RANGE = re.compile(r"^c\.(\d+)_(\d+)del([ACGTN]+)?$", re.IGNORECASE)
RE_DEL_SINGLE = re.compile(r"^c\.(\d+)del([ACGTN]+)?$", re.IGNORECASE)
RE_DUP       = re.compile(r"^c\.(\d+)_(\d+)dup([ACGTN]+)$", re.IGNORECASE)

def parse_hgvs(hgvs: str):
    """Interpreta HGVS de indel e retorna cds_pos, del_len, ins_seq ou None"""
    m = RE_DELINS.match(hgvs)
    if m:
        pos1 = int(m.group(1))
        del_seq = m.group(3).upper()
        ins_seq = m.group(4).upper()
        return pos1, len(del_seq), ins_seq

    m = RE_DUP.match(hgvs)
    if m:
        pos2 = int(m.group(2))
        dup_seq = m.group(3).upper()
        return pos2, 0, dup_seq

    m = RE_INS.match(hgvs)
    if m:
        pos1 = int(m.group(1))
        ins_seq = m.group(3).upper()
        return pos1, 0, ins_seq

    m = RE_DEL_RANGE.match(hgvs)
    if m:
        pos1 = int(m.group(1))
        pos2 = int(m.group(2))
        del_seq = (m.group(3) or "").upper()
        del_len = len(del_seq) if del_seq else abs(pos2 - pos1) + 1
        return pos1, del_len, ""

    m = RE_DEL_SINGLE.match(hgvs)
    if m:
        pos1 = int(m.group(1))
        del_seq = (m.group(2) or "").upper()
        del_len = len(del_seq) if del_seq else 1
        return pos1, del_len, ""

    return None

def apply_indel(nucleotides: str, orf: int, cds_pos: int, del_len: int, ins_seq: str) -> str:
    """Aplica indel no transcript, cds_pos 1-based relativo ao inicio do CDS"""
    abs_pos = orf + cds_pos - 1  # 0-based absolute position in transcript
    seqs = list(nucleotides)
    seqs[abs_pos:abs_pos + del_len] = list(ins_seq)
    return "".join(seqs)

def resolve_output_path(category, path_value):
    if os.path.dirname(path_value):
        out_path = path_value
    else:
        out_path = os.path.join("output", category, path_value)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    return out_path

def main():
    if len(sys.argv) not in (4, 5, 6):
        print("Uso: python add_frameshift.py <frameshift_nm.txt> <transcritos> <saida> [saida_inframe] [saida_ext]")
        sys.exit(1)

    frameshift_nm, transcritos, saida = sys.argv[1:4]
    saida_inframe = sys.argv[4] if len(sys.argv) >= 5 else None
    saida_ext = sys.argv[5] if len(sys.argv) == 6 else None
    saida = resolve_output_path("frameshift", saida)
    if saida_inframe:
        saida_inframe = resolve_output_path("frameshift", saida_inframe)
    if saida_ext:
        saida_ext = resolve_output_path("frameshift", saida_ext)

    nm_cds = load_nm_cds_positions(Path("data/nm_cds_positions.tsv"))

    transcripts_dict = {}
    current_nm = None

    try:
        inframe_file = open(saida_inframe, "w") if saida_inframe else None
        ext_file = open(saida_ext, "w") if saida_ext else None
        with open(transcritos, "r") as transcripts_file:
            for line in transcripts_file:
                line = line.strip().replace('\r', '')
                if not line:
                    continue
                if line.startswith('>'):
                    header = line[1:].split()
                    if header:
                        current_nm = strip_version(header[0])
                        if current_nm not in transcripts_dict:
                            transcripts_dict[current_nm] = ""
                else:
                    if current_nm is not None:
                        transcripts_dict[current_nm] += line

        with open(frameshift_nm, "r") as src, open(saida, "w") as dst:
            for raw in src:
                line = raw.strip().replace('\r', '')
                if not line:
                    continue

                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                left = parts[0]
                var_id = parts[1]
                if ":" not in left or not left.startswith("NM_"):
                    continue

                nm_versioned, hgvs = left.split(":", 1)
                nm_id = strip_version(nm_versioned)
                rs_id = f"rs{var_id}"

                parsed = parse_hgvs(hgvs)
                if parsed is None:
                    continue

                cds_pos, del_len, ins_seq = parsed
                is_inframe = (len(ins_seq) - del_len) % 3 == 0

                if is_inframe and not inframe_file:
                    continue

                if nm_id not in transcripts_dict or nm_id not in nm_cds:
                    continue

                np_id, cds_start, cds_end = nm_cds[nm_id]
                orf = cds_start - 1
                ref_protein_len = (cds_end - cds_start) // 3

                if orf <= 0:
                    continue

                sequence = apply_indel(transcripts_dict[nm_id], orf, cds_pos, del_len, ins_seq)

                if not sequence:
                    continue

                protein = ""
                for i in range(orf, len(sequence) - 2, 3):
                    codon = sequence[i:i+3].upper()
                    if codon in CODON_TABLE:
                        if CODON_TABLE[codon] == "STOP":
                            break
                        protein += CODON_TABLE[codon]
                    else:
                        break

                # extensao por frameshift proteina mutada mais longa que a referencia
                is_ext = (not is_inframe) and ext_file and (len(protein) > ref_protein_len)

                if is_inframe:
                    dst_file = inframe_file
                elif is_ext:
                    dst_file = ext_file
                else:
                    dst_file = dst

                if dst_file is None:
                    continue

                dst_file.write(f">{np_id}|{rs_id}\n")
                dst_file.write(protein + "\n")

    except FileNotFoundError as e:
        print(f"Arquivo nao encontrado: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Erro: {e}")
        sys.exit(1)
    finally:
        if inframe_file:
            inframe_file.close()
        if ext_file:
            ext_file.close()

if __name__ == "__main__":
    main()
