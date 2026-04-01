import sys
import os


def resolve_output_path(category, path_value):
    if os.path.dirname(path_value):
        out_path = path_value
    else:
        out_path = os.path.join("output", category, path_value)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    return out_path

if len(sys.argv) < 4:
    print("Uso: python script.py <dbSNP> <transcritos> <saida>")
    sys.exit(1)

arqentrada = sys.argv[1]
arqentrada2 = sys.argv[2]  
arqtmp = sys.argv[3]
arqtmp = resolve_output_path("frameshift", arqtmp)

try:
    entrada = open(arqentrada, 'r')
    entrada2 = open(arqentrada2, 'r')
    saida = open(arqtmp, 'w')
except FileNotFoundError as e:
    print(f"Arquivo nao encontrado: {e}")
    sys.exit(1)

aminoacids = {
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

seq_transcrits = {}
orfs = {}
id_nm = None


for lin in entrada2:
    lin = lin.strip().replace('\r', '')
    
    if lin.startswith('>'):
        # Linha de cabeçalho FASTA
        lin = lin[1:]  # Remove o > do início
        head = lin.split('|')
        
        if len(head) >= 3:  # Verifica se tem campos suficientes
            id_nm = head[0]
            id_np = head[1]
            
            # Armazena informações do ORF
            orfs[id_nm] = f"{id_np}\t{head[2]}"
            
            # Inicializa entrada no dicionário de sequências
            seq_transcrits[id_nm] = ""
    else:
        # Concatena as linhas de sequência
        if id_nm is not None and id_nm in seq_transcrits:
            seq_transcrits[id_nm] += lin

# Fecha o arquivo após leitura
entrada2.close()

# Processa as variações
for lin in entrada:
    lin = lin.strip().replace('\r', '')
    if not lin:  # Pula linhas vazias
        continue
        
    linhas = lin.split('\t')
    
    if len(linhas) < 6:  # Verifica número mínimo de colunas
        continue
    
    chave_nm = linhas[0]
    indel = linhas[1]
    snp = linhas[-1]  # Último elemento
    
    if chave_nm not in orfs or chave_nm not in seq_transcrits:
        continue
        
    vet = orfs[chave_nm].split('\t')
    if len(vet) < 2:
        continue
        
    chave_np = vet[0]
    try:
        orf = int(vet[1])
    except ValueError:
        continue
    
    if orf > 0:
        sequence = None
        
        if indel in ["ins", "dup"]:
            if len(linhas) < 4:
                continue
            insertion = linhas[2]
            try:
                position = int(linhas[3]) + orf
            except ValueError:
                continue
            
            seqs = list(seq_transcrits[chave_nm])
            # Insere a mutação
            if position <= len(seqs):
                seqs[position:position] = list(insertion)
                sequence = "".join(seqs)
            else:
                sequence = seq_transcrits[chave_nm]  # Mantém original se posição inválida
            
        else:
            # Lógica para deleções
            if len(linhas) < 6:
                continue
                
            var = linhas[3]
            position = orf
            insertion = ""
            num_del = 0
            
            try:
                if var == "ins":
                    position += int(linhas[5])
                    insertion = linhas[4]
                else:
                    position += int(linhas[3])
                    
                # Determina número de bases a deletar
                if linhas[2].isdigit():
                    num_del = int(linhas[2])
                else:
                    num_del = len(linhas[2])
                    
                seqs = list(seq_transcrits[chave_nm])
                # Remove e insere (splice)
                start_pos = max(0, position-1)
                end_pos = start_pos + num_del
                
                if start_pos < len(seqs):
                    seqs[start_pos:end_pos] = list(insertion)
                    sequence = "".join(seqs)
                else:
                    sequence = seq_transcrits[chave_nm]
                    
            except (ValueError, IndexError):
                sequence = seq_transcrits[chave_nm]
        
        # Escreve apenas se a sequence foi modificada
        if sequence:
            saida.write(f">{chave_np}|{snp}\n")
            
            # Traduz a sequência mutada
            protein_sequence = ""
            for i in range(orf, len(sequence) - 2, 3):  
                codon = sequence[i:i+3].upper()
                
                if codon in aminoacids:
                    if aminoacids[codon] == "STOP":
                        break
                    protein_sequence += aminoacids[codon]
                    
            saida.write(protein_sequence + "\n")

# Fecha os arquivos
entrada.close()
saida.close()
