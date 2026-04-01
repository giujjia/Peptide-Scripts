import sys
import re
import os

# Dicionário de tradução de códons para aminoácidos
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

def orf(nucleotideos):
    """
    Localiza o ORF (Open Reading Frame) e traduz a sequência de nucleotídeos para proteína.
    Retorna (orf, proteina) onde orf é a posição de início do ORF e proteina é a sequência traduzida.
    """
    # Procura pelo códon de início ATG
    start_codon = "ATG"
    orf_pos = -1
    
    # Procura o primeiro ATG na sequência
    for i in range(len(nucleotideos) - 2):
        codon = nucleotideos[i:i+3].upper()
        if codon == start_codon:
            orf_pos = i
            break
    
    # Se não encontrou ATG, retorna None
    if orf_pos == -1:
        return (None, None)
    
    # Traduz a sequência a partir do ORF
    proteina = ""
    for i in range(orf_pos, len(nucleotideos) - 2, 3):
        codon = nucleotideos[i:i+3].upper()
        if codon in aminoacids:
            if aminoacids[codon] == "STOP":
                break
            proteina += aminoacids[codon]
        else:
            # Se o códon não está no dicionário, para a tradução
            break
    
    return (orf_pos, proteina)


def resolve_output_path(category, path_value):
    if os.path.dirname(path_value):
        out_path = path_value
    else:
        out_path = os.path.join("output", category, path_value)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    return out_path

def main():
    if len(sys.argv) != 6:
        print("Uso: python add_loststop.py <transcritos> <dbsnp_nm> <dbsnp_np> <dbsaida> <dbpepmutref>")
        sys.exit(1)
    
    transcritos, dbsnp_nm, dbsnp_np, dbsaida, dbpepmutref = sys.argv[1:6]
    dbsaida = resolve_output_path("stoploss", dbsaida)
    dbpepmutref = resolve_output_path("stoploss", dbpepmutref)
    
    try:
        # Abre os arquivos
        with open(transcritos, 'r') as TRANSCRITOS, \
             open(dbsnp_np, 'r') as DBSNP_NP, \
             open(dbsnp_nm, 'r') as DBSNP_NM, \
             open(dbsaida, 'w') as DBSAIDA, \
             open(dbpepmutref, 'w') as DBRELACAO:
            
            # Dicionário para armazenar sequências de transcritos
            hash_transcritos = {}
            idnm = None
            orfs = {}
            
            # Processa o arquivo de transcritos em fasta
            # O ORF será encontrado automaticamente procurando pelo primeiro códon ATG
            for linha in TRANSCRITOS:
                linha = linha.strip()
                
                # Verifica se é linha de cabeçalho FASTA (começa com >)
                if linha.startswith('>'):
                    # Formato FASTA: extrai NM_ do cabeçalho
                    # Remove o > e pega o primeiro campo (ID)
                    linha_sem_gt = linha[1:].strip()
                    campos = linha_sem_gt.split()
                    if campos:
                        idnm = campos[0]
                        # Remove extensão após o ponto
                        idnm = re.sub(r'\..+', '', idnm)
                        # ORF será encontrado automaticamente pela função orf()
                        orfs[idnm] = None
                        # Inicializa a sequência
                        if idnm not in hash_transcritos:
                            hash_transcritos[idnm] = ""
                else:
                    # Linha de sequência (nucleotídeos)
                    if linha != "" and idnm is not None:
                        if idnm not in hash_transcritos:
                            hash_transcritos[idnm] = ""
                        hash_transcritos[idnm] += linha
            
            # Dicionário para armazenar proteínas mutadas
            proteinas = {}
            
            # Processa o arquivo DBSNP_NM
            for linha in DBSNP_NM:
                linha = linha.strip()
                campos = linha.split('\t')
                
                if len(campos) < 4:
                    continue
                
                np = campos[0]
                nm = campos[1]
                # Remove extensão após o ponto
                nm = re.sub(r'\..+', '', nm)
                info = campos[2]  # Contém: c.427_429delTAAinsCAT
                snp = "rs" + campos[3]
                
                # Remove "c." do início
                info = re.sub(r'^c\.', '', info)
                
                if nm in hash_transcritos:
                    nucleotideos = hash_transcritos[nm]
                    orf_val = orfs.get(nm)
                    proteina = None
                    
                    # Se ORF não está definido ou é inválido, busca o ORF
                    if orf_val is None or orf_val < 1:
                        orf_val, proteina = orf(nucleotideos)
                        if orf_val is None:
                            continue
                    else:
                        # Se ORF está definido, traduz a partir dele
                        proteina = ""
                        for i in range(orf_val, len(nucleotideos) - 2, 3):
                            codon = nucleotideos[i:i+3].upper()
                            if codon in aminoacids:
                                if aminoacids[codon] == "STOP":
                                    break
                                proteina += aminoacids[codon]
                            else:
                                break
                    
                    if orf_val is None:
                        continue
                    
                    # Processa diferentes tipos de mutações
                    if 'del' in info and 'ins' in info:  # del e ins
                        mutacao = info
                        # Remove números, underscore e "del"
                        mutacao = re.sub(r'\d+|_|del', '', mutacao)
                        # Divide por letras minúsculas (del, ins, etc)
                        partes = re.split(r'[a-z]+', mutacao)
                        del_seq = partes[0] if len(partes) > 0 and partes[0] else ""
                        ins_seq = partes[1] if len(partes) > 1 and partes[1] else ""
                        
                        # Extrai posições
                        pos_match = re.match(r'(\d+)_(\d+)', info)
                        if pos_match:
                            pos1 = int(pos_match.group(1))
                            pos2_str = pos_match.group(2)
                            # Remove não-dígitos de pos2
                            pos2 = int(re.sub(r'\D+', '', pos2_str))
                        else:
                            continue
                        
                        pos1 += orf_val
                        dels = len(del_seq)
                        
                        # Aplica a mutação (substitui del por ins)
                        nucleotideos_list = list(nucleotideos)
                        if pos1 - 1 < len(nucleotideos_list):
                            # Remove del_seq e insere ins_seq
                            nucleotideos_list[pos1-1:pos1-1+dels] = list(ins_seq)
                            nucleotideos = ''.join(nucleotideos_list)
                        
                        # Traduz novamente
                        orf_val, proteina = orf(nucleotideos)
                        if proteina:
                            if np not in proteinas:
                                proteinas[np] = {}
                            if snp not in proteinas[np]:
                                proteinas[np][snp] = []
                            proteinas[np][snp].append(proteina)
                    
                    elif 'del' not in info and 'ins' in info:  # Apenas ins
                        # Extrai a sequência de inserção (remove números, underscore e letras minúsculas)
                        ins = info
                        ins = re.sub(r'\d+|_|[a-z]+', '', ins)
                        
                        # Extrai posições
                        pos_match = re.match(r'(\d+)_(\d+)', info)
                        if pos_match:
                            pos1 = int(pos_match.group(1))
                            pos2_str = pos_match.group(2)
                            # Remove não-dígitos de pos2
                            pos2 = int(re.sub(r'\D+', '', pos2_str))
                        else:
                            continue
                        
                        pos1 += orf_val
                        
                        # Aplica a inserção (insere na posição pos1-1)
                        nucleotideos_list = list(nucleotideos)
                        if pos1 - 1 <= len(nucleotideos_list):
                            # Insere a sequência completa na posição
                            for i, base in enumerate(ins):
                                nucleotideos_list.insert(pos1 - 1 + i, base)
                            nucleotideos = ''.join(nucleotideos_list)
                        
                        # Traduz novamente
                        orf_val, proteina = orf(nucleotideos)
                        if proteina:
                            if np not in proteinas:
                                proteinas[np] = {}
                            if snp not in proteinas[np]:
                                proteinas[np][snp] = []
                            proteinas[np][snp].append(proteina)
                    
                    elif 'del' in info and 'ins' not in info:  # Apenas del
                        # Extrai a sequência de deleção (remove números, underscore e letras minúsculas)
                        mutacao = info
                        mutacao = re.sub(r'\d+|_|[a-z]', '', mutacao)
                        
                        # Extrai posições
                        pos_match = re.match(r'(\d+)_(\d+)', info)
                        if pos_match:
                            pos1 = int(pos_match.group(1))
                            pos2_str = pos_match.group(2)
                            # Remove não-dígitos de pos2
                            pos2 = int(re.sub(r'\D+', '', pos2_str))
                        else:
                            continue
                        
                        pos1 += orf_val
                        dels = len(mutacao)
                        
                        if dels > 1:
                            # Aplica a deleção
                            nucleotideos_list = list(nucleotideos)
                            if pos1 - 1 < len(nucleotideos_list):
                                # Remove dels caracteres
                                del nucleotideos_list[pos1-1:pos1-1+dels]
                                nucleotideos = ''.join(nucleotideos_list)
                            
                            # Traduz novamente
                            orf_val, proteina = orf(nucleotideos)
                            if proteina:
                                if np not in proteinas:
                                    proteinas[np] = {}
                                if snp not in proteinas[np]:
                                    proteinas[np][snp] = []
                                proteinas[np][snp].append(proteina)
                    
                    else:  # SNP
                        mutacao = info
                        # Extrai posição
                        pos_str = mutacao
                        pos_str = re.sub(r'\D+', '', pos_str)
                        if not pos_str:
                            continue
                        pos = int(pos_str)
                        
                        # Remove dígitos para obter ref>alt
                        mutacao = re.sub(r'\d+', '', mutacao)
                        # Divide por >
                        partes = mutacao.split('>')
                        if len(partes) != 2:
                            continue
                        ref = partes[0]
                        alt = partes[1]
                        
                        pos += orf_val
                        
                        if pos <= len(nucleotideos):
                            # Aplica a mutação
                            nucleotideos_list = list(nucleotideos)
                            if pos - 1 < len(nucleotideos_list):
                                nucleotideos_list[pos-1] = alt
                                nucleotideos = ''.join(nucleotideos_list)
                            
                            # Traduz novamente
                            orf_val, proteina = orf(nucleotideos)
                            if proteina:
                                if np not in proteinas:
                                    proteinas[np] = {}
                                if snp not in proteinas[np]:
                                    proteinas[np][snp] = []
                                proteinas[np][snp].append(proteina)
            
            for linha in DBSNP_NP:
                linha = linha.strip()
                campos = linha.split('\t')
                
                if len(campos) < 4:
                    continue
                
                np = campos[0]
                snp = campos[1]
                try:
                    posicao = int(campos[3])
                except (ValueError, IndexError):
                    continue
                sitiopos = None
                
                if np in proteinas and snp in proteinas[np]:
                    for sequencia in proteinas[np][snp]:
                        aminoacidos = sequencia
                        aminoacidos_list = list(sequencia)
                        sitiopos = None
                        
                        # Encontrar o último R ou K da sequência antes da posição
                        limite = min(posicao - 1, len(aminoacidos_list))
                        for i in range(limite):
                            if aminoacidos_list[i] == "R" or aminoacidos_list[i] == "K":
                                sitiopos = i
                        
                        if sitiopos is not None:
                            # Calcula pepref da posição sitiopos até (posicao - sitiopos) - 1
                            pepref = sequencia[sitiopos:sitiopos + (posicao - sitiopos) - 1]
                            sequencia_sub = sequencia[sitiopos:]
                            
                            # Aplica regex para encontrar peptídeo de referência
                            pepref_match = re.search(r'([^RK]+(R|K|$))', pepref)
                            if pepref_match:
                                pepref = pepref_match.group(1)
                            
                            # Aplica regex para encontrar peptídeos trípticos
                            pattern = re.compile(r'([^RK]+(R|K|$))')
                            # Simula o comportamento /gc do Perl usando finditer
                            matches = list(pattern.finditer(sequencia_sub))
                            if matches:
                                # Pega o primeiro match
                                match = matches[0]
                                pep = match.group(1)
                                tam_pep = len(pep)
                                sitiopos_final = aminoacidos.find(pep)
                                
                                if tam_pep >= 7 and tam_pep <= 35:
                                    DBSAIDA.write(f"{np}\n")
                                    DBSAIDA.write(f"{pep}\n")
                                    DBRELACAO.write(f"{np}\t{snp}\t{sitiopos_final}\n{pepref}\t{pep}\n")
    
    except FileNotFoundError as e:
        print(f"Arquivo nao encontrado: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Erro: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()