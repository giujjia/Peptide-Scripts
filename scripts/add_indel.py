import sys
import re
import os


def resolve_output_path(category, path_value):
    if os.path.dirname(path_value):
        out_path = path_value
    else:
        out_path = os.path.join("output", category, path_value)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    return out_path

amino_hash = {
    "Ala": "a", "Arg": "r", "Asn": "n", "Asp": "d", "Cys": "c", "Gln": "q",
    "Glu": "e", "Gly": "g", "His": "h", "Ile": "i", "Leu": "l", "Lys": "k",
    "Met": "m", "Phe": "f", "Pro": "p", "Ser": "s", "Thr": "t", "Trp": "w",
    "Tyr": "y", "Val": "v", "Ter": "z"
}

def stop_codon(peptideo):
    stop_pos = peptideo.find('z')
    if stop_pos == -1:
        return None

    pepstop = peptideo[:stop_pos]
    if 7 <= len(pepstop) <= 35:
        return pepstop
    else:
        return None

def sitio_triptico(peptideo):
    peptideo_triplico = ""
    for match in re.finditer(r'([^RK]+[RK]?)', peptideo, re.IGNORECASE):
        pep = match.group(0) 
        if 7 <= len(pep) <= 35:
            peptideo_triplico += pep
    return peptideo_triplico

def processar_proteinas(proteinas):
    hash_proteinas = {}
    idnp = None
    
    with open(proteinas, "r") as proteinas:
        for lin in proteinas:
            lin = lin.strip()
           
            if lin.startswith(">"):
                idnp = lin.split('|')[3]
                hash_proteinas[idnp] = ""
            else:
                if idnp is not None:    
                    hash_proteinas[idnp] += lin
         
    return hash_proteinas

def process_mutation(hash_proteinas, mutacao, dbsaida, dbpepmutref, dbfinal):
    prev_id = None
    concat_peptideo = ""
    try:
        with open(mutacao, "r") as DBSNP, \
            open(dbsaida, "w") as DBSAIDA, \
            open(dbpepmutref, "w") as DBRELACAO, \
            open(dbfinal, 'w') as DBFINAL:
        
            for lin in DBSNP:
                lin = lin.strip().replace('\r', '')
                linhas = lin.split('\t')

                if (len(linhas) < 5):
                    print("Linha ignorada: não contém informações suficientes.")
                    continue

                np = linhas[0]
                snp_id = linhas[1]
                ref = linhas[2]
                pos1 = int(linhas[3])
                if (len(linhas) > 4):
                    var = linhas[4]
                    if var == "del":
                        var = "del" 
                    elif var.isdigit():
                        var = int(var)
                    else:
                        continue
                else:
                    var = None
                alt = linhas[5] if len(linhas) > 5 else None
                indel = linhas[6] if len(linhas) > 6 else None
                tipo = linhas[7] if len(linhas) > 7 else None

                if prev_id and prev_id != np:
                    DBFINAL.write(f">{prev_id}\n{concat_peptideo}\n")
                    concat_peptideo = ""

                if np not in hash_proteinas:
                    continue

                aminoacidos = hash_proteinas[np]

                # Caso 1: Deleção simples (ex: NP_001159474.2 rs139292 Asn 15 del)
                if var == "del" and alt is None:
                    aminoacidos = aminoacidos[:pos1-1] + "x" + aminoacidos[pos1:]
                    for match in re.finditer(r"([^RK]+[RK]?)", aminoacidos):
                        pep = match.group(1)

                        if re.search(r"[a-z]+", pep):
                            pep_pos = aminoacidos.index(pep)
                            pepref = pep
                            aa = amino_hash.get(ref, ref[0].upper())
                            pepref = pepref.replace("x", aa)
                            pep = pep.replace("x", "")

                            if pos1 == 1: # n-terminal
                                pep += pep[1:]

                            tam_pep = len(pep)
                            if 7 <= tam_pep <= 35:
                                DBSAIDA.write(f"{np}\n{pep}\n")
                                DBRELACAO.write(f"{np}\t{snp_id}\t{pep_pos}\n{pepref}\t{pep}\n")
                                concat_peptideo += pep
                # Caso 2: Deleção e inserção (ex: NP_115829.1 rs3081552 Tyr 302 del ins Cys Asn)
                elif var == "del" and alt == "ins":
                    pos2 = pos1 
                    referencia = aminoacidos 
                    seq1 = aminoacidos[:pos1 - 1]
                    seq2 = aminoacidos[pos2:]

                    num_del = 1
                    num_ins = 0

                    for index in range(6, len(linhas)):
                        seq1 += amino_hash.get(linhas[index], linhas[index])
                        num_ins += 1
                    
                    aminoacidos = seq1 + seq2

                    for match in re.finditer(r"([^RK]+[RK]?)", aminoacidos):
                        pep = match.group(1)

                        if re.search(r"[a-z]+", pep):
                            if "z" in pep:
                                pep = stop_codon(pep)
                                if not pep:
                                    continue
                                
                            if re.search(r"[rk]", pep):
                                pep = sitio_triptico(pep)

                            if pos1 == 1: # N-terminal
                                pep += pep[1:]

                            # Reportando peptídeo de referência
                            peptideos = list(referencia)
                            peptideos[pos1 - 1] = peptideos[pos1 - 1].lower()
                            referencia = "".join(peptideos)
                            pepref = ""
                            for ref_match in re.finditer(r"([^RK]+[RK]?)", referencia):
                                ref_pep = ref_match.group(1)
                                if re.search(r"[a-z]", ref_pep):
                                    pepref = ref_pep

                            pep_pos = 0
                            tam_pep = len(pep)
                            if pep:
                                if 7 <= tam_pep <= 35:
                                    DBSAIDA.write(f"{np}\n{pep}\n")
                                    for p_match in re.finditer(r"([^RK]+[RK]?)", pep, re.IGNORECASE):
                                        p = p_match.group(1)
                                        if pos1 == 1:
                                            pep_pos = 0
                                        else:
                                            pep_pos = aminoacidos.find(p)
                                        DBRELACAO.write(f"{np}\t{snp_id}\t{pep_pos}\n{pepref}\t{p}\n")
                                    concat_peptideo += pep
                # Caso 3: Apenas inserção (ex: NP_001003891.1 rs361923 Gln 262 263 Ala ins Gln)
                elif indel == "ins":
                    pos2 = var
                    seq1 = aminoacidos[:pos1]  
                    seq2 = aminoacidos[pos2-1:]
                    num_alts = 0

                    for index in range(7, len(linhas)):
                        seq1 += amino_hash.get(linhas[index], linhas[index])
                        num_alts += 1
                    
                    aminoacidos = seq1 + seq2

                    for match in re.finditer(r"([^RK]+[RK]?)", aminoacidos):
                        pep = match.group(1)

                        if re.search(r"[a-z]", pep):
                            pepref = pep  
                            pepref = re.sub(r"[a-z]+", "", pepref)

                            if "z" in pep:
                                pep = stop_codon(pep)
                                if not pep:
                                    continue
                            if re.search(r"[rk]", pep):
                                pep = sitio_triptico(pep)

                            if pos1 == 1: # N-terminal
                                pep += pep[1:]

                            pep_pos = 0
                            tam_pep = len(pep)
                            if pep:
                                if 7 <= tam_pep <= 35:
                                    DBSAIDA.write(f"{np}\n{pep}\n")
                                    for p_match in re.finditer(r"([^RK]+[RK]?)", pep, re.IGNORECASE):
                                        p = p_match.group(1)
                                        if pos1 == 1:
                                            pep_pos = 0
                                        else:
                                            pep_pos = aminoacidos.find(p)

                                        DBRELACAO.write(f"{np}\t{snp_id}\t{pep_pos}\n{pepref}\t{p}\n")
                                    concat_peptideo += pep
                # Caso 4: Deleção de intervalo (ex: NP_001157937.2 rs2005802 Leu 396 401 Arg del)
                elif indel == "del" and tipo is None:
                    pos2 = var
                    referencia = aminoacidos 
                    seq1 = aminoacidos[:pos1 - 1]  # Parte antes da posição inicial de deleção
                    seq2 = aminoacidos[pos2:]  # Parte após a posição final de deleção
                    num_del = (pos2 - pos1) + 1
                    seq1 += "x"
                    aminoacidos = seq1 + seq2

                    for match in re.finditer(r"([^RK]+[RK]?)", aminoacidos):
                        pep = match.group(1)

                        if re.search(r"[a-z]+", pep):
                            pep_pos = aminoacidos.index(pep)  # Posição do peptídeo na sequência
                            pep = pep.replace("x", "")  # Remove o "x" da sequência do peptídeo

                            # Reportando o peptídeo de referência
                            peptideos = list(referencia)  # Converte a sequência de referência em lista
                            peptideos[pos1 - 1] = peptideos[pos1 - 1].lower()  # Marca mutação em minúsculo
                            referencia = "".join(peptideos)

                            pepref = ""
                            for ref_match in re.finditer(r"([^RK]+[RK]?)", referencia):
                                ref = ref_match.group(1)
                                if re.search(r"[a-z]", ref):
                                    pepref = ref

                            if pos1 == 1: # N-terminal
                                pep = pep[1:]

                            tam_pep = len(pep)
                            if 7 <= tam_pep <= 35:
                                DBSAIDA.write(f"{np}\n{pep}\n")
                                DBRELACAO.write(f"{np}\t{snp_id}\t{pep_pos}\n{pepref}\t{pep}\n")
                                concat_peptideo += pep
                # caso 5: Indel, deleção e inserção
                else:
                    pos2 = var
                    referencia = aminoacidos
                    seq1 = aminoacidos[:pos1 - 1]
                    seq2 = aminoacidos[pos2:]
                    num_del = (pos2 - pos1) + 1
                    num_ins = 0

                    for index in range(8, len(linhas)):
                        seq1 += amino_hash.get(linhas[index], linhas[index])  # Concatena aminoácidos do dbSNP
                        num_ins += 1
                    aminoacidos = seq1 + seq2

                    for match in re.finditer(r"([^RK]+[RK]?)", aminoacidos):
                        pep = match.group(1)

                        if re.search(r"[a-z]+", pep):
                          
                            if "z" in pep:
                                pep = stop_codon(pep)
                                if not pep:
                                    continue
                            
                            if re.search(r"[rk]", pep):
                                pep = sitio_triptico(pep)
                                
                            if pos1 == 1: # n-terminal
                                pep += pep[1:]

                            # Reportando o peptídeo referência
                            peptideos = list(referencia)
                            peptideos[pos1 - 1] = peptideos[pos1 - 1].lower()
                            referencia = "".join(peptideos)

                            pepref = ""
                            for ref_match in re.finditer(r"([^RK]+[RK]?)", referencia):
                                ref = ref_match.group(1)
                                if re.search(r"[a-z]", ref):
                                    pepref = ref 

                            pep_pos = 0
                            tam_pep = len(pep)
                            if pep:
                                if 7 <= tam_pep <= 35:
                                    DBSAIDA.write(f"{np}\n{pep}\n")
                                    for p_match in re.finditer(r"([^RK]+[RK]?)", pep, re.IGNORECASE):
                                        p = p_match.group(1)
                                        if pos1 == 1:
                                            pep_pos = 0
                                        else:
                                            pep_pos = aminoacidos.find(p)
                                        DBRELACAO.write(f"{np}\t{snp_id}\t{pep_pos}\n{pepref}\t{p}\n")
                                    concat_peptideo += pep
                prev_id = np          
            if prev_id:
                DBFINAL.write(f">{prev_id}\n{concat_peptideo}\n")
         
    except FileNotFoundError:
        print(f"Erro: O arquivo '{mutacao}' não foi encontrado.")
    except PermissionError:
        print(f"Erro: Permissão negada para acessar um dos arquivos.")
    except ValueError as e:
        print(f"Erro de valor: {e}")
    except FileNotFoundError as e:
        print(f"Erro: {e}")
        
def main():
    if len(sys.argv) != 6:
        print("Uso: python script.py <arquivo_proteinas> <arquivo_mutacao> <arquivo_saida> <arquivo_relacao> <arquivo_final>")
        sys.exit(1)

    proteinas, mutacao, dbsaida, dbpepmutref, dbfinal= sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]
    dbsaida = resolve_output_path("indel", dbsaida)
    dbpepmutref = resolve_output_path("indel", dbpepmutref)
    dbfinal = resolve_output_path("indel", dbfinal)

    hash_proteinas = processar_proteinas(proteinas)
    process_mutation(hash_proteinas, mutacao, dbsaida, dbpepmutref, dbfinal)

if __name__ == "__main__":
    main()