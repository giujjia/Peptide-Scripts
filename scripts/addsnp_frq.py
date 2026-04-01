import sys
import re
import os
from collections import defaultdict
from Mutacao import Mutacao


def resolve_output_path(category, path_value):
    if os.path.dirname(path_value):
        out_path = path_value
    else:
        out_path = os.path.join("output", category, path_value)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    return out_path

def process_protein_file(dbproteinas):
    hash_proteinas = {}
    with open(dbproteinas, "r") as proteinas:
        idnp = ""
        for lin in proteinas:
            lin = lin.strip()
            if lin.startswith(">"):
                head = lin.split("|")
                idnp = head[3]
                idnp = re.sub(r'\.\d+$', '', idnp)
                hash_proteinas[idnp] = ""
            else:
                hash_proteinas[idnp] += lin
    return hash_proteinas

#Gera todas as combinações possíveis de mutações e retorna uma lista de combinações em formato binário
def generate_mutation_combinations(total_mutations): 
    combinations = []
    valor = (2 ** total_mutations) - 1
    
    for i in range(1, valor + 1):
        binarios = list(reversed(list(bin(i)[2:])))
        # Completar com zeros à direita se necessário
        while len(binarios) < total_mutations:
            binarios.append('0')
        combinations.append(binarios)
    
    return combinations

#Aplica uma combinação específica de mutações ao peptídeo de referência
#Retorna o peptídeo mutado e a lista de SNPs aplicados
def apply_mutations(peptide_reference, alterations, positions, snp_ids, binary_combination):
    pepmutado = list(peptide_reference)
    snps = []
    
    for pos, idx in enumerate(binary_combination):
        if pos < len(positions):
            posmut = positions[pos]
            if 0 <= posmut < len(pepmutado):
                mut = alterations[pos]
                if idx == '1':
                    pepmutado[posmut] = mut
                    snp = snp_ids[pos]
                    if snp not in snps:
                        snps.append(snp)
    
    return "".join(pepmutado), snps


#Trata casos N-terminal, sítio tríptico, stopcodon
#Retorna o peptídeo processado e um flag indicando se é um stopcodon (para add no dbfinal.fasta)
def handle_special_peptide_cases(pepmutado, key, alterations_length, pos, mutacao):
    is_stopcodon = False
    
    # Verificar condições para peptideos em posição n-terminal
    if key == 1 and alterations_length > 1 and pos >= 1:
        print("cond 1 - n-terminal")
        print(f"chave: {key}\npos:{pos}\ntam: {alterations_length}")
        pepmutado = pepmutado + pepmutado[1:]
        print(pepmutado)
    
    # Verifica se contém r ou k (sítio tríptico)
    if re.search(r'[r|k]', pepmutado): 
        print("cond 2 - sitio triptico")
        print(f"chave: {key}\npos:{pos}\ntam: {alterations_length}")
        try:
            resultado = mutacao.sitiotriptico(pepmutado)
            if resultado:
                pepmutado = resultado
        except Exception as e:
            print(f"Erro ao processar sitiotriptico: {e}")
    
    # Mutacao nonsense (STOPCODON)
    if 'z' in pepmutado:
        print("cond 3 - stopcodon")
        print(f"chave: {key}\npos:{pos}\ntam: {alterations_length}")
        try:
            pepstop = mutacao.stopcodon(pepmutado)
            if pepstop:
                pepmutado = pepstop
                is_stopcodon = True
            else:
                return None, False
        except Exception as e:
            print(f"Erro ao processar stopcodon: {e}")
            return None, False

    return pepmutado, is_stopcodon

def write_dbrelacao_entry(id_anterior, peptide_reference, mutation_peptide, snp_info, pos_sitio, dbrelacao):
    multipep_match = re.search(r"([^R|K]+(R|K|$))+([^R|K]+(R|K|$))", mutation_peptide, re.IGNORECASE)
    if multipep_match:
        pep1 = multipep_match.group(1)
        pep2 = multipep_match.group(3)
        dbrelacao.write(f"{id_anterior}\t{snp_info}\t{pos_sitio}\n{peptide_reference}\t{pep1}\n")
        dbrelacao.write(f"{id_anterior}\t{snp_info}\t{pos_sitio}\n{peptide_reference}\t{pep2}\n")
    else:
        dbrelacao.write(f"{id_anterior}\t{snp_info}\t{pos_sitio}\n{peptide_reference}\t{mutation_peptide}\n")

def process_previous_protein(id_anterior, hash_proteinas, pepsalterados, pepsreferencia, indices, idsnps, dbsaida, dbrelacao, dbfinal, mutacao):
    print(f"\n=== Entrando em process_previous_protein ===")
    print(f"id_anterior: {id_anterior}")
    print(f"hash_proteinas[{id_anterior}]: {hash_proteinas.get(id_anterior, 'N/A')}")
    print(f"indices: {indices}")
    print(f"pepsalterados: {pepsalterados}")
    print(f"pepsreferencia: {pepsreferencia}")
    print(f"idsnps: {idsnps}")

    peptideos_normais = []
    peptideos_stopcodon = []
    
    for chave in sorted(indices.keys()):
        print(chave)

        totalmutacoes = len(indices[chave])
        print(f"totalmutacoes: {totalmutacoes}")

        valor = (2 ** totalmutacoes) - 1
        print(f"valor: {valor}")
        peptideos = []
        snp_pep = {}

        peptideo = pepsreferencia.get(chave, "")
        if len(peptideo) <= 7 or len(peptideo) >= 35:
            continue

        # Gerar todas as combinações possíveis de mutações
        combinacoes = generate_mutation_combinations(totalmutacoes)
        
        for i, binarios in enumerate(combinacoes, 1):
            print(f"\nProcessando combinação {i} de {valor}")
            print(f"binários: {binarios}")
            
            # Aplicar as mutações ao peptídeo
            pepmutado, snps_aplicados = apply_mutations(
                pepsreferencia[chave], 
                pepsalterados[chave], 
                indices[chave], 
                idsnps[chave], 
                binarios
            )
            
            print(f"pepmutado após aplicação de mutações: {pepmutado}")

            # Tratar casos especiais (N-terminal, sítio tríptico, stopcodon)
            pepmutado, is_stopcodon = handle_special_peptide_cases(
                pepmutado, 
                chave, 
                len(pepsalterados[chave]), 
                len(binarios), 
                mutacao
            )
            
            if pepmutado is None:
                continue

            if pepmutado and len(pepmutado) >= 7 and pepmutado not in peptideos:
                print("cond 4 - 1° ocorrencia da mutação")
                print(f"chave: {chave}\npos:{len(binarios)}\ntam: {len(pepsalterados[chave])}")
                if pepmutado not in snp_pep:
                    snp_pep[pepmutado] = []
                snp_pep[pepmutado].append("\t".join(snps_aplicados))
                peptideos.append(pepmutado)
                
                # Adicionar à lista apropriada para o arquivo final
                if is_stopcodon:
                    peptideos_stopcodon.append(pepmutado)
                else:
                    peptideos_normais.append(pepmutado)

        print(f"\nPeptídeos para {id_anterior}: {peptideos}")
        dbsaida.write(f">{id_anterior}\n")
        for pep in peptideos:
            peptideoreferencia = pepsreferencia[chave]
            dbsaida.write(f"{pep}\n")
            print(f"{pep}\n")
            
            for match in re.finditer(r"([^R|K]+(R|K|$))", peptideoreferencia, re.IGNORECASE):
                pepref = match.group(1)
                pos_sitio = hash_proteinas[id_anterior].find(pepref)
                
                if len(pepref) >= 7:
                    for snp_info in snp_pep[pep]:
                        write_dbrelacao_entry(id_anterior, pepref, pep, snp_info, pos_sitio, dbrelacao)
    
    # Escrever para o arquivo dbfinal.fasta
    if peptideos_normais:
        dbfinal.write(f">{id_anterior}\n")
        dbfinal.write(f"{''.join(peptideos_normais)}\n")
    
    # Adicionar peptídeos stopcodon separadamente
    if peptideos_stopcodon:
        for pep_stop in peptideos_stopcodon:
            dbfinal.write(f">{id_anterior}\n")
            dbfinal.write(f"{pep_stop}\n")

def process_mutations(dbmutacao, hash_proteinas, mutacao, dbsaida, dbrelacao, dbfinal):
    id_anterior = ""
    pepsalterados = defaultdict(list)
    pepsreferencia = {}
    indices = defaultdict(list)
    idsnps = defaultdict(list)

    with open(dbmutacao, "r") as dbsnp:
        for lin in dbsnp:
            lin = lin.strip()
            linhas = lin.split("\t")
            id = re.sub(r'\.\d+$', '', linhas[0])
            
            if id in hash_proteinas:
                print(f"\n--- Processando proteína {id} ---", file=sys.stderr)
                aminoacidos = hash_proteinas[id]

                if id_anterior == "" or id_anterior == id:
                    id_anterior = id
                    mutacao.polimorfismo(aminoacidos, id, pepsalterados, pepsreferencia, indices, idsnps, linhas)
                    print(f"DEBUG IF: Chamada polimorfismo -> ID: {id}, Linhas: {linhas}", file=sys.stderr)
                else:
                    process_previous_protein(id_anterior, hash_proteinas, pepsalterados, pepsreferencia, indices, idsnps, dbsaida, dbrelacao, dbfinal, mutacao)
                    id_anterior = id
                    pepsalterados = defaultdict(list)
                    pepsreferencia = {}
                    indices = defaultdict(list)
                    idsnps = defaultdict(list)
                    mutacao.polimorfismo(aminoacidos, id, pepsalterados, pepsreferencia, indices, idsnps, linhas)
                    print(f"DEBUG ELSE: Chamada polimorfismo -> ID: {id}, Linhas: {linhas}", file=sys.stderr)

        if indices:
            print("IF INDICES")
            process_previous_protein(id_anterior, hash_proteinas, pepsalterados, pepsreferencia, indices, idsnps, dbsaida, dbrelacao, dbfinal, mutacao)

def main():
    if len(sys.argv) != 6:
        print("Uso: python script.py dbproteinas dbmutacao dbsaida dbpepmutref dbfinal")
        sys.exit(1)

    dbproteinas, dbmutacao, dbsaida, dbpepmutref, dbfinal = sys.argv[1:6]
    dbsaida = resolve_output_path("snp", dbsaida)
    dbpepmutref = resolve_output_path("snp", dbpepmutref)
    dbfinal = resolve_output_path("snp", dbfinal)

    with open(dbsaida, "w") as dbsaida, \
         open(dbpepmutref, "w") as dbrelacao, \
         open(dbfinal, "w") as dbfinal:

        mutacao = Mutacao()
        hash_proteinas = process_protein_file(dbproteinas)
        process_mutations(dbmutacao, hash_proteinas, mutacao, dbsaida, dbrelacao, dbfinal)

if __name__ == "__main__":
    main()
