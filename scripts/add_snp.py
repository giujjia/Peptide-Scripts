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

if len(sys.argv) != 5:
    print("Uso: python script.py <proteinas> <mutacao> <dbsaida> <dbpepmutref>")
    sys.exit(1)

proteinas, mutacao, dbsaida, dbpepmutref = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
dbsaida = resolve_output_path("snp", dbsaida)
dbpepmutref = resolve_output_path("snp", dbpepmutref)
dbfinal = resolve_output_path("snp", "dbf.txt")

# Definindo o dicionário de aminoácidos
amino = {
    "Ala": "a", "Arg": "r", "Asn": "n", "Asp": "d", "Cys": "c", "Gln": "q", "Glu": "e", "Gly": "g", "His": "h",
    "Ile": "i", "Leu": "l", "Lys": "k", "Met": "m", "Phe": "f", "Pro": "p", "Ser": "s", "Thr": "t", "Trp": "w",
    "Tyr": "y", "Val": "v", "Ter": "z"
}

idnp = None
hash_proteinas = {}
concat_peptideo = ""  
prev_id = None 

try:
    with open(proteinas, 'r') as PROTEINAS, \
         open(mutacao, 'r') as DBSNP, \
         open(dbsaida, 'w') as DBSAIDA, \
         open(dbpepmutref, 'w') as DBRELACAO, \
         open(dbfinal, 'w') as DBFINAL: 

        # Processando o arquivo de proteínas
        for lin in PROTEINAS:
            lin = lin.strip().replace('\r', '')
            if lin.startswith(">"):
                head = lin.split("|")
                idnp = head[3].split('.')[0]
                hash_proteinas[idnp] = ""
            else:
                if idnp is not None:
                    hash_proteinas[idnp] += lin

        # Processando o arquivo de mutações
        for lin in DBSNP:
            lin = lin.strip().replace('\r', '')
            linhas = lin.split('\t')

            # Verifica se a linha tem o número esperado de colunas
            if len(linhas) < 5:
                continue

            id_ = linhas[0].split('.')[0]
            snp = linhas[1]
            ref = linhas[2]
            pos = int(linhas[3])
            alt = linhas[4]
            mutacao = f"p.{ref}{pos}{alt}"

            if id_ in hash_proteinas:
                aminoacidos = hash_proteinas[id_]
                aminoacidos = aminoacidos[:pos - 1] + amino[alt] + aminoacidos[pos:]

                pattern = re.compile(r'([^RK]+(R|K|$))')
                for match in pattern.finditer(aminoacidos):
                    pepmutado = match.group(1)
                    tam_pep = len(pepmutado)
                    if 7 <= tam_pep <= 35 and re.search(r'[a-z]', pepmutado):
                        aminoref = amino[ref]
                        aminomut = amino[alt]
                        pepref = pepmutado.replace(aminomut, aminoref)

                        sitiopos = aminoacidos.find(pepmutado)

                        if sitiopos == 0:
                            pepmutado = pepmutado + pepmutado[1:]
                        
                        if re.search(r'[r|k]', pepmutado):  # Verifica se um novo peptídeo tríptico foi criado
                            peptriptico = ""
                            pattern = re.compile(r'([^RK]+(R|K|$))', re.IGNORECASE)  # Busca novamente fragmentos trípticos
                            for match in pattern.finditer(pepmutado):
                                pep = match.group(1)
                                if len(pep) >= 7:
                                    peptriptico += pep  # Concatena fragmentos com tamanho adequado
                            pepmutado = peptriptico 
                        
                        if 'z' in pepmutado:
                            stop = 'z'
                            pos_stop = pepmutado.index(stop)
                            pepstop = pepmutado[:pos_stop]
                            if len(pepstop) >= 7:
                                pepmutado = pepstop
                        
                        if pepmutado: 
                            DBSAIDA.write(f">{id_}\n{pepmutado}\n")


                            if prev_id and prev_id != id_: #verificação para caso o id mudar
                                DBFINAL.write(f"{prev_id}\n{concat_peptideo}\n")
                                concat_peptideo = ""

                            prev_id = id_
                            concat_peptideo += pepmutado #concatena o pep

                            pattern = re.compile(r'([^RK]+(R|K|$))', re.IGNORECASE)
                            for match in pattern.finditer(pepmutado):
                                pepmut = match.group(1)
                                sitiopos = aminoacidos.find(pepmut)
                                DBRELACAO.write(f">{id_}\t{snp}\t{sitiopos}\t{mutacao}\t{pepref}\t{pepmut}\n")

                            if sitiopos == 0 and alt == "Ter":
                                sitiopos = 1  
                                pepmutado = pepmutado[1:]
                                DBRELACAO.write(f">{id_}\t{snp}\t{sitiopos}\t{mutacao}\t{pepref}\t{pepmut}\n")

       #grava o último ID e sequência concatenada
        if prev_id:
            DBFINAL.write(f">{prev_id}\n{concat_peptideo}\n")

except FileNotFoundError as e:
    print(f"Erro: {e}")
