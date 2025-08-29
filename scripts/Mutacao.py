class Mutacao:
    def __init__(self):
        self.amino = {
            'Ala': 'a', 'Arg': 'r', 'Asn': 'n', 'Asp': 'd', 'Cys': 'c',
            'Gln': 'q', 'Glu': 'e', 'Gly': 'g', 'His': 'h', 'Ile': 'i',
            'Leu': 'l', 'Lys': 'k', 'Met': 'm', 'Phe': 'f', 'Pro': 'p',
            'Ser': 's', 'Thr': 't', 'Trp': 'w', 'Tyr': 'y', 'Val': 'v',
            'Ter': 'z'
        }

    def polimorfismo(self, aminoacidos, id_proteina, pepsalterados, pepsreferencia, indices, idsnps, linhas):
        import re
        
        snp = linhas[1]
        ref = linhas[2]
        pos = int(linhas[3])
        alt = linhas[4]
                
        # Create amino acid sequence with mutation
        aminomutado = aminoacidos[:pos-1] + self.amino[alt] + aminoacidos[pos:]
        pattern = re.compile(r'[^RK]+(R|K|$)')
        start_pos = 0
        
        for match in pattern.finditer(aminomutado):
            pep = match.group(0)
            tam_pep = len(pep)
            
            # Find mutation in peptide
            mut_match = re.search(r'[a-z]', pep)
            if tam_pep >= 7 and mut_match:
                mut = mut_match.group(0)
                pos_sitio = aminomutado.find(pep, start_pos)  # Position of mutated peptide in full sequence
                seqref = aminoacidos[pos_sitio:pos_sitio+tam_pep]  # Reference peptide                
                pos_sitio += 1  # Add 1 to avoid index 0
                pos_mut = pep.find(mut)  # Position of mutation in peptide
                
                # Update data structures
                if pos_sitio not in pepsalterados:
                    pepsalterados[pos_sitio] = []
                    indices[pos_sitio] = []
                    idsnps[pos_sitio] = []
                
                pepsalterados[pos_sitio].append(mut)
                indices[pos_sitio].append(pos_mut)
                idsnps[pos_sitio].append(snp)
                pepsreferencia[pos_sitio] = seqref
                
                start_pos = pos_sitio  # Update start position for next search

    def stopcodon(self, peptideo):
        import re
        
        match = re.search(r'(z)', peptideo)
        if match:
            stop = match.group(1)
            pos = peptideo.find(stop)
            pepstop = peptideo[:pos]
            
            if len(pepstop) >= 7 and len(pepstop) <= 35:
                return pepstop
        return None

    def sitiotriptico(self, peptideo):
        import re
        
        peptriptico = ""
        pattern = re.compile(r'[^RK]+(R|K|$)', re.IGNORECASE)
        
        for match in pattern.finditer(peptideo):
            pep = match.group(0)
            if len(pep) >= 7 and len(pep) <= 35:
                peptriptico += pep
        
        return peptriptico