""" Exemple basique, juste sur l'ADN aléatoire """
""" ------------------------------------------ """


import random

#  Dico des mots
def construire_dico_mots(seq, w):
    dico_mots = {}
    nb_mots = len(seq) - w + 1
    for i in range(nb_mots):
        mot = seq[i:i+w]
        if mot not in dico_mots:
            dico_mots[mot] = []
        dico_mots[mot].append(i)
    return dico_mots


def find_hits(dico_mots, database, w):
    hits = []
    for seq_id, seq in enumerate(database):
        nb_mots_db = len(seq) - w + 1
        # Parcours des mots de la séquence de la database
        for j in range(nb_mots_db):
            mot_db = seq[j:j+w]
            if mot_db in dico_mots:
                for pos in dico_mots[mot_db]:
                    hits.append((seq_id, pos, j))  # (id_sequence, position_seq1, position_seq2)
    return hits


def extend_hit(hit, sequence, database, w, matrice_substitution, seuil):
    data_id, pos_seq, pos_data = hit
    score = sum(matrice_substitution[sequence[pos_seq + i]][database[data_id][pos_data + i]] for i in range(w))
    max_score = score

    qL, qR = pos_seq, pos_seq + w - 1
    sL, sR = pos_data, pos_data + w - 1
    
    # Extension à droite
    pos_seq = qR
    pos_data = sR

    while max_score - score <= seuil and pos_seq + 1 < len(sequence) and pos_data + 1 < len(database[data_id]):
        pos_seq += 1
        pos_data += 1
        score += matrice_substitution[sequence[pos_seq]][database[data_id][pos_data]]
        if score > max_score:
            max_score = score
            qR = pos_seq
            sR = pos_data
    
    # Extension à gauche
    pos_seq = qL
    pos_data = sL
    score = max_score
    
    while max_score - score <= seuil and pos_seq - 1 >= 0 and pos_data - 1 >= 0:
        pos_seq -= 1
        pos_data -= 1
        score += matrice_substitution[sequence[pos_seq]][database[data_id][pos_data]]
        if score > max_score:
            max_score = score
            qL = pos_seq
            sL = pos_data

    return (qL, sL, qR, sR, max_score)


if __name__ == "__main__":
    
    # Paramètres

    w = 11  # taille des mots
    n = 500 # taille de la séquence
    nb_mots = n - w + 1
    seuil = 50 # seuil de score pour garder un alignement
    taille_database = 1000 # taille de la base de données
    database = [''.join(random.choices('ACGT', k=random.randint(1000, 5000))) for _ in range(taille_database)]

    sequence = ''.join(random.choices('ACGT', k=n))

    matrice_substitution = {
        'A': {'A': 5, 'C': -4, 'G': -4, 'T': -4},
        'C': {'A': -4, 'C': 5, 'G': -4, 'T': -4},
        'G': {'A': -4, 'C': -4, 'G': 5, 'T': -4},
        'T': {'A': -4, 'C': -4, 'G': -4, 'T': 5},
    }




    dico_mots = construire_dico_mots(sequence, w)
    hits = find_hits(dico_mots, database, w)
    
    print(f"Nombre de hits trouvés : {len(hits)}")
    
    alignements = {}
    
    for hit in hits:
        alignment = extend_hit(hit, sequence, database, w, matrice_substitution, seuil)
        qL, sL, qR, sR, score = alignment
        if score not in alignements:
            alignements[score] = []
        alignements[score].append((qL, sL, qR, sR, hit[0]))

        if score >= seuil:
            print(f"Alignement trouvé dans la séquence {hit[0]} avec un score de {score}")
            print(f"Sequence 1: {sequence[qL:qR+1]} (positions {qL} à {qR})")
            print(f"Sequence 2: {database[hit[0]][sL:sR+1]} (positions {sL} à {sR})")
            print()


    alignements = dict(sorted(alignements.items(), reverse=True))
    print(alignements)