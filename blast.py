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
    hits = set()
    for seq_id, seq in enumerate(database):
        nb_mots_db = len(seq) - w + 1
        # Parcours des mots de la séquence de la database
        for j in range(nb_mots_db):
            mot_db = seq[j:j+w]
            if mot_db in dico_mots:
                for pos in dico_mots[mot_db]:
                    hits.add((seq_id, pos, j))  # (id_sequence, position_seq1, position_seq2)
    return hits


def extend_hit(hit, sequence, database, w, M, x_drop=20):
    data_id, pos_q0, pos_s0 = hit

    # score du seed (w-mer)
    seed_score = sum(M[sequence[pos_q0+i]][database[data_id][pos_s0+i]] for i in range(w))

    # bornes initiales = seed
    qL = pos_q0
    sL = pos_s0
    qR = pos_q0 + w - 1
    sR = pos_s0 + w - 1

    # --------- Extension à droite (X-drop local) ----------
    score = seed_score
    best_r = seed_score
    best_qR, best_sR = qR, sR

    q = qR; s = sR
    while q + 1 < len(sequence) and s + 1 < len(database[data_id]):
        q += 1; s += 1
        score += M[sequence[q]][database[data_id][s]]

        if score > best_r:
            best_r = score
            best_qR, best_sR = q, s

        if best_r - score > x_drop:
            break

    # --------- Extension à gauche (X-drop local) ----------
    score = seed_score
    best_l = seed_score
    best_qL, best_sL = qL, sL

    q = qL; s = sL
    while q - 1 >= 0 and s - 1 >= 0:
        q -= 1; s -= 1
        score += M[sequence[q]][database[data_id][s]]

        if score > best_l:
            best_l = score
            best_qL, best_sL = q, s

        if best_l - score > x_drop:
            break

    # --------- Combinaison des gains ----------
    # total_best = seed + gain_droite + gain_gauche
    total_best = seed_score + (best_r - seed_score) + (best_l - seed_score)

    return (best_qL, best_sL, best_qR, best_sR, total_best)


def run_adn_ungapped(w, n, x_drop, seuil, taille_database):
    nb_mots = n - w + 1
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
        alignment = extend_hit(hit, sequence, database, w, matrice_substitution, x_drop)
        qL, sL, qR, sR, score = alignment
        if score >= seuil:

            if score not in alignements:
                alignements[score] = set()
            alignements[score].add((qL, sL, qR, sR, hit[0]))

            print(f"Alignement trouvé dans la séquence {hit[0]} avec un score de {score}")
            print(f"Sequence 1: {sequence[qL:qR+1]} (positions {qL} à {qR})")
            print(f"Sequence 2: {database[hit[0]][sL:sR+1]} (positions {sL} à {sR})")
            print()


    alignements = dict(sorted(alignements.items(), reverse=True))
    print(alignements)
    
 


# Protéines
# Ungapped
def construire_mots_initiaux(sequence, w):
    mots_initiaux = {}
    nb_mots = len(sequence) - w + 1
    for i in range(nb_mots):
        mot = sequence[i:i+w]
        if mot not in mots_initiaux:
            mots_initiaux[mot] = []
        mots_initiaux[mot].append(i)
    return mots_initiaux

def generate_neighbors(mot, matrice_substitution, seuil_t):
    score = sum(matrice_substitution[lettre][lettre] for lettre in mot)
    voisins = []
    w = len(mot)
    if score >= seuil_t:
        for i in range(w):
            score_ligne = matrice_substitution[mot[i]]
            score_ligne = dict(sorted(score_ligne.items(), key=lambda item: item[1], reverse=True))
            for lettre, score_lettre in score_ligne.items():
                if lettre != mot[i]:
                    nouveau_mot = mot[:i] + lettre + mot[i+1:]
                    nouveau_score = score - matrice_substitution[mot[i]][mot[i]] + score_lettre
                    if nouveau_score >= seuil_t:
                        if nouveau_mot not in voisins:
                            voisins.append(nouveau_mot)
                    else:
                        break
    return voisins

def construire_liste_voisins(sequence, w, matrice_substitution, seuil_t):
    mots_initiaux = construire_mots_initiaux(sequence, w)
    voisins = {}
    for mot in mots_initiaux:
        pos = mots_initiaux[mot]
        liste_voisins = generate_neighbors(mot, matrice_substitution, seuil_t)
        for voisin in liste_voisins:
            if voisin not in voisins:
                voisins[voisin] = set()
            voisins[voisin].update(pos)

    return voisins
   
def run_proteines_ungapped(w, n, x_drop, seuil_t, seuil, taille_database):
    alphabete_proteines = 'ACDEFGHIKLMNPQRSTVWY'
    nb_mots = n - w + 1
    database = [''.join(random.choices(alphabete_proteines, k=random.randint(1000, 5000))) for _ in range(taille_database)]

    sequence = ''.join(random.choices(alphabete_proteines, k=n))

    # TODO : import matrices de substitutions
    PAM250 = {
        'A': {'A':  2, 'R': -2, 'N':  0, 'D':  0, 'C': -2, 'Q':  0, 'E':  0, 'G':  1, 'H': -1, 'I': -1,
            'L': -2, 'K': -1, 'M': -1, 'F': -3, 'P':  1, 'S':  1, 'T':  1, 'W': -6, 'Y': -3, 'V':  0},
        'R': {'A': -2, 'R':  6, 'N':  0, 'D': -1, 'C': -4, 'Q':  1, 'E': -1, 'G': -3, 'H':  2, 'I': -2,
            'L': -3, 'K':  3, 'M':  0, 'F': -4, 'P':  0, 'S':  0, 'T': -1, 'W':  2, 'Y': -4, 'V': -2},
        'N': {'A':  0, 'R':  0, 'N':  2, 'D':  2, 'C': -4, 'Q':  1, 'E':  1, 'G':  0, 'H':  2, 'I': -2,
            'L': -3, 'K':  1, 'M': -2, 'F': -3, 'P':  0, 'S':  1, 'T':  0, 'W': -4, 'Y': -2, 'V': -2},
        'D': {'A':  0, 'R': -1, 'N':  2, 'D':  4, 'C': -5, 'Q':  2, 'E':  3, 'G':  1, 'H':  1, 'I': -2,
            'L': -4, 'K':  0, 'M': -3, 'F': -6, 'P': -1, 'S':  0, 'T':  0, 'W': -7, 'Y': -4, 'V': -2},
        'C': {'A': -2, 'R': -4, 'N': -4, 'D': -5, 'C': 12, 'Q': -5, 'E': -5, 'G': -3, 'H': -3, 'I': -2,
            'L': -6, 'K': -5, 'M': -5, 'F': -4, 'P': -3, 'S':  0, 'T': -2, 'W': -8, 'Y':  0, 'V': -2},
        'Q': {'A':  0, 'R':  1, 'N':  1, 'D':  2, 'C': -5, 'Q':  4, 'E':  2, 'G': -1, 'H':  3, 'I': -2,
            'L': -2, 'K':  1, 'M': -1, 'F': -5, 'P':  0, 'S': -1, 'T': -1, 'W': -5, 'Y': -4, 'V': -2},
        'E': {'A':  0, 'R': -1, 'N':  1, 'D':  3, 'C': -5, 'Q':  2, 'E':  4, 'G':  0, 'H':  1, 'I': -2,
            'L': -3, 'K':  0, 'M': -2, 'F': -5, 'P': -1, 'S':  0, 'T':  0, 'W': -7, 'Y': -4, 'V': -2},
        'G': {'A':  1, 'R': -3, 'N':  0, 'D':  1, 'C': -3, 'Q': -1, 'E':  0, 'G':  5, 'H': -2, 'I': -3,
            'L': -4, 'K': -2, 'M': -3, 'F': -5, 'P':  0, 'S':  1, 'T':  0, 'W': -7, 'Y': -5, 'V': -1},
        'H': {'A': -1, 'R':  2, 'N':  2, 'D':  1, 'C': -3, 'Q':  3, 'E':  1, 'G': -2, 'H':  6, 'I': -2,
            'L': -2, 'K':  0, 'M': -2, 'F': -2, 'P':  0, 'S': -1, 'T': -1, 'W': -3, 'Y':  0, 'V': -2},
        'I': {'A': -1, 'R': -2, 'N': -2, 'D': -2, 'C': -2, 'Q': -2, 'E': -2, 'G': -3, 'H': -2, 'I':  5,
            'L':  2, 'K': -2, 'M':  2, 'F':  1, 'P': -2, 'S': -1, 'T':  0, 'W': -5, 'Y': -1, 'V':  4},
        'L': {'A': -2, 'R': -3, 'N': -3, 'D': -4, 'C': -6, 'Q': -2, 'E': -3, 'G': -4, 'H': -2, 'I':  2,
            'L':  6, 'K': -3, 'M':  4, 'F':  2, 'P': -3, 'S': -3, 'T': -2, 'W': -2, 'Y': -1, 'V':  2},
        'K': {'A': -1, 'R':  3, 'N':  1, 'D':  0, 'C': -5, 'Q':  1, 'E':  0, 'G': -2, 'H':  0, 'I': -2,
            'L': -3, 'K':  5, 'M':  0, 'F': -5, 'P': -1, 'S':  0, 'T':  0, 'W': -3, 'Y': -4, 'V': -2},
        'M': {'A': -1, 'R':  0, 'N': -2, 'D': -3, 'C': -5, 'Q': -1, 'E': -2, 'G': -3, 'H': -2, 'I':  2,
            'L':  4, 'K':  0, 'M':  6, 'F':  0, 'P': -2, 'S': -2, 'T': -1, 'W': -4, 'Y': -2, 'V':  2},
        'F': {'A': -3, 'R': -4, 'N': -3, 'D': -6, 'C': -4, 'Q': -5, 'E': -5, 'G': -5, 'H': -2, 'I':  1,
            'L':  2, 'K': -5, 'M':  0, 'F':  9, 'P': -5, 'S': -3, 'T': -3, 'W':  0, 'Y':  7, 'V': -1},
        'P': {'A':  1, 'R':  0, 'N':  0, 'D': -1, 'C': -3, 'Q':  0, 'E': -1, 'G':  0, 'H':  0, 'I': -2,
            'L': -3, 'K': -1, 'M': -2, 'F': -5, 'P':  6, 'S':  1, 'T':  0, 'W': -6, 'Y': -5, 'V': -1},
        'S': {'A':  1, 'R':  0, 'N':  1, 'D':  0, 'C':  0, 'Q': -1, 'E':  0, 'G':  1, 'H': -1, 'I': -1,
            'L': -3, 'K':  0, 'M': -2, 'F': -3, 'P':  1, 'S':  2, 'T':  1, 'W': -2, 'Y': -3, 'V': -1},
        'T': {'A':  1, 'R': -1, 'N':  0, 'D':  0, 'C': -2, 'Q': -1, 'E':  0, 'G':  0, 'H': -1, 'I':  0,
            'L': -2, 'K':  0, 'M': -1, 'F': -3, 'P':  0, 'S':  1, 'T':  3, 'W': -5, 'Y': -3, 'V':  0},
        'W': {'A': -6, 'R':  2, 'N': -4, 'D': -7, 'C': -8, 'Q': -5, 'E': -7, 'G': -7, 'H': -3, 'I': -5,
            'L': -2, 'K': -3, 'M': -4, 'F':  0, 'P': -6, 'S': -2, 'T': -5, 'W': 17, 'Y':  0, 'V': -6},
        'Y': {'A': -3, 'R': -4, 'N': -2, 'D': -4, 'C':  0, 'Q': -4, 'E': -4, 'G': -5, 'H':  0, 'I': -1,
            'L': -1, 'K': -4, 'M': -2, 'F':  7, 'P': -5, 'S': -3, 'T': -3, 'W':  0, 'Y': 10, 'V': -2},
        'V': {'A':  0, 'R': -2, 'N': -2, 'D': -2, 'C': -2, 'Q': -2, 'E': -2, 'G': -1, 'H': -2, 'I':  4,
            'L':  2, 'K': -2, 'M':  2, 'F': -1, 'P': -1, 'S': -1, 'T':  0, 'W': -6, 'Y': -2, 'V':  4}
    }
    
    liste = construire_liste_voisins(sequence, w, PAM250, seuil_t)
    print(f"Nombre de voisins générés : {len(liste)}")
    
    hits = find_hits(liste, database, w)
    print(f"Nombre de hits trouvés : {len(hits)}")
    
    alignements = {}
        
    for hit in hits:
        alignment = extend_hit(hit, sequence, database, w, PAM250, x_drop)
        qL, sL, qR, sR, score = alignment
        if score >= seuil:

            if score not in alignements:
                alignements[score] = set()
            alignements[score].add((qL, sL, qR, sR, hit[0]))

            print(f"Alignement trouvé dans la séquence {hit[0]} avec un score de {score}")
            print(f"Sequence 1: {sequence[qL:qR+1]} (positions {qL} à {qR})")
            print(f"Sequence 2: {database[hit[0]][sL:sR+1]} (positions {sL} à {sR})")
            print()


    alignements = dict(sorted(alignements.items(), reverse=True))
    print(len(alignements))
    for score in alignements:
        print(f"Score: {score}, Nombre d'alignements: {len(alignements[score])}")
        
# Gapped
def find_double_hits(dico_mots, database, w, A):
    doubles_hits = set()
    
    for seq_id, seq in enumerate(database):
        last_hit = {}
        for pos_in_seq in range(len(seq) - w + 1):
            mot = seq[pos_in_seq:pos_in_seq + w]
            if mot in dico_mots:
                for pos_in_querry in dico_mots[mot]:
                    diagonal = pos_in_querry - pos_in_seq
                    if diagonal in last_hit:
                        if w <= pos_in_querry - last_hit[diagonal] <= A:
                            doubles_hits.add((seq_id, pos_in_querry, pos_in_seq)) # On ajoute que le second hits, que l'on va étendre
                    last_hit[diagonal] = pos_in_querry
    return doubles_hits

def find_hits_need_extension(list_hit, sequence, database, w, M, Sg, x_drop, seed_window_length = 11):
    hits_nead_extension = set()
    for hit in list_hit:
        seq_id, _, _ = hit
        best_qL, best_sL, best_qR, best_sR, total_best = extend_hit(hit, sequence, database, w, M, x_drop) 
        if total_best >= Sg:
            best_seed_q = best_qL
            best_seed_s = best_sL
            qL = best_qL
            sL = best_sL
            if best_qR - best_qL + 1 < seed_window_length:
                seed_q = (best_qL + best_qR) // 2
                seed_s = (best_sL + best_sR) // 2
                hits_nead_extension.add((seq_id, seed_q, seed_s))
            score = sum(M[sequence[best_qL + i]][database[seq_id][best_sL + i]] for i in range(seed_window_length))
            best_score = score
            for i in range(best_qR - best_qL - seed_window_length + 1):
                qL = qL + 1
                sL = sL + 1
                score = score - M[sequence[qL - 1]][database[seq_id][sL - 1]] + M[sequence[qL + seed_window_length - 1]][database[seq_id][sL + seed_window_length - 1]]
                if score > best_score:
                    best_score = score
                    best_seed_q = qL
                    best_seed_s = sL
            seed_s = best_seed_s + seed_window_length // 2
            seed_q = best_seed_q + seed_window_length // 2
            if best_score >= Sg:
                hits_nead_extension.add((seq_id, seed_q, seed_s))
    return hits_nead_extension


def run_proteines_gapped(w, n, taille_database, seuil_t, A, x_drop, Sg):
    alphabete_proteines = 'ACDEFGHIKLMNPQRSTVWY'
    database = [''.join(random.choices(alphabete_proteines, k=random.randint(1000, 5000))) for _ in range(taille_database)]

    sequence = ''.join(random.choices(alphabete_proteines, k=n))

    # TODO : import matrices de substitutions
    PAM250 = {
        'A': {'A':  2, 'R': -2, 'N':  0, 'D':  0, 'C': -2, 'Q':  0, 'E':  0, 'G':  1, 'H': -1, 'I': -1,
            'L': -2, 'K': -1, 'M': -1, 'F': -3, 'P':  1, 'S':  1, 'T':  1, 'W': -6, 'Y': -3, 'V':  0},
        'R': {'A': -2, 'R':  6, 'N':  0, 'D': -1, 'C': -4, 'Q':  1, 'E': -1, 'G': -3, 'H':  2, 'I': -2,
            'L': -3, 'K':  3, 'M':  0, 'F': -4, 'P':  0, 'S':  0, 'T': -1, 'W':  2, 'Y': -4, 'V': -2},
        'N': {'A':  0, 'R':  0, 'N':  2, 'D':  2, 'C': -4, 'Q':  1, 'E':  1, 'G':  0, 'H':  2, 'I': -2,
            'L': -3, 'K':  1, 'M': -2, 'F': -3, 'P':  0, 'S':  1, 'T':  0, 'W': -4, 'Y': -2, 'V': -2},
        'D': {'A':  0, 'R': -1, 'N':  2, 'D':  4, 'C': -5, 'Q':  2, 'E':  3, 'G':  1, 'H':  1, 'I': -2,
            'L': -4, 'K':  0, 'M': -3, 'F': -6, 'P': -1, 'S':  0, 'T':  0, 'W': -7, 'Y': -4, 'V': -2},
        'C': {'A': -2, 'R': -4, 'N': -4, 'D': -5, 'C': 12, 'Q': -5, 'E': -5, 'G': -3, 'H': -3, 'I': -2,
            'L': -6, 'K': -5, 'M': -5, 'F': -4, 'P': -3, 'S':  0, 'T': -2, 'W': -8, 'Y':  0, 'V': -2},
        'Q': {'A':  0, 'R':  1, 'N':  1, 'D':  2, 'C': -5, 'Q':  4, 'E':  2, 'G': -1, 'H':  3, 'I': -2,
            'L': -2, 'K':  1, 'M': -1, 'F': -5, 'P':  0, 'S': -1, 'T': -1, 'W': -5, 'Y': -4, 'V': -2},
        'E': {'A':  0, 'R': -1, 'N':  1, 'D':  3, 'C': -5, 'Q':  2, 'E':  4, 'G':  0, 'H':  1, 'I': -2,
            'L': -3, 'K':  0, 'M': -2, 'F': -5, 'P': -1, 'S':  0, 'T':  0, 'W': -7, 'Y': -4, 'V': -2},
        'G': {'A':  1, 'R': -3, 'N':  0, 'D':  1, 'C': -3, 'Q': -1, 'E':  0, 'G':  5, 'H': -2, 'I': -3,
            'L': -4, 'K': -2, 'M': -3, 'F': -5, 'P':  0, 'S':  1, 'T':  0, 'W': -7, 'Y': -5, 'V': -1},
        'H': {'A': -1, 'R':  2, 'N':  2, 'D':  1, 'C': -3, 'Q':  3, 'E':  1, 'G': -2, 'H':  6, 'I': -2,
            'L': -2, 'K':  0, 'M': -2, 'F': -2, 'P':  0, 'S': -1, 'T': -1, 'W': -3, 'Y':  0, 'V': -2},
        'I': {'A': -1, 'R': -2, 'N': -2, 'D': -2, 'C': -2, 'Q': -2, 'E': -2, 'G': -3, 'H': -2, 'I':  5,
            'L':  2, 'K': -2, 'M':  2, 'F':  1, 'P': -2, 'S': -1, 'T':  0, 'W': -5, 'Y': -1, 'V':  4},
        'L': {'A': -2, 'R': -3, 'N': -3, 'D': -4, 'C': -6, 'Q': -2, 'E': -3, 'G': -4, 'H': -2, 'I':  2,
            'L':  6, 'K': -3, 'M':  4, 'F':  2, 'P': -3, 'S': -3, 'T': -2, 'W': -2, 'Y': -1, 'V':  2},
        'K': {'A': -1, 'R':  3, 'N':  1, 'D':  0, 'C': -5, 'Q':  1, 'E':  0, 'G': -2, 'H':  0, 'I': -2,
            'L': -3, 'K':  5, 'M':  0, 'F': -5, 'P': -1, 'S':  0, 'T':  0, 'W': -3, 'Y': -4, 'V': -2},
        'M': {'A': -1, 'R':  0, 'N': -2, 'D': -3, 'C': -5, 'Q': -1, 'E': -2, 'G': -3, 'H': -2, 'I':  2,
            'L':  4, 'K':  0, 'M':  6, 'F':  0, 'P': -2, 'S': -2, 'T': -1, 'W': -4, 'Y': -2, 'V':  2},
        'F': {'A': -3, 'R': -4, 'N': -3, 'D': -6, 'C': -4, 'Q': -5, 'E': -5, 'G': -5, 'H': -2, 'I':  1,
            'L':  2, 'K': -5, 'M':  0, 'F':  9, 'P': -5, 'S': -3, 'T': -3, 'W':  0, 'Y':  7, 'V': -1},
        'P': {'A':  1, 'R':  0, 'N':  0, 'D': -1, 'C': -3, 'Q':  0, 'E': -1, 'G':  0, 'H':  0, 'I': -2,
            'L': -3, 'K': -1, 'M': -2, 'F': -5, 'P':  6, 'S':  1, 'T':  0, 'W': -6, 'Y': -5, 'V': -1},
        'S': {'A':  1, 'R':  0, 'N':  1, 'D':  0, 'C':  0, 'Q': -1, 'E':  0, 'G':  1, 'H': -1, 'I': -1,
            'L': -3, 'K':  0, 'M': -2, 'F': -3, 'P':  1, 'S':  2, 'T':  1, 'W': -2, 'Y': -3, 'V': -1},
        'T': {'A':  1, 'R': -1, 'N':  0, 'D':  0, 'C': -2, 'Q': -1, 'E':  0, 'G':  0, 'H': -1, 'I':  0,
            'L': -2, 'K':  0, 'M': -1, 'F': -3, 'P':  0, 'S':  1, 'T':  3, 'W': -5, 'Y': -3, 'V':  0},
        'W': {'A': -6, 'R':  2, 'N': -4, 'D': -7, 'C': -8, 'Q': -5, 'E': -7, 'G': -7, 'H': -3, 'I': -5,
            'L': -2, 'K': -3, 'M': -4, 'F':  0, 'P': -6, 'S': -2, 'T': -5, 'W': 17, 'Y':  0, 'V': -6},
        'Y': {'A': -3, 'R': -4, 'N': -2, 'D': -4, 'C':  0, 'Q': -4, 'E': -4, 'G': -5, 'H':  0, 'I': -1,
            'L': -1, 'K': -4, 'M': -2, 'F':  7, 'P': -5, 'S': -3, 'T': -3, 'W':  0, 'Y': 10, 'V': -2},
        'V': {'A':  0, 'R': -2, 'N': -2, 'D': -2, 'C': -2, 'Q': -2, 'E': -2, 'G': -1, 'H': -2, 'I':  4,
            'L':  2, 'K': -2, 'M':  2, 'F': -1, 'P': -1, 'S': -1, 'T':  0, 'W': -6, 'Y': -2, 'V':  4}
    }
    
    liste = construire_liste_voisins(sequence, w, PAM250, seuil_t)
    print(f"Nombre de voisins générés : {len(liste)}")
    
    hits_simple = find_hits(liste, database, w)
    print(f"Nombre de hits simples trouvés : {len(hits_simple)}")
    
    hits_double = find_double_hits(liste, database, w, A)
    print(f"Nombre de doubles hits trouvés : {len(hits_double)}")

    hits_need_extension = find_hits_need_extension(hits_double, sequence, database, w, PAM250, Sg, x_drop)
    print(f"Nombre de hits à étendre trouvés : {len(hits_need_extension)}")
if __name__ == "__main__":
    # ADN aléatoire
    # Paramètres

    w_adn_ungapped = 13  # taille des mots
    n_adn_ungapped = 250 # taille de la séquence
    nb_mots = n_adn_ungapped - w_adn_ungapped + 1
    x_drop_adn_ungapped = 20
    seuil_adn_ungapped = 50 # seuil de score pour garder un alignement
    taille_database_adn_ungapped = 1000 # taille de la base de données


    # Protéines aléatoires
    # Paramètres
    w_proteines_ungapped = 3     # taille des mots
    n_proteines_ungapped = 250 # taille de la séquence
    seuil_t_proteines_ungapped = 13 # seuil de score pour garder un alignement
    seuil_proteines_ungapped = 80 # seuil de score pour garder un alignement
    x_drop_proteines_ungapped = 20
    taille_database_proteines_ungapped = 1000 # taille de la base de données


    # Gapped
    w_proteines_gapped = 3    # taille des mots
    n_proteines_gapped = 250 # taille de la séquence
    seuil_t_proteines_gapped = 11 # seuil de score pour garder un alignement
    seuil_proteines_gapped = 80 # seuil de score pour garder un alignement
    A_gapped = 40
    x_drop_proteines_gapped = 16
    taille_database_proteines_gapped = 1000 # taille de la base de données
    Sg = 42 # seuil de score pour trigger une gapped extension


    run_proteines_gapped(w_proteines_gapped, n_proteines_gapped, taille_database_proteines_gapped, seuil_t_proteines_gapped, A_gapped, x_drop_proteines_gapped, Sg)