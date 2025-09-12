""" Exemple basique, juste sur l'ADN aléatoire """
""" ------------------------------------------ """


import random
import time
import matplotlib.pyplot as plt


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
    
    
    
# Fréquences de Robinson & Robinson (1991)
aa_freqs = {
    "A": 0.07805,  # Alanine
    "C": 0.01925,  # Cystéine
    "D": 0.05364,  # Aspartate
    "E": 0.06295,  # Glutamate
    "F": 0.03856,  # Phénylalanine
    "G": 0.07377,  # Glycine
    "H": 0.02199,  # Histidine
    "I": 0.05142,  # Isoleucine
    "K": 0.05744,  # Lysine
    "L": 0.09019,  # Leucine
    "M": 0.02243,  # Méthionine
    "N": 0.04487,  # Asparagine
    "P": 0.05203,  # Proline
    "Q": 0.04264,  # Glutamine
    "R": 0.05130,  # Arginine
    "S": 0.07120,  # Sérine
    "T": 0.05841,  # Thréonine
    "V": 0.06441,  # Valine
    "W": 0.01330,  # Tryptophane
    "Y": 0.03216   # Tyrosine
}

alphabet_proteines = ''.join(aa_freqs.keys())
weights = list(aa_freqs.values())

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

    database = [
        ''.join(random.choices(alphabet_proteines, weights=weights, k=random.randint(1000, 5000)))
        for _ in range(taille_database)
    ]
    
    sequence = ''.join(random.choices(alphabet_proteines, weights=weights, k=n))

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

def find_hits_need_extension(list_hit, query, database, w, M, Sg, x_drop, seed_window_length = 11):
    hits_nead_extension = set()
    for hit in list_hit:
        seq_id, _, _ = hit
        best_qL, best_sL, best_qR, best_sR, total_best = extend_hit(hit, query, database, w, M, x_drop) 
        if total_best >= Sg:
            best_seed_q = best_qL
            best_seed_s = best_sL
            qL = best_qL
            sL = best_sL
            if best_qR - best_qL + 1 < seed_window_length:
                seed_q = (best_qL + best_qR) // 2
                seed_s = (best_sL + best_sR) // 2
                hits_nead_extension.add((seq_id, seed_q, seed_s))
                continue
            score = sum(M[query[best_qL + i]][database[seq_id][best_sL + i]] for i in range(seed_window_length))
            best_score = score
            for i in range(best_qR - best_qL - seed_window_length + 1):
                qL = qL + 1
                sL = sL + 1
                score = score - M[query[qL - 1]][database[seq_id][sL - 1]] + M[query[qL + seed_window_length - 1]][database[seq_id][sL + seed_window_length - 1]]
                if score > best_score:
                    best_score = score
                    best_seed_q = qL
                    best_seed_s = sL
            seed_s = best_seed_s + seed_window_length // 2
            seed_q = best_seed_q + seed_window_length // 2
            if best_score >= Sg:
                hits_nead_extension.add((seq_id, seed_q, seed_s))
    return hits_nead_extension

def extend_gapped(hit, query, database, M, x_drop, ouverture=10, extension=1, direction="right"):
    """
    Extension gappée type BLAST II avec X-drop.
    
    Args:
        hit : (seq_id, qi, si) positions du seed
        query, database : séquences
        M : matrice de substitution, indexée par lettres
        x_drop : seuil d'élagage
        ouverture, extension : coûts des gaps
        direction : "right" ou "left"
    
    Returns:
        (score_max, q_pos, s_pos)
    """
    seq_id, qi, si = hit
    seq_db = database[seq_id]

    # longueur dispo selon la direction
    if direction == "right":
        n, m = len(query) - qi, len(seq_db) - si
        q_start, s_start = qi, si

    elif direction == "left":
        # inverser les séquences
        seq_db = seq_db[::-1]
        query = query[::-1]

        # positions du seed dans les séquences inversées
        qi = len(query) - 1 - qi
        si = len(seq_db) - 1 - si

        # longueurs dispo
        n, m = len(query) - qi, len(seq_db) - si

        q_start, s_start = qi, si
    else:
        raise ValueError("direction must be 'right' or 'left'")

    # états : dictionnaires "colonne -> valeur"
    H_prev, E_prev, F_prev = {}, {}, {}
    H_curr, E_curr, F_curr = {}, {}, {}

    # init au seed (0,0 dans repère local)
    H_prev[0] = M[query[q_start]][seq_db[s_start]]
    E_prev[0] = float("-inf")
    F_prev[0] = float("-inf")

    max_score = H_prev[0]
    best_q, best_s = qi, si

    active_cols = set({0, 1})
    

    for i in range(1, n):
        H_curr.clear()
        E_curr.clear()
        F_curr.clear()
        new_active_cols = set()

        for j in active_cols:
            if j+1 > m:
                continue

            # gap vertical (insertion dans DB → avance query)
            e_val = max(H_prev.get(j, float("-inf")) - ouverture,
                        E_prev.get(j, float("-inf")) - extension)
            E_curr[j] = e_val

            # gap horizontal (insertion dans query → avance DB)
            f_val = max(
                H_curr.get(j-1, float("-inf")) - ouverture,
                F_curr.get(j-1, float("-inf")) - extension
            )
            F_curr[j] = f_val

            # match / substitution
            if j-1 in H_prev:
                q_idx = q_start + i
                s_idx = s_start + j

                sub = H_prev[j-1] + M[query[q_idx]][seq_db[s_idx]]
            else:
                sub = float("-inf")

            h_val = max(sub, e_val, f_val)
            H_curr[j] = h_val

            if h_val > max_score:
                max_score = h_val
                best_q = q_start + i
                best_s = s_start + j

            # X-drop pruning
            if h_val >= max_score - x_drop:
                new_active_cols.add(j)
                if j+1 <= m:
                    new_active_cols.add(j+1)

        H_prev, E_prev, F_prev = H_curr.copy(), E_curr.copy(), F_curr.copy()
        active_cols = new_active_cols
    
        if not active_cols:
            break
        
    if direction == "left":
        best_q = len(query) - 1 - best_q
        best_s = len(seq_db) - 1 - best_s

    return max_score, best_q, best_s

def find_alignments_gapped(hits, query, database, M, seuil_accept, x_drop):
    results = set()
    for hit in hits:
        score_right, best_q_right, best_s_right = extend_gapped(hit, query, database, M, x_drop, 10, 1, "right")
        score_left, best_q_left, best_s_left = extend_gapped(hit, query, database, M, x_drop, 10, 1, "left")
        total_score = score_left + score_right - M[query[hit[1]]][database[hit[0]][hit[2]]]
        if total_score >= seuil_accept:
            results.add((total_score, hit[0], best_q_left, best_q_right, best_s_left, best_s_right))
    return results

def run_proteines_gapped(w, n, taille_database, seuil_t, seuil_proteines_gapped, A, x_drop, x_drop_gapped, Sg):
    
    database = [
        ''.join(random.choices(alphabet_proteines, weights=weights, k=random.randint(1000, 5000)))
        for _ in range(taille_database)
    ]
    
    sequence = ''.join(random.choices(alphabet_proteines, weights=weights, k=n))
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
    
    alignements = find_alignments_gapped(hits_need_extension, sequence, database, PAM250, seuil_proteines_gapped, x_drop_gapped)
    print(f"Nombre d'alignements gappés trouvés : {len(alignements)}")
    alignements = sorted(alignements, reverse=True)
    sc = {}
    for score, seq_id, qL, qR, sL, sR in alignements:
        if score not in sc:
            sc[score] = set()
        sc[score].add((seq_id))
    sc = dict(sorted(sc.items(), reverse=True))
    for s in sc:
        print(f"Score: {s}, Nombre d'alignements: {len(sc[s])}")
        



# -------- TEST 1 : Hits simples / doubles selon T --------
def test_hits_vs_T(w=3, n=250, taille_db=200,
                   seuils_T=[9, 10, 11, 12, 13, 14, 15], A=40, delta=3):
    database = [
        ''.join(random.choices(alphabet_proteines, weights=weights, k=random.randint(1000, 3000)))
        for _ in range(taille_db)
    ]
    query = ''.join(random.choices(alphabet_proteines, weights=weights, k=n))

    results = []
    for T in seuils_T:
        voisins = construire_liste_voisins(query, w, PAM250, T)
        hits_simple = find_hits(voisins, database, w)
        hits_double = find_double_hits(voisins, database, w, A)

        # Hits simples avec T décalé (pour comparaison BLAST I vs II)
        voisins_shift = construire_liste_voisins(query, w, PAM250, T + delta)
        hits_simple_shift = find_hits(voisins_shift, database, w)

        results.append((T, len(hits_simple), len(hits_double), len(hits_simple_shift)))

    # Affichage des ratios
    print("Facteurs de variation (par rapport au seuil précédent) :")
    for i, (T, h, d, h_shift) in enumerate(results):
        ratio_double_simple = d / max(h, 1)
        ratio_double_shift = d / max(h_shift, 1)
        ratio_simple_shift = h / max(h_shift, 1)

        if i == 0:
            print(f"T={T} : Hits simples={h}, Doubles hits={d}, "
                  f"Hits simples (T+{delta})={h_shift}, "
                  f"Ratio double/simple={ratio_double_simple:.3f}, "
                  f"Ratio double/(simple T+{delta})={ratio_double_shift:.3f} , "
                  f"Ratio simple/(simple T+{delta})={ratio_simple_shift:.3f}")
        else:
            T_prev, h_prev, d_prev, h_shift_prev = results[i-1]
            red_hits = h / max(h_prev, 1)
            red_doubles = d / max(d_prev, 1)
            red_hits_shift = h_shift / max(h_shift_prev, 1)

            print(f"T={T} : "
                  f"Hits simples ×{red_hits:.2f}, "
                  f"Doubles hits ×{red_doubles:.2f}, "
                  f"Hits simples (T+{delta}) ×{red_hits_shift:.2f}, "
                  f"Ratio double/simple={ratio_double_simple:.3f}, "
                  f"Ratio double/(simple T+{delta})={ratio_double_shift:.3f}, "
                  f"Ratio simple/(simple T+{delta})={ratio_simple_shift:.3f}")

    # Plot
    T_vals, hits_vals, doubles_vals, hits_shift_vals = zip(*results)
    plt.figure(figsize=(8, 5))
    plt.plot(T_vals, hits_vals, marker="s", label="Hits simples (BLAST I et II)")
    plt.plot(T_vals, hits_shift_vals, marker="o", linestyle="--",
             label=f"Hits simples (T+{delta}, BLAST I)")
    plt.plot(T_vals, doubles_vals, marker="^", label="Doubles hits (BLAST II)")
    plt.xlabel("Seuil T")
    plt.ylabel("Nombre de hits (log)")
    plt.yscale("log")
    plt.title("Impact du seuil T sur hits simples/doubles (BLAST I vs BLAST II)")
    plt.legend()
    plt.grid(True, which="both", ls="--", alpha=0.7)
    plt.show()


def compare_acceptance(w=3, n=250, taille_db=3000,
                       T_blast1=13, T_blast2=11,
                       A=30, x_drop=16, x_drop_gapped=40, Sg=42,
                       S_accept_range=range(60, 81, 2)):
    results_b1, results_b2 = [], []

    # Génération DB + query
    database = [
        ''.join(random.choices(alphabet_proteines, weights=weights, k=random.randint(1000, 3000)))
        for _ in range(taille_db)
    ]
    query = ''.join(random.choices(alphabet_proteines, weights=weights, k=n))

    print(f"\n=== Test comparaison BLAST I vs BLAST II ===")
    print(f"Taille DB = {taille_db} séquences, longueur query = {n}")
    print(f"S_accept testé de {min(S_accept_range)} à {max(S_accept_range)} (Sg={Sg})\n")

    # ---- BLAST I ----
    voisins1 = construire_liste_voisins(query, w, PAM250, T_blast1)
    hits_simple = find_hits(voisins1, database, w)
    align_blast1_all = [(hit[0], extend_hit(hit, query, database, w, PAM250, x_drop)[-1])
                        for hit in hits_simple]

    # ---- BLAST II ----
    voisins2 = construire_liste_voisins(query, w, PAM250, T_blast2)
    hits_double = find_double_hits(voisins2, database, w, A)
    hits_ext = find_hits_need_extension(hits_double, query, database, w, PAM250, Sg, x_drop)
    align_blast2_all = find_alignments_gapped(hits_ext, query, database, PAM250, Sg, x_drop_gapped)
    # -> on récupère tous les alignements gappés déclenchés avec Sg

    # ---- Boucle sur les seuils ----
    for S_accept in S_accept_range:
        # BLAST I filtré
        retained_b1 = [(sid, sc) for sid, sc in align_blast1_all if sc >= S_accept]
        results_b1.append(len(retained_b1))

        # BLAST II filtré
        retained_b2 = [(aln[1], aln[0]) for aln in align_blast2_all if aln[0] >= S_accept]
        results_b2.append(len(retained_b2))

        # Impression terminal
        print(f"S_accept = {S_accept}")
        print(f"  BLAST I : {len(retained_b1)} alignements retenus")
        if retained_b1:
            print("    Top 5 scores :", sorted(retained_b1, key=lambda x: -x[1])[:5])
        print(f"  BLAST II : {len(retained_b2)} alignements retenus")
        if retained_b2:
            print("    Top 5 scores :", sorted(retained_b2, key=lambda x: -x[1])[:5])
        print()

    # ---- Plot ----
    plt.figure(figsize=(8,5))
    plt.plot(S_accept_range, results_b1, marker="o", label=f"BLAST I (T={T_blast1})")
    plt.plot(S_accept_range, results_b2, marker="s", label=f"BLAST II (T={T_blast2}, double-hit, Sg={Sg})")
    plt.xlabel("Seuil d'acceptance $S_{accept}$")
    plt.ylabel("Nombre d'alignements avec score ≥ $S_{accept}$")
    plt.title("Comparaison BLAST I vs BLAST II en fonction de $S_{accept}$")
    plt.yscale("log")
    plt.grid(True, which="both", ls="--", alpha=0.7)
    plt.legend()
    plt.show()

    return results_b1, results_b2


import time
import matplotlib.pyplot as plt

# -------- TEST 3 : Temps d’exécution vs taille base --------
def test_time_vs_dbsize(
    w=3, n=250,
    T1=13, T2=11,    # T différents pour BLAST I et II (sensibilité équivalente)
    A=40,
    x_drop=20, x_drop_gapped=40,
    Sg=42, S_accept=60,
    db_sizes=[500, 1000, 1500, 2000, 2500, 3000]
):
    query = ''.join(random.choices(alphabet_proteines, weights=weights, k=n))

    times_blast1, times_blast2 = [], []

    for size in db_sizes:
        database = [
            ''.join(random.choices(alphabet_proteines, weights=weights, k=random.randint(1000, 3000)))
            for _ in range(size)
        ]

        # ---------- BLAST I ----------
        t0 = time.time()
        voisins1 = construire_liste_voisins(query, w, PAM250, T1)
        t1 = time.time()
        hits_simple = find_hits(voisins1, database, w)
        t2 = time.time()
        _ = [extend_hit(hit, query, database, w, PAM250, x_drop) for hit in hits_simple]
        t3 = time.time()

        total1 = t3 - t0
        step1 = t1 - t0
        step2 = t2 - t1
        step3 = t3 - t2
        times_blast1.append(total1)

        print(f"\n[BLAST I] DB size={size}, total={total1:.3f}s")
        print(f"  Step1 (voisins) : {step1:.3f}s ({100*step1/total1:.1f}%)")
        print(f"  Step2 (hits)    : {step2:.3f}s ({100*step2/total1:.1f}%)")
        print(f"  Step3 (extend)  : {step3:.3f}s ({100*step3/total1:.1f}%)")

        # ---------- BLAST II ----------
        t0 = time.time()
        voisins2 = construire_liste_voisins(query, w, PAM250, T2)
        t1 = time.time()
        hits_double = find_double_hits(voisins2, database, w, A)
        t2 = time.time()
        hits_ext = find_hits_need_extension(hits_double, query, database, w, PAM250, Sg, x_drop)
        align_blast2 = find_alignments_gapped(hits_ext, query, database, PAM250, S_accept, x_drop_gapped)
        t3 = time.time()

        total2 = t3 - t0
        step1 = t1 - t0
        step2 = t2 - t1
        step3 = t3 - t2
        times_blast2.append(total2)

        print(f"\n[BLAST II] DB size={size}, total={total2:.3f}s")
        print(f"  Step1 (voisins) : {step1:.3f}s ({100*step1/total2:.1f}%)")
        print(f"  Step2 (hits)    : {step2:.3f}s ({100*step2/total2:.1f}%)")
        print(f"  Step3 (extend)  : {step3:.3f}s ({100*step3/total2:.1f}%)")

    # Plot comparatif global
    plt.figure(figsize=(8, 5))
    plt.plot(db_sizes, times_blast1, marker="o", label=f"BLAST I (T={T1})")
    plt.plot(db_sizes, times_blast2, marker="s", label=f"BLAST II (T={T2}, A={A}, Sg={Sg})")
    plt.xlabel("Taille de la base (nombre de séquences)")
    plt.ylabel("Temps total (s)")
    plt.title("Temps d’exécution BLAST I vs BLAST II")
    plt.legend()
    plt.grid(True)
    plt.show()

    return times_blast1, times_blast2

if __name__ == "__main__":
    import argparse
    import sys
    import random

    # --- Profils par défaut (reprennent tes valeurs) ---
    defaults = {
        "adn": dict(
            w=13, n=250, x_drop=20, seuil=50, taille_database=1000
        ),
        "prot_ungapped": dict(
            w=3, n=250, x_drop=20, seuil_t=13, seuil=80, taille_database=1000
        ),
        "prot_gapped": dict(
            w=3, n=250, seuil_t=11, seuil=80, A=40,
            x_drop=20, x_drop_gapped=40, taille_database=1000, Sg=42
        ),
        "test1": dict(
            w=3, n=250, taille_db=200, seuils_T=[9,10,11,12,13,14,15], A=40, delta=3
        ),
        "test2": dict(
            w=3, n=250, taille_db=3000, T_blast1=13, T_blast2=11, A=30,
            x_drop=16, x_drop_gapped=40, Sg=42, S_accept_range=range(60, 81, 2)
        ),
        "test3": dict(
            w=3, n=250, T1=13, T2=11, A=40, x_drop=20, x_drop_gapped=40,
            Sg=42, S_accept=60, db_sizes=[500,1000,1500,2000,2500,3000]
        ),
    }

    parser = argparse.ArgumentParser(
        description="Runner BLAST I / BLAST II (tests & modes)."
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    # --- Sous-commandes : modes simples ---
    p_adn = sub.add_parser("adn", help="ADN ungapped")
    p_adn.add_argument("--w", type=int, default=defaults["adn"]["w"])
    p_adn.add_argument("--n", type=int, default=defaults["adn"]["n"])
    p_adn.add_argument("--x-drop", type=int, default=defaults["adn"]["x_drop"])
    p_adn.add_argument("--seuil", type=int, default=defaults["adn"]["seuil"])
    p_adn.add_argument("--db", type=int, default=defaults["adn"]["taille_database"])

    p_pu = sub.add_parser("prot-ungapped", help="Protéines ungapped")
    p_pu.add_argument("--w", type=int, default=defaults["prot_ungapped"]["w"])
    p_pu.add_argument("--n", type=int, default=defaults["prot_ungapped"]["n"])
    p_pu.add_argument("--x-drop", type=int, default=defaults["prot_ungapped"]["x_drop"])
    p_pu.add_argument("--T", type=int, default=defaults["prot_ungapped"]["seuil_t"])
    p_pu.add_argument("--seuil", type=int, default=defaults["prot_ungapped"]["seuil"])
    p_pu.add_argument("--db", type=int, default=defaults["prot_ungapped"]["taille_database"])

    p_pg = sub.add_parser("prot-gapped", help="Protéines gapped (double-hit + X-drop gappé)")
    p_pg.add_argument("--w", type=int, default=defaults["prot_gapped"]["w"])
    p_pg.add_argument("--n", type=int, default=defaults["prot_gapped"]["n"])
    p_pg.add_argument("--T", type=int, default=defaults["prot_gapped"]["seuil_t"])
    p_pg.add_argument("--seuil", type=int, default=defaults["prot_gapped"]["seuil"])
    p_pg.add_argument("--A", type=int, default=defaults["prot_gapped"]["A"])
    p_pg.add_argument("--x-drop", type=int, default=defaults["prot_gapped"]["x_drop"])
    p_pg.add_argument("--x-drop-gapped", type=int, default=defaults["prot_gapped"]["x_drop_gapped"])
    p_pg.add_argument("--db", type=int, default=defaults["prot_gapped"]["taille_database"])
    p_pg.add_argument("--Sg", type=int, default=defaults["prot_gapped"]["Sg"])

    # --- Sous-commandes : tests 1/2/3 ---
    p_t1 = sub.add_parser("test1", help="Hits simples/doubles vs T (BLAST I vs II)")
    p_t1.add_argument("--w", type=int, default=defaults["test1"]["w"])
    p_t1.add_argument("--n", type=int, default=defaults["test1"]["n"])
    p_t1.add_argument("--db", type=int, default=defaults["test1"]["taille_db"])
    p_t1.add_argument("--A", type=int, default=defaults["test1"]["A"])
    p_t1.add_argument("--delta", type=int, default=defaults["test1"]["delta"])
    p_t1.add_argument("--seuils-T", type=int, nargs="+", default=defaults["test1"]["seuils_T"])

    p_t2 = sub.add_parser("test2", help="Nb d'alignements vs S_accept (BLAST I vs II)")
    p_t2.add_argument("--w", type=int, default=defaults["test2"]["w"])
    p_t2.add_argument("--n", type=int, default=defaults["test2"]["n"])
    p_t2.add_argument("--db", type=int, default=defaults["test2"]["taille_db"])
    p_t2.add_argument("--T1", type=int, default=defaults["test2"]["T_blast1"])
    p_t2.add_argument("--T2", type=int, default=defaults["test2"]["T_blast2"])
    p_t2.add_argument("--A", type=int, default=defaults["test2"]["A"])
    p_t2.add_argument("--x-drop", type=int, default=defaults["test2"]["x_drop"])
    p_t2.add_argument("--x-drop-gapped", type=int, default=defaults["test2"]["x_drop_gapped"])
    p_t2.add_argument("--Sg", type=int, default=defaults["test2"]["Sg"])
    p_t2.add_argument("--S-accept", type=int, nargs="+", default=list(defaults["test2"]["S_accept_range"]))

    p_t3 = sub.add_parser("test3", help="Temps total vs taille base (BLAST I vs II)")
    p_t3.add_argument("--w", type=int, default=defaults["test3"]["w"])
    p_t3.add_argument("--n", type=int, default=defaults["test3"]["n"])
    p_t3.add_argument("--T1", type=int, default=defaults["test3"]["T1"])
    p_t3.add_argument("--T2", type=int, default=defaults["test3"]["T2"])
    p_t3.add_argument("--A", type=int, default=defaults["test3"]["A"])
    p_t3.add_argument("--x-drop", type=int, default=defaults["test3"]["x_drop"])
    p_t3.add_argument("--x-drop-gapped", type=int, default=defaults["test3"]["x_drop_gapped"])
    p_t3.add_argument("--Sg", type=int, default=defaults["test3"]["Sg"])
    p_t3.add_argument("--S-accept", type=int, default=defaults["test3"]["S_accept"])
    p_t3.add_argument("--db-sizes", type=int, nargs="+", default=defaults["test3"]["db_sizes"])

    # Options communes
    parser.add_argument("--seed", type=int, default=None, help="Graine RNG pour reproductibilité")

    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    # --- Dispatch ---
    cmd = args.cmd
    if cmd == "adn":
        run_adn_ungapped(
            w=args.w, n=args.n, x_drop=args.__dict__["x_drop"],
            seuil=args.seuil, taille_database=args.db
        )

    elif cmd == "prot-ungapped":
        run_proteines_ungapped(
            w=args.w, n=args.n, x_drop=args.__dict__["x_drop"],
            seuil_t=args.T, seuil=args.seuil, taille_database=args.db
        )

    elif cmd == "prot-gapped":
        run_proteines_gapped(
            w=args.w, n=args.n, taille_database=args.db,
            seuil_t=args.T, seuil_proteines_gapped=args.seuil,
            A=args.A,
            x_drop=args.__dict__["x_drop"],
            x_drop_gapped=args.__dict__["x_drop_gapped"],
            Sg=args.Sg
        )

    elif cmd == "test1":
        test_hits_vs_T(
            w=args.w, n=args.n, taille_db=args.db,
            seuils_T=args.__dict__["seuils_T"], A=args.A, delta=args.delta
        )

    elif cmd == "test2":
        # reconstruit un range propre depuis la liste si besoin
        s_list = sorted(set(args.__dict__["S_accept"]))
        results_b1, results_b2 = compare_acceptance(
            w=args.w, n=args.n, taille_db=args.db,
            T_blast1=args.T1, T_blast2=args.T2,
            A=args.A, x_drop=args.__dict__["x_drop"],
            x_drop_gapped=args.__dict__["x_drop_gapped"],
            Sg=args.Sg, S_accept_range=range(s_list[0], s_list[-1]+1, max(1, (s_list[-1]-s_list[0]) // max(len(s_list)-1,1)))
        )
        # Optionnel: print synthèse
        print("Résumé test2:", {"BLAST I (total)": sum(results_b1), "BLAST II (total)": sum(results_b2)})

    elif cmd == "test3":
        test_time_vs_dbsize(
            w=args.w, n=args.n, T1=args.T1, T2=args.T2, A=args.A,
            x_drop=args.__dict__["x_drop"], x_drop_gapped=args.__dict__["x_drop_gapped"],
            Sg=args.Sg, S_accept=args.__dict__["S_accept"],
            db_sizes=args.__dict__["db_sizes"]
        )

    else:
        print("Commande inconnue.", file=sys.stderr)
        sys.exit(1)