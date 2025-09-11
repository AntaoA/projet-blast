# Projet BLAST

Ce projet implémente différentes variantes de l’algorithme **BLAST** (Basic Local Alignment Search Tool), sur ADN et protéines, en versions *ungapped* et *gapped*.

---

## Installation

Pour utiliser ce projet, vous devez d'abord installer **[uv](https://github.com/astral-sh/uv)** sur votre machine.  
Ensuite, vous pouvez cloner ce dépôt et installer les dépendances avec :

```bash
git clone https://github.com/AntaoA/projet-blast
cd projet-blast
uv sync
```

## Run

### On peut lancer des tests (paramètres non nécessaires et modifiables) :

**Test 1 – Hits vs T**
```bash
uv run blast.py test1 --w 3 --n 250 --db 300 --A 40 --delta 3
```

**Test 2 - S_accept**
```bash
uv run blast.py test2 --T1 13 --T2 11 --A 30 --x-drop 16 --x-drop-gapped 40 --Sg 42 --S-accept 60 62 64 66 68 70 --db 2000
```

**Test 3 - temps vs DB** :
```bash
uv run blast.py test3 --db-sizes 500 1000 2000 3000
```

### Ou tester différentes implémenations de blast (paramètres non nécessaires et modifiables) :
**ungapped ADN (BLAST I)**
```bash
uv run blast.py adn --w 13 --n 250 --x-drop 20 --seuil 50 --db 1000
```

**ungapped protéines (BLAST I)**
```bash
uv run blast.py prot-ungapped --T 13 --seuil 80 --db 1000
```

**gapped protéines (BLAST II)**
```bash
uv run blast.py prot-gapped --T 11 --A 40 --Sg 42 --x-drop 20 --x-drop-gapped 40 --seuil 80 --db 1000
```