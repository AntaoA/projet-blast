Pour utiliser ce projet, vous devez d'abord installer uv sur votre machine. Ensuite, vous pouvez cloner ce dépôt et utiliser les commandes uv pour gérer l'environnement et les dépendances.


git clone https://github.com/AntaoA/projet-blast
cd projet-blast
uv sync



On peut lancer des tests (paramètres non nécessaires et modifiables) :
Test 1 (hits vs T) :
uv run blast.py test1 --w 3 --n 250 --db 300 --A 40 --delta 3 --seed 1337

Test 2 (S_accept) :
uv run blast.py test2 --T1 13 --T2 11 --A 30 --x-drop 16 --x-drop-gapped 40 --Sg 42 --S-accept 60 62 64 66 68 70 --db 2000

Test 3 (temps vs DB) :
uv run blast.py test3 --db-sizes 500 1000 2000 3000 --seed 1


Ou tester différentes implémenations de blast :
Lancer ungapped ADN :
python blast.py adn --w 13 --n 250 --x-drop 20 --seuil 50 --db 1000 --seed 7

Lancer ungapped protéines :
python blast.py prot-ungapped --T 13 --seuil 80 --db 1000

Lancer gapped protéines (BLAST II) :
python blast.py prot-gapped --T 11 --A 40 --Sg 42 --x-drop 20 --x-drop-gapped 40 --seuil 80 --db 1000