# TP_ALGO_NUM — Travaux pratiques d'algorithmes numériques

Ce dépôt contient des implémentations pédagogiques en Python pour des méthodes numériques courantes : résolution de systèmes linéaires, décompositions matricielles, interpolation, méthodes itératives et résolution d'équations différentielles ordinaires.

## Contexte
- **But :** fournir des codes simples, lisibles et utilisables en TD/TP, avec des variantes pure‑Python (pour la pédagogie) et `numpy` (pour la performance).
- **Langage :** Python 3 (tests effectués sous Python 3.13).

## Structure du dépôt
- `TPMTH300/` : code des travaux et exemples.
  - `TPMTH300/EquationLineaire/` : méthodes directes et itératives (Crout, Doolittle, Gauss, Cholesky, pivot partiel/total, etc.).
  - `TPMTH300/INTERPOLATION/` : Lagrange, Newton, moindres carrés.
  - `TPMTH300/EQUADIFF/` : Euler explicite, Runge‑Kutta 2/4.
  - `TPMTH300/ITERACTIVES/` : Jacobi, Gauss‑Seidel.

## Ce qui a été fait (récapitulatif)
- Ajout d'un `.gitignore`.
- Ajout de blocs `if __name__ == '__main__'` pour des tests rapides dans plusieurs modules (Gauss, pivot partiel/total, Gauss‑Jordan, Doolittle, etc.).
- Harmonisation partielle des API (certaines fonctions acceptent `numpy.ndarray`, d'autres `list`).

---

## Détails mathématiques et implantation (méthode par méthode)

Ci‑dessous, pour chaque méthode implémentée, tu trouveras : la formule mathématique essentielle, l'idée algorithmique et les choix concrets d'implémentation en Python.

### 1) Élimination de Gauss (sans pivot)

- Principe mathématique : transformer la matrice $A$ en une matrice triangulaire supérieure $U$ par opérations élémentaires de ligne, ce qui produit un second membre modifié $\tilde b$ tel que $U x = \tilde b$.
- Étapes algorithmiques :
  1. Pour chaque colonne k de 0 à n−2 :
     - prendre pivot $p=A_{k,k}$ ;
     - pour chaque ligne i>k, calculer factor = $A_{i,k}/p$ puis mettre à jour la ligne i : $A_{i,j} \leftarrow A_{i,j} - factor\cdot A_{k,j}$ pour j≥k ; $b_i \leftarrow b_i - factor\cdot b_k$.
  2. Substitution arrière pour calculer x à partir de U et $\tilde b$.
- Implémentation en Python :
  - Variante pédagogique : listes de listes pour `A` et liste pour `b`, boucles `for` imbriquées, copies de travail pour ne pas modifier les entrées originales.
  - Vérification de pivot nul (raise ValueError). Complexité O(n^3).

### 2) Gauss avec pivot partiel

- Principe : à chaque étape, choisir la ligne dont l'élément dans la colonne k a la plus grande valeur absolue ; échanger les lignes pour limiter les erreurs numériques puis éliminer.
- Implémentation : fonction `PivotPartiel(A,k,n)` qui renvoie l'indice de la meilleure ligne, échange simple de lignes (`A[i], A[p] = A[p], A[i]` pour listes; slicing pour `numpy`).

### 3) Pivot total

- Principe : choisir le coefficient maximal en valeur absolue dans la sous‑matrice restante et permuter lignes et colonnes ; enregistrer la permutation de colonnes pour reconstruire la solution.
- Implémentation : parcours double boucle pour trouver l'indice (i*,j*), swap de colonnes et lignes, vecteur `col_order` pour réordonner la solution finale.

### 4) Gauss‑Jordan

- Principe : normaliser la ligne pivot et éliminer dans toutes les autres lignes, obtenant directement la solution sur la colonne augmentée.
- Implémentation : construire la matrice augmentée `M = [A|b]`, normaliser ligne pivot `M[k,k:] /= M[k,k]`, éliminer `M[i,k:] -= M[i,k] * M[k,k:]` pour i≠k.

### 5) Décomposition LU — Doolittle

- Formules clefs :
  - $U_{k,j} = A_{k,j} - \sum_{r=0}^{k-1} L_{k,r} U_{r,j}$ (pour $j\ge k$),
  - $L_{i,k} = \dfrac{1}{U_{k,k}}\left(A_{i,k} - \sum_{r=0}^{k-1} L_{i,r} U_{r,k}\right)$ (pour $i>k$).
- Implémentation : initialiser `L = I`, `U = zeros`, boucles k→ calcul U puis L. Résolution via substitution avant/arrière.

### 6) Décomposition LU — Crout

- Formules clefs (Crout) :
  - $L_{i,k} = A_{i,k} - \sum_{r=0}^{k-1} L_{i,r} U_{r,k}$ pour $i\ge k$,
  - $U_{k,j} = \dfrac{1}{L_{k,k}}\left(A_{k,j} - \sum_{r=0}^{k-1} L_{k,r} U_{r,j}\right)$ pour $j>k$.
- Implémentation : `U` avec diagonale unité, `L` stocke les pivots ; option pivot partiel disponible dans la version pédagogique.

### 7) Cholesky

- Principe : pour $A$ symétrique définie positive, $A = L L^T$ avec
  - $L_{k,k} = \sqrt{A_{k,k} - \sum_{r=0}^{k-1} L_{k,r}^2}$,
  - $L_{i,k} = \dfrac{1}{L_{k,k}}\left(A_{i,k} - \sum_{r=0}^{k-1} L_{i,r}L_{k,r}\right)$ pour $i>k$.
- Implémentation : calcul en place, vérifier symétrie et positivité définie (test sur $L_{k,k}^2$ non négatif). Plus efficace qu'une LU générale (≈ moitié du coût).

### 8) Méthodes itératives (Jacobi, Gauss‑Seidel)

- Formules : décomposer $A=D+L+U$ (D diag),
  - Jacobi : $x^{(k+1)} = D^{-1}(b - (L+U)x^{(k)})$,
  - Gauss‑Seidel : $x^{(k+1)} = (D+L)^{-1}(b - U x^{(k)})$.
- Implémentation : boucle d'itérations, calculer norme du résidu pour la convergence, tolérance et nombre maximal d'itérations en paramètres.

### 9) Interpolation (Lagrange, Newton)

- Lagrange : $P(x)=\sum_i y_i \ell_i(x)$ avec $\ell_i(x)=\prod_{j\ne i} \dfrac{x-x_j}{x_i-x_j}$ ; implémenté en évaluant les polynômes de base pour chaque point d'évaluation.
- Newton : calcul des différences divisées (tableau triangulaire), coefficients pour la forme de Newton, évaluation par Horner généralisé.

### 10) Moindres carrés (régression polynomiale)

- Résoudre $\min_c ||V c - y||_2$ où $V$ est la matrice de Vandermonde ; équations normales $V^T V c = V^T y$.
- Implémentation : construction de `V` puis résolution via décomposition (LU) ou méthode plus stable (QR) si nécessaire.

### 11) Équations différentielles (Euler, RK2, RK4)

- Euler explicite : $y_{n+1}=y_n+h f(t_n,y_n)$ — implémenté dans `EQUADIFF/euler_explicite.py` ; attention stabilité pour pas trop grands.
- RK2 (point milieu / Heun) : deux évaluations de `f` par pas ; RK4 : quatre évaluations et combinaison pondérée (ordre 4).

## Notes transversales d'implémentation

- Types acceptés : certaines fonctions acceptent `numpy.ndarray`, d'autres `list` (les versions pure‑Python vérifient les types). Si besoin, je peux homogénéiser pour accepter les deux (`np.asarray` ou conversion interne en listes).
- Gestion des erreurs : l'implémentation lève des exceptions claires (ValueError) quand elle détecte des cas singuliers (pivot nul, matrice non symétrique pour Cholesky, etc.).
- Performances : les versions `numpy` sont vectorisées et plus rapides ; les versions pédagogiques privilégient la lisibilité (boucles explicites).

---

Si tu veux, je peux :
- extraire ces descriptions en fichiers `docs/` (un fichier par méthode),
- créer des notebooks pédagogiques montrant étapes numériques pas à pas,
- committer ces changements sur la branche `poly`.

---

**Licence** : voir `LICENSE`.

---

## Cours d'immersion détaillé

Cette section est conçue comme un mini‑cours pour chaque méthode numérique implémentée dans le dépôt. Chaque sous‑section présente :
- le contexte et l'objectif ;
- les équations mathématiques essentielles (formalisme) ;
- l'algorithme étape par étape (pseudocode) ;
- remarques d'implémentation et points numériques importants.

Remarque : les équations sont écrites en notation LaTeX (KaTeX) pour une lecture mathématique claire.

### A. Résolution directe : élimination de Gauss

Contexte : résoudre $A x = b$ avec $A\in\mathbb R^{n\times n}$ non singulière.

Formulation : appliquer des opérations élémentaires de lignes pour obtenir une matrice triangulaire supérieure $U$ et un vecteur $\tilde b$ tel que
$$U x = \tilde b.$$

Pseudocode (Gauss sans pivot) :

```
Input: A (n×n), b (n)
for k in 0..n-2:
  p = A[k,k]
  if |p| < eps: raise PivotError
  for i in k+1..n-1:
    factor = A[i,k] / p
    for j in k..n-1:
      A[i,j] -= factor * A[k,j]
    b[i] -= factor * b[k]
# substitution arrière
x = zeros(n)
for i in n-1..0:
  s = sum(A[i,j]*x[j] for j in i+1..n-1)
  x[i] = (b[i] - s) / A[i,i]
return x
```

Remarques : complexité O(n^3). Sans pivot, la méthode peut échouer numériquement si un pivot est proche de zéro.

### B. Pivot partiel et pivot total

Pivot partiel : à l'étape k, choisir la ligne p≥k telle que $|A_{p,k}|$ est maximal, puis échanger les lignes k et p avant l'élimination. Ceci limite l'amplification d'erreurs.

Pivot total : chercher $(p,q)$ qui maximise $|A_{p,q}|$ dans la sous‑matrice restante et échanger lignes et colonnes. Il faut mémoriser la permutation des colonnes pour reconstruire la solution.

Pseudocode (Pivot partiel) :

```
for k in 0..n-2:
  p = argmax_{i=k..n-1} |A[i,k]|
  swap_rows(A,k,p); swap(b,k,p)
  pivot = A[k,k]
  ... (élimination comme Gauss) ...
```

Pseudocode (Pivot total) :

```
col_order = [0..n-1]
for k in 0..n-2:
  (p,q) = argmax_{i=k..n-1, j=k..n-1} |A[i,j]|
  swap_rows(A,k,p); swap(b,k,p)
  swap_cols(A,k,q); swap(col_order[k], col_order[q])
  pivot = A[k,k]
  ... (élimination) ...
# après substitution arrière, permuter x selon col_order inverse
```

Remarque : pivot total est plus coûteux mais utile pour matrices pathologiques.

### C. Décomposition LU (Doolittle et Crout)

Objectif : factoriser $A = L U$.

Doolittle (L avec diagonale unité) :
$$U_{k,j} = A_{k,j} - \sum_{r=0}^{k-1} L_{k,r} U_{r,j}, \quad j\ge k$$
$$L_{i,k} = \frac{1}{U_{k,k}}\left(A_{i,k} - \sum_{r=0}^{k-1} L_{i,r} U_{r,k}\right), \quad i>k$$

Pseudocode (Doolittle) :

```
L = I; U = zeros(n,n)
for k in 0..n-1:
  for j in k..n-1:
    U[k,j] = A[k,j] - sum(L[k,r]*U[r,j] for r in 0..k-1)
  for i in k+1..n-1:
    L[i,k] = (A[i,k] - sum(L[i,r]*U[r,k] for r in 0..k-1)) / U[k,k]
```

Crout (U avec diagonale unité) :
$$L_{i,k} = A_{i,k} - \sum_{r=0}^{k-1} L_{i,r} U_{r,k}, \quad i\ge k$$
$$U_{k,j} = \frac{1}{L_{k,k}}\left(A_{k,j} - \sum_{r=0}^{k-1} L_{k,r} U_{r,j}\right), \quad j>k$$

Remarques d'implémentation :
- Faire attention aux divisions par zéro (tester $U_{k,k}$ ou $L_{k,k}$).\
- Après décomposition, résoudre $L y = b$ puis $U x = y$ (substitution avant/arrière).

### D. Cholesky

Applicable si $A$ symétrique définie positive. $A = L L^T$ avec :
$$L_{k,k} = \sqrt{A_{k,k} - \sum_{r=0}^{k-1} L_{k,r}^2},\qquad L_{i,k} = \frac{1}{L_{k,k}}\left(A_{i,k} - \sum_{r=0}^{k-1} L_{i,r}L_{k,r}\right).$$

Pseudocode :

```
L = zeros(n,n)
for k in 0..n-1:
  sum_sq = sum(L[k,r]**2 for r in 0..k-1)
  L[k,k] = sqrt(A[k,k] - sum_sq)
  for i in k+1..n-1:
    s = sum(L[i,r]*L[k,r] for r in 0..k-1)
    L[i,k] = (A[i,k] - s) / L[k,k]
```

Notes : moitié du coût d'une LU générale en opérations; vérifier la positivité définie via $A$ ou via étapes (si racine négative, A n'est pas SPD).

### E. Méthodes itératives : Jacobi et Gauss‑Seidel

Formulation : $A = D + L + U$.

Jacobi :
$$x^{(k+1)} = D^{-1}(b - (L+U)x^{(k)}).$$

Gauss‑Seidel :
$$x^{(k+1)} = (D+L)^{-1}(b - U x^{(k)}),$$
où on utilise les composantes mises à jour immédiatement.

Pseudocode (Jacobi) :

```
x_new = zeros(n)
while not converged:
  for i in 0..n-1:
    s = sum(A[i,j]*x_old[j] for j in 0..i-1) + sum(A[i,j]*x_old[j] for j in i+1..n-1)
    x_new[i] = (b[i] - s) / A[i,i]
  check convergence; x_old = x_new
```

Pseudocode (Gauss‑Seidel) : on remplace x_old par x (mise à jour en place) et on utilise les valeurs récentes dans la boucle.

Convergence : garanties sous conditions (diagonale dominante, ou matrice symétrique définie positive pour certaines variantes).

### F. Interpolation — Lagrange et Newton (différences divisées)

Lagrange :
$$P(x) = \sum_{i=0}^{n} y_i \ell_i(x),\qquad \ell_i(x)=\prod_{j\ne i} \frac{x-x_j}{x_i-x_j}.$$ 

Newton : construire le tableau des différences divisées $f[x_0], f[x_0,x_1], ...$ et écrire
$$P(x)=a_0+a_1(x-x_0)+a_2(x-x_0)(x-x_1)+\cdots$$
avec coefficients $a_k$ les premières entrées de chaque ligne du tableau.

Implémentation : Lagrange est direct mais coûteux si on évalue pour beaucoup de points ; Newton est plus pratique pour ajouter des points.

### G. Moindres carrés (régression polynomiale)

Construire la matrice de Vandermonde $V$ telle que $V_{i,j}=x_i^j$. Résoudre :
$$V^T V c = V^T y$$
ou préférer une factorisation QR de $V$ pour une meilleure stabilité numérique.

Pseudocode (équations normales) :

```
V = vandermonde(x, deg)
G = V.T @ V
h = V.T @ y
c = solve(G, h)  # par LU ou autre
```

### H. EDO — Euler, RK2, RK4

Euler explicite :
$$y_{n+1} = y_n + h f(t_n,y_n).$$

RK2 (point milieu) :
$$k_1 = f(t_n,y_n),\quad k_2 = f(t_n+h/2, y_n + h k_1/2),\quad y_{n+1} = y_n + h k_2.$$ 

RK4 :
$$k_1=f(t_n,y_n),\ k_2=f(t_n+h/2,y_n+h k_1/2),\ k_3=f(t_n+h/2,y_n+h k_2/2),\ k_4=f(t_n+h,y_n+h k_3),$$
$$y_{n+1}=y_n + \frac{h}{6}(k_1 + 2k_2 + 2k_3 + k_4).$$

Implémentation : stocker t et y sur la grille, boucle de pas, évaluer f aux temps requis.
