# Stage-ETE24-Geometrie-Spectrale
Stage effectué sous la supervision de [Maxime Fortier Bourque](https://dms.umontreal.ca/fr/repertoire-departement/professeurs/portrait/fortierboum) pendant l'été 2024.<br>

## Exécution
Les fichiers `.edp` doivent être exécutés avec une version modifiée du logiciel [FreeFEM++](https://github.com/FreeFem/FreeFem-sources/). Il faut recompiler le logiciel en suivant [ce guide](https://doc.freefem.org/introduction/installation.html) et mettre en commentaires les lignes 973 `ffassert(kkk++ < 10);` et 992 `ffassert(l++ < 10);`. Attention, les liens vers les dépendances additionnelles dans le fichier CMAKE du projet sont pour la plupart brisés. Utiliser [ce dépôt](https://github.com/FreeFem/FreeFEM-3rdparties) pour trouver les fichiers.

## Principales références
Cook, J. (2021, 26 août). [Properties of Eigenvalues on Riemann Surfaces with Large Symmetry Groups.](http://arxiv.org/abs/2108.11825) arXiv.<br>
Kao, C.-Y., Osting, B. et Oudet, E. (2023). [Computational approaches for extremal geometric eigenvalue problems.](https://doi.org/10.1016/bs.hna.2022.08.001) Dans Handbook of Numerical Analysis (vol. 24, p. 377‑406). Elsevier.<br>
Karcher, H. et Weber, M. J. (1999). [The Geometry of Klein’s Riemann Surface.](https://www.semanticscholar.org/paper/The-Geometry-of-Klein%27s-Riemann-Surface-Karcher-Weber/79bea9fe125d8ef75ed545f93ae0e9c9a71fb618)<br>
Ros, A. (2021, 5 mai). [On the first eigenvalue of the laplacian on compact surfaces of genus three.](https://doi.org/10.48550/arXiv.2010.14857) arXiv.
