import os
import numpy as np
import random
import math

# Pour éventuellement afficher le facteur conforme
import matplotlib.pyplot as plt


# Fonction approximant un nombre avec un nombre
# prédéfini de chiffres significatifs.
def approx(nbr, nbr_chiffres_significatifs):
    partie_entiere = 0
    partie_decimale = 0.0
    nbr_chiffres_partie_decimale = 0

    decomposition_partie_entiere = []
    decomposition_partie_decimale = []

    resultat = 0

    if nbr > 0:
        partie_entiere = math.floor(nbr)
        partie_decimale = nbr - math.floor(nbr)

        while partie_decimale != 0:
            decomposition_partie_decimale.append(int(math.floor(partie_decimale * 10)))
            partie_decimale = partie_decimale * 10 - decomposition_partie_decimale[nbr_chiffres_partie_decimale]
            nbr_chiffres_partie_decimale += 1

    else:
        partie_entiere = math.ceil(nbr)
        partie_decimale = nbr - math.ceil(nbr)

        while partie_decimale != 0:
            decomposition_partie_decimale.append(int(math.ceil(partie_decimale * 10)))
            partie_decimale = partie_decimale * 10 - decomposition_partie_decimale[nbr_chiffres_partie_decimale]
            nbr_chiffres_partie_decimale += 1

    while partie_entiere != 0:
        decomposition_partie_entiere.append(partie_entiere % 10)
        partie_entiere = (partie_entiere - partie_entiere % 10) // 10

    decomposition_partie_entiere.reverse()

    if len(decomposition_partie_entiere) == 0:
        position = 0

        while decomposition_partie_decimale[position] == 0:
            position += 1

        for i in range(position, position + min(len(decomposition_partie_decimale), nbr_chiffres_significatifs)):
            if i == position + nbr_chiffres_significatifs - 1:
                decomposition_partie_decimale[i] = round(float(decomposition_partie_decimale[i]) + float(decomposition_partie_decimale[i + 1]) / 10)

            resultat += decomposition_partie_decimale[i] * pow(10, - i - 1)
    else:
        for i in range(min(len(decomposition_partie_entiere), nbr_chiffres_significatifs)):
            if i == nbr_chiffres_significatifs - 1 and len(decomposition_partie_entiere) != nbr_chiffres_significatifs:
                decomposition_partie_entiere[-i - 1] = round(float(decomposition_partie_entiere[-i - 1]) + float(decomposition_partie_entiere[-i - 2]) / 10)

            resultat += decomposition_partie_entiere[-i - 1] * pow(10, len(decomposition_partie_entiere) - i - 1)

        for i in range(min(len(decomposition_partie_decimale), nbr_chiffres_significatifs - len(decomposition_partie_entiere))):
            if i + min(len(decomposition_partie_entiere), nbr_chiffres_significatifs) == nbr_chiffres_significatifs - 1:
                decomposition_partie_decimale[i] = round(float(decomposition_partie_decimale[i]) + float(decomposition_partie_decimale[i + 1]) / 10)

            resultat += decomposition_partie_decimale[i] * pow(10, - i - 1)

    return resultat


# Nom du fichier .edp
nom = "14triangles"

# Chemin du fichier .edp
chemin = "path"

# Chemin de l'exécutable modifié, si besoin, de FreeFem++
cheminFF = "path"

n = 14

taille = 1
espacement = 2

a = np.pi / 6 + 0.1

# Nombre de gaussiennes en x / y.
divX = 20
divY = 20

# Liste des coefficients des gaussiennes.
coeff = []

try:
    # Possibilité de charger les coefficients depuis un fichier.
    fichier = open("path", "r")
    texte = fichier.read()

    if len(texte) != 0:
        for elem in texte.split("\n"):
            coeff.append(float(elem) / 2)

    fichier.close()

except:
    print("Le fichier n'existe pas")

    coeff = [random.random() for i in range((divX + 1) * (divY + 1))]

def genererFacteurConforme(divX, divY, coeff):
    """
    Fonction génèrant le facteur conforme pour le fichier .edp
    """

    pos = 0

    facteurConforme = ""

    for i in range(divX + 1):
        for j in range(divY + 1):
            facteurConforme += f"{coeff[pos]} * exp(- ( (x - taille * {i / divX} - floor(x))^2 + (y - taille * {j / divY})^2  ) / 0.00001)"
            
            if pos != (divX + 1) * (divY + 1) - 1:
                facteurConforme += " + "

            pos += 1

    return facteurConforme

def valeurPropre(coeff):
    """
    Fonction générant le fichier .edp puis l'exécutant.
    La fonction renvoie un tuple valeur propre normalisée / gradient.
    """

    # Une formule pour calculer le gradient est présente sur cette ligne.
    # Voir arXiv:2008.03385v2.
    texte = f'int divX = {divX};\nint divY = {divY};' + '\n\nreal taille = 0.5;\nreal espacement = 1;\nreal alpha = 1.0455119243681479;\nreal beta = 1.0505688048534976;\n\nint n = 25;\nint nbrVP = 10;\n\nreal t1 = 0;\nreal pente = 0;\n\nif (beta > pi / 2)\n{\n    t1 = - taille * tan(pi - beta) / (tan(alpha) - tan(pi - beta));\n    pente = tan(pi - beta);\n    n *= -1;\n}\n\nelse\n{\n    t1 = taille * tan(beta) / (tan(alpha) + tan(beta));\n    pente = - tan(beta);\n}\n\nborder T1C1 (t = espacement * 0, espacement * 0 + taille) {x = t; y = 0; label = 1;};\nborder T1C2 (t = espacement * 0, espacement * 0 + t1) {x = t; y = tan(alpha) * (t - 0 * espacement); label = 2;};\nborder T1C3 (t = 0 * espacement + taille, 0 * espacement + t1) {x = t; y = pente * (t - espacement * 0 - taille); label = 3;};\nborder T2C1 (t = espacement * 1, espacement * 1 + taille) {x = t; y = 0; label = 4;};\nborder T2C2 (t = espacement * 1, espacement * 1 + t1) {x = t; y = tan(alpha) * (t - 1 * espacement); label = 5;};\nborder T2C3 (t = 1 * espacement + taille, 1 * espacement + t1) {x = t; y = pente * (t - espacement * 1 - taille); label = 6;};\nborder T3C1 (t = espacement * 2, espacement * 2 + taille) {x = t; y = 0; label = 7;};\nborder T3C2 (t = espacement * 2, espacement * 2 + t1) {x = t; y = tan(alpha) * (t - 2 * espacement); label = 8;};\nborder T3C3 (t = 2 * espacement + taille, 2 * espacement + t1) {x = t; y = pente * (t - espacement * 2 - taille); label = 9;};\nborder T4C1 (t = espacement * 3, espacement * 3 + taille) {x = t; y = 0; label = 10;};\nborder T4C2 (t = espacement * 3, espacement * 3 + t1) {x = t; y = tan(alpha) * (t - 3 * espacement); label = 11;};\nborder T4C3 (t = 3 * espacement + taille, 3 * espacement + t1) {x = t; y = pente * (t - espacement * 3 - taille); label = 12;};\nborder T5C1 (t = espacement * 4, espacement * 4 + taille) {x = t; y = 0; label = 13;};\nborder T5C2 (t = espacement * 4, espacement * 4 + t1) {x = t; y = tan(alpha) * (t - 4 * espacement); label = 14;};\nborder T5C3 (t = 4 * espacement + taille, 4 * espacement + t1) {x = t; y = pente * (t - espacement * 4 - taille); label = 15;};\nborder T6C1 (t = espacement * 5, espacement * 5 + taille) {x = t; y = 0; label = 16;};\nborder T6C2 (t = espacement * 5, espacement * 5 + t1) {x = t; y = tan(alpha) * (t - 5 * espacement); label = 17;};\nborder T6C3 (t = 5 * espacement + taille, 5 * espacement + t1) {x = t; y = pente * (t - espacement * 5 - taille); label = 18;};\nborder T7C1 (t = espacement * 6, espacement * 6 + taille) {x = t; y = 0; label = 19;};\nborder T7C2 (t = espacement * 6, espacement * 6 + t1) {x = t; y = tan(alpha) * (t - 6 * espacement); label = 20;};\nborder T7C3 (t = 6 * espacement + taille, 6 * espacement + t1) {x = t; y = pente * (t - espacement * 6 - taille); label = 21;};\nborder T8C1 (t = espacement * 7, espacement * 7 + taille) {x = t; y = 0; label = 22;};\nborder T8C2 (t = espacement * 7, espacement * 7 + t1) {x = t; y = tan(alpha) * (t - 7 * espacement); label = 23;};\nborder T8C3 (t = 7 * espacement + taille, 7 * espacement + t1) {x = t; y = pente * (t - espacement * 7 - taille); label = 24;};\nborder T9C1 (t = espacement * 8, espacement * 8 + taille) {x = t; y = 0; label = 25;};\nborder T9C2 (t = espacement * 8, espacement * 8 + t1) {x = t; y = tan(alpha) * (t - 8 * espacement); label = 26;};\nborder T9C3 (t = 8 * espacement + taille, 8 * espacement + t1) {x = t; y = pente * (t - espacement * 8 - taille); label = 27;};\nborder T10C1 (t = espacement * 9, espacement * 9 + taille) {x = t; y = 0; label = 28;};\nborder T10C2 (t = espacement * 9, espacement * 9 + t1) {x = t; y = tan(alpha) * (t - 9 * espacement); label = 29;};\nborder T10C3 (t = 9 * espacement + taille, 9 * espacement + t1) {x = t; y = pente * (t - espacement * 9 - taille); label = 30;};\nborder T11C1 (t = espacement * 10, espacement * 10 + taille) {x = t; y = 0; label = 31;};\nborder T11C2 (t = espacement * 10, espacement * 10 + t1) {x = t; y = tan(alpha) * (t - 10 * espacement); label = 32;};\nborder T11C3 (t = 10 * espacement + taille, 10 * espacement + t1) {x = t; y = pente * (t - espacement * 10 - taille); label = 33;};\nborder T12C1 (t = espacement * 11, espacement * 11 + taille) {x = t; y = 0; label = 34;};\nborder T12C2 (t = espacement * 11, espacement * 11 + t1) {x = t; y = tan(alpha) * (t - 11 * espacement); label = 35;};\nborder T12C3 (t = 11 * espacement + taille, 11 * espacement + t1) {x = t; y = pente * (t - espacement * 11 - taille); label = 36;};\nborder T13C1 (t = espacement * 12, espacement * 12 + taille) {x = t; y = 0; label = 37;};\nborder T13C2 (t = espacement * 12, espacement * 12 + t1) {x = t; y = tan(alpha) * (t - 12 * espacement); label = 38;};\nborder T13C3 (t = 12 * espacement + taille, 12 * espacement + t1) {x = t; y = pente * (t - espacement * 12 - taille); label = 39;};\nborder T14C1 (t = espacement * 13, espacement * 13 + taille) {x = t; y = 0; label = 40;};\nborder T14C2 (t = espacement * 13, espacement * 13 + t1) {x = t; y = tan(alpha) * (t - 13 * espacement); label = 41;};\nborder T14C3 (t = 13 * espacement + taille, 13 * espacement + t1) {x = t; y = pente * (t - espacement * 13 - taille); label = 42;};\n\n\nmesh Th = buildmesh(T1C1(n) + T1C2(-n) + T1C3(n) + T2C1(n) + T2C2(-n) + T2C3(n) + T3C1(n) + T3C2(-n) + T3C3(n) + T4C1(n) + T4C2(-n) + T4C3(n) + T5C1(n) + T5C2(-n) + T5C3(n) + T6C1(n) + T6C2(-n) + T6C3(n) + T7C1(n) + T7C2(-n) + T7C3(n) + T8C1(n) + T8C2(-n) + T8C3(n) + T9C1(n) + T9C2(-n) + T9C3(n) + T10C1(n) + T10C2(-n) + T10C3(n) + T11C1(n) + T11C2(-n) + T11C3(n) + T12C1(n) + T12C2(-n) + T12C3(n) + T13C1(n) + T13C2(-n) + T13C3(n) + T14C1(n) + T14C2(-n) + T14C3(n));\nfespace Vh (Th, P2, periodic = [[1, x + espacement], [4, x],\n[5, y], [8, y],\n[6, y], [21, y],\n[7, x + espacement], [10, x],\n[11, y], [14, y],\n[12, y], [27, y],\n[13, x + espacement], [16, x],\n[17, y], [20, y],\n[18, y], [33, y],\n[19, x + espacement], [22, x],\n[23, y], [26, y],\n[24, y], [39, y],\n[25, x + espacement], [28, x],\n[29, y], [32, y],\n[30, y], [3, y],\n[31, x + espacement], [34, x],\n[35, y], [38, y],\n[36, y], [9, y],\n[37, x + espacement], [40, x],\n[41, y], [2, y],\n[42, y], [15, y]]);\n\nVh u1, u2;\nreal sigma = 0.00001;\n\n//func facteurConforme = 1 / (1 - (x - floor(x))^2 - y^2)^2;\n//func facteurConforme = exp(- ((x - floor(x) - taille / 2)^2 + (y - 0.1)^2));\n\nfunc facteurConforme = 1;\n\n' + f'varf op (u1, u2) = int2d(Th)(dx(u1) * dx(u2) + dy(u1) * dy(u2) - sigma * u1 * u2 * ({genererFacteurConforme(divX, divY, coeff)}));\nvarf bm ([u1], [u2]) = int2d(Th)(u1 * u2 * ({genererFacteurConforme(divX, divY, coeff)}));' + '\n\nmatrix OP = op(Vh, Vh, solver = Crout, factorize = 1);\nmatrix B = bm(Vh, Vh, solver = CG, eps = 1e-20);\n\nreal[int] ev(nbrVP);\nVh[int] eV(nbrVP);\n\nint k = EigenValue(OP, B, sym = true, sigma = sigma, value = ev, vector = eV, tol=1e-10, maxit = 0, ncv = 0);\nreal val = 0;' + f'\nreal aire = int2d(Th)({genererFacteurConforme(divX, divY, coeff)});' + '\nint nbr = 0;\nint pos = 0;\n\nfor (int i = 0; i < k; i++)\n{\n    if (abs(ev[i]) > 0.01)\n    {\n        val = ev[i] * aire;\n        nbr = i;\n        break;\n    }\n}\n\ncout << "Retour : " << val << "|" << endl;\n\ncout << "GRADIENT [" << endl;\n\n' + f'real denom = int2d(Th)(({genererFacteurConforme(divX, divY, coeff)}) * eV[nbr]^2);' + '\n\nfor (int i = 0; i < divX + 1; i++)\n{\n    for (int j = 0; j < divY + 1; j++)\n    {\n        cout << -ev[nbr] * int2d(Th)((aire * eV[nbr]^2 / (denom) - 1) * exp(- ((x - taille * i / divX - floor(x))^2 + (y - taille * j / divY + 1)^2 ) / 0.00001)) << endl;\n        pos++;\n    }\n}\n\ncout << "]" << endl;\ncout << "Aire : " << aire << endl;'

    fichier = open(chemin + nom + ".edp", "w")
    fichier.write(texte)
    fichier.close()

    if cheminFF != "":
        texte = os.popen(f"{cheminFF} {chemin}{nom}.edp").read()
    else:
        texte = os.popen(f"{chemin}{nom}.edp").read()

    if len(texte.split("Exec error ")) == 1:
        val = float(texte.split("Retour : ")[-1].split("|")[0])
        aire = float(texte.split("Aire : ")[-1].split("\n")[0])

        print(f"vpn = {val}, vp = {val / aire}, aire = {aire}")

        return (val, [float(texte.split("GRADIENT [")[-1].split("]")[0].split("\n")[i]) for i in range((divX + 1) * (divY + 1) + 1) if texte.split("GRADIENT [")[-1].split("]")[0].split("\n")[i] != ''])

    return (0, [0 for i in range((divX + 1) * (divY + 1) + 1)])


def vp(coeff):
    """
    Recherche linéaire :
    Fonction générant le fichier .edp p l'exécutant.
    La fonction renvoie la valeur propre normalisée.
    """

    texte = f'int divX = {divX};\nint divY = {divY};' + '\n\nreal taille = 0.5;\nreal espacement = 1;\nreal alpha = 1.0455119243681479;\nreal beta = 1.0505688048534976;\n\nint n = 25;\nint nbrVP = 10;\n\nreal t1 = 0;\nreal pente = 0;\n\nif (beta > pi / 2)\n{\n    t1 = - taille * tan(pi - beta) / (tan(alpha) - tan(pi - beta));\n    pente = tan(pi - beta);\n    n *= -1;\n}\n\nelse\n{\n    t1 = taille * tan(beta) / (tan(alpha) + tan(beta));\n    pente = - tan(beta);\n}\n\nborder T1C1 (t = espacement * 0, espacement * 0 + taille) {x = t; y = 0; label = 1;};\nborder T1C2 (t = espacement * 0, espacement * 0 + t1) {x = t; y = tan(alpha) * (t - 0 * espacement); label = 2;};\nborder T1C3 (t = 0 * espacement + taille, 0 * espacement + t1) {x = t; y = pente * (t - espacement * 0 - taille); label = 3;};\nborder T2C1 (t = espacement * 1, espacement * 1 + taille) {x = t; y = 0; label = 4;};\nborder T2C2 (t = espacement * 1, espacement * 1 + t1) {x = t; y = tan(alpha) * (t - 1 * espacement); label = 5;};\nborder T2C3 (t = 1 * espacement + taille, 1 * espacement + t1) {x = t; y = pente * (t - espacement * 1 - taille); label = 6;};\nborder T3C1 (t = espacement * 2, espacement * 2 + taille) {x = t; y = 0; label = 7;};\nborder T3C2 (t = espacement * 2, espacement * 2 + t1) {x = t; y = tan(alpha) * (t - 2 * espacement); label = 8;};\nborder T3C3 (t = 2 * espacement + taille, 2 * espacement + t1) {x = t; y = pente * (t - espacement * 2 - taille); label = 9;};\nborder T4C1 (t = espacement * 3, espacement * 3 + taille) {x = t; y = 0; label = 10;};\nborder T4C2 (t = espacement * 3, espacement * 3 + t1) {x = t; y = tan(alpha) * (t - 3 * espacement); label = 11;};\nborder T4C3 (t = 3 * espacement + taille, 3 * espacement + t1) {x = t; y = pente * (t - espacement * 3 - taille); label = 12;};\nborder T5C1 (t = espacement * 4, espacement * 4 + taille) {x = t; y = 0; label = 13;};\nborder T5C2 (t = espacement * 4, espacement * 4 + t1) {x = t; y = tan(alpha) * (t - 4 * espacement); label = 14;};\nborder T5C3 (t = 4 * espacement + taille, 4 * espacement + t1) {x = t; y = pente * (t - espacement * 4 - taille); label = 15;};\nborder T6C1 (t = espacement * 5, espacement * 5 + taille) {x = t; y = 0; label = 16;};\nborder T6C2 (t = espacement * 5, espacement * 5 + t1) {x = t; y = tan(alpha) * (t - 5 * espacement); label = 17;};\nborder T6C3 (t = 5 * espacement + taille, 5 * espacement + t1) {x = t; y = pente * (t - espacement * 5 - taille); label = 18;};\nborder T7C1 (t = espacement * 6, espacement * 6 + taille) {x = t; y = 0; label = 19;};\nborder T7C2 (t = espacement * 6, espacement * 6 + t1) {x = t; y = tan(alpha) * (t - 6 * espacement); label = 20;};\nborder T7C3 (t = 6 * espacement + taille, 6 * espacement + t1) {x = t; y = pente * (t - espacement * 6 - taille); label = 21;};\nborder T8C1 (t = espacement * 7, espacement * 7 + taille) {x = t; y = 0; label = 22;};\nborder T8C2 (t = espacement * 7, espacement * 7 + t1) {x = t; y = tan(alpha) * (t - 7 * espacement); label = 23;};\nborder T8C3 (t = 7 * espacement + taille, 7 * espacement + t1) {x = t; y = pente * (t - espacement * 7 - taille); label = 24;};\nborder T9C1 (t = espacement * 8, espacement * 8 + taille) {x = t; y = 0; label = 25;};\nborder T9C2 (t = espacement * 8, espacement * 8 + t1) {x = t; y = tan(alpha) * (t - 8 * espacement); label = 26;};\nborder T9C3 (t = 8 * espacement + taille, 8 * espacement + t1) {x = t; y = pente * (t - espacement * 8 - taille); label = 27;};\nborder T10C1 (t = espacement * 9, espacement * 9 + taille) {x = t; y = 0; label = 28;};\nborder T10C2 (t = espacement * 9, espacement * 9 + t1) {x = t; y = tan(alpha) * (t - 9 * espacement); label = 29;};\nborder T10C3 (t = 9 * espacement + taille, 9 * espacement + t1) {x = t; y = pente * (t - espacement * 9 - taille); label = 30;};\nborder T11C1 (t = espacement * 10, espacement * 10 + taille) {x = t; y = 0; label = 31;};\nborder T11C2 (t = espacement * 10, espacement * 10 + t1) {x = t; y = tan(alpha) * (t - 10 * espacement); label = 32;};\nborder T11C3 (t = 10 * espacement + taille, 10 * espacement + t1) {x = t; y = pente * (t - espacement * 10 - taille); label = 33;};\nborder T12C1 (t = espacement * 11, espacement * 11 + taille) {x = t; y = 0; label = 34;};\nborder T12C2 (t = espacement * 11, espacement * 11 + t1) {x = t; y = tan(alpha) * (t - 11 * espacement); label = 35;};\nborder T12C3 (t = 11 * espacement + taille, 11 * espacement + t1) {x = t; y = pente * (t - espacement * 11 - taille); label = 36;};\nborder T13C1 (t = espacement * 12, espacement * 12 + taille) {x = t; y = 0; label = 37;};\nborder T13C2 (t = espacement * 12, espacement * 12 + t1) {x = t; y = tan(alpha) * (t - 12 * espacement); label = 38;};\nborder T13C3 (t = 12 * espacement + taille, 12 * espacement + t1) {x = t; y = pente * (t - espacement * 12 - taille); label = 39;};\nborder T14C1 (t = espacement * 13, espacement * 13 + taille) {x = t; y = 0; label = 40;};\nborder T14C2 (t = espacement * 13, espacement * 13 + t1) {x = t; y = tan(alpha) * (t - 13 * espacement); label = 41;};\nborder T14C3 (t = 13 * espacement + taille, 13 * espacement + t1) {x = t; y = pente * (t - espacement * 13 - taille); label = 42;};\n\n\nmesh Th = buildmesh(T1C1(n) + T1C2(-n) + T1C3(n) + T2C1(n) + T2C2(-n) + T2C3(n) + T3C1(n) + T3C2(-n) + T3C3(n) + T4C1(n) + T4C2(-n) + T4C3(n) + T5C1(n) + T5C2(-n) + T5C3(n) + T6C1(n) + T6C2(-n) + T6C3(n) + T7C1(n) + T7C2(-n) + T7C3(n) + T8C1(n) + T8C2(-n) + T8C3(n) + T9C1(n) + T9C2(-n) + T9C3(n) + T10C1(n) + T10C2(-n) + T10C3(n) + T11C1(n) + T11C2(-n) + T11C3(n) + T12C1(n) + T12C2(-n) + T12C3(n) + T13C1(n) + T13C2(-n) + T13C3(n) + T14C1(n) + T14C2(-n) + T14C3(n));\nfespace Vh (Th, P2, periodic = [[1, x + espacement], [4, x],\n[5, y], [8, y],\n[6, y], [21, y],\n[7, x + espacement], [10, x],\n[11, y], [14, y],\n[12, y], [27, y],\n[13, x + espacement], [16, x],\n[17, y], [20, y],\n[18, y], [33, y],\n[19, x + espacement], [22, x],\n[23, y], [26, y],\n[24, y], [39, y],\n[25, x + espacement], [28, x],\n[29, y], [32, y],\n[30, y], [3, y],\n[31, x + espacement], [34, x],\n[35, y], [38, y],\n[36, y], [9, y],\n[37, x + espacement], [40, x],\n[41, y], [2, y],\n[42, y], [15, y]]);\n\nVh u1, u2;\nreal sigma = 0.00001;\n\n//func facteurConforme = 1 / (1 - (x - floor(x))^2 - y^2)^2;\n//func facteurConforme = exp(- ((x - floor(x) - taille / 2)^2 + (y - 0.1)^2));\n\nfunc facteurConforme = 1;\n\n' + f'varf op (u1, u2) = int2d(Th)(dx(u1) * dx(u2) + dy(u1) * dy(u2) - sigma * u1 * u2 * ({genererFacteurConforme(divX, divY, coeff)}));\nvarf bm ([u1], [u2]) = int2d(Th)(u1 * u2 * ({genererFacteurConforme(divX, divY, coeff)}));' + '\n\nmatrix OP = op(Vh, Vh, solver = Crout, factorize = 1);\nmatrix B = bm(Vh, Vh, solver = CG, eps = 1e-20);\n\nreal[int] ev(nbrVP);\nVh[int] eV(nbrVP);\n\nint k = EigenValue(OP, B, sym = true, sigma = sigma, value = ev, vector = eV, tol=1e-10, maxit = 0, ncv = 0);\nreal val = 0;' + f'\nreal aire = int2d(Th)({genererFacteurConforme(divX, divY, coeff)});' + '\nint nbr = 0;\nint pos = 0;\n\nfor (int i = 0; i < k; i++)\n{\n    if (abs(ev[i]) > 0.01)\n    {\n        val = ev[i] * aire;\n        nbr = i;\n        break;\n    }\n}\n\ncout << "Retour : " << val << "|" << endl;\ncout << "Aire : " << aire << endl;'

    fichier = open(chemin + nom + ".edp", "w")
    fichier.write(texte)
    fichier.close()

    if cheminFF != "":
        texte = os.popen(f"{cheminFF} {chemin}{nom}.edp").read()
    else:
        texte = os.popen(f"{chemin}{nom}.edp").read()

    if len(texte.split("Exec error ")) == 1:
        return float(texte.split("Retour : ")[-1].split("|")[0])

    return 0

# Borne de A. Ros. Voir arXiv:2010.14857v3
cible = (16 * (4 - np.sqrt(7)) * np.pi)

temp = []

def rechercheLineaire(x, gradient):
    """
    Fonction effectuant la recherche linéaire selon la direction
    donnée par le gradient. Le but est de trouver le meilleur pas
    (pasOptimal) tel que f(x - pasOptimal * gradient) > f(x).
    La fonction renvoie le pas optimal.
    """

    # Pour éviter de faire trop de calculs, on fixe certains pas
    # selon une échelle de grandeur.
    valPermises = [-1000, -500, -200, -100, -50, -10, -5, -2, -1, -0.1, -0.01, -0.0001, 0.0001, 0.01, 0.1, 1, 2, 5, 10, 50, 100, 200, 500, 1000]
    pasOptimal = []

    resultats = []

    for i in range(len(valPermises)):
        temp = list(x)

        for j in range(len(coeff)):
            if gradient[j] != 0 and approx(abs(gradient[j] * 100), 1) != 0:
                temp[j] = temp[j] - valPermises[i] * gradient[j] / approx(abs(gradient[j] * 100), 1)

        resultats.append(vp(temp))

    iMax = 0
    resultatMax = 0

    for i in range(len(valPermises)):
        if resultatMax < resultats[i]:
            resultatMax = resultats[i]
            iMax = i

    resultats = []

    # Section pour raffiner la selection du pas.
    for i in np.linspace(valPermises[iMax] * 0.75, valPermises[iMax] * 1.25, 10):
        pasOptimal.append(i)

        temp = list(x)

        for j in range(len(coeff)):
            if gradient[j] != 0 and approx(abs(gradient[j] * 100), 1) != 0:
                temp[j] = temp[j] - i * gradient[j] / approx(abs(gradient[j] * 100), 1)

        resultats.append(vp(temp))

    iMax = 0
    resultatMax = 0

    for i in range(len(pasOptimal)):
        if resultatMax < resultats[i]:
            resultatMax = resultats[i]
            iMax = i
    
    return pasOptimal[iMax]

base = (0, [0] * len(coeff))

while abs(cible - base[0]) > 0.01:
    base = valeurPropre(coeff)

    eta = rechercheLineaire(coeff, base[1])

    for i in range(len(coeff)):
        if base[1][i] != 0 and approx(abs(base[1][i] * 100), 1) != 0:
            coeff[i] = coeff[i] - eta * base[1][i] / approx(abs(base[1][i] * 100), 1)

    # Sauvegarde éventuelle des coefficients et de la valeur propre normalisée.
    fichier = open("path", "a")
    fichier.write(str([base[0], coeff]))
    fichier.write("\n")
    fichier.close()