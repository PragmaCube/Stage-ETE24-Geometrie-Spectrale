import numpy as np
import math

import os

import random

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator

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

###########################################################
# Nombre de valeurs propres à calculer
nbrVP = 5

# Nombre de divisions par côté
precision = 25

# Nom du fichier .edp
nom = "klein"

# Chemin du fichier .edp
chemin = "path"

# Chemin de l'exécutable modifié, si besoin, de FreeFem++
cheminFF = "C:/Users/pragm/OneDrive/Bureau/FreeFEM++.exe"

divX = 26
divY = 26

# Il faudra mettre à jour les coefficients d'une manière
# avant de créer le facteur conforme...

coeff = []

try:
    fichier = open("path", "r")
    texte = fichier.read()

    if len(texte) != 0:
        for elem in texte.split("\n"):
            coeff.append(float(elem))

    fichier.close()

except:
    print("Le fichier n'existe pas")

    for i in range(divX + 1):
        for j in range(divY + 1):
            x = -2 * i / divX + 1
            y = -2 * j / divY + 1

            if x ** 2 + y ** 2 == 1:
                coeff.append(0)
                continue

            coeff.append((1 / (1 - (x ** 2 + y ** 2)) ** 2))

def genererFacteurConforme(divX, divY, coeff):
    pos = 0

    facteurConforme = ""

    for i in range(divX + 1):
        for j in range(divY + 1):
            facteurConforme += f"{coeff[pos]} * exp(- ( (x - {2 * i / divX} + 1)^2 + (y - {2 * j / divY} + 1)^2  ) / 0.001)"
            
            if pos != (divX + 1) * (divY + 1) - 1:
                facteurConforme += " + "

            pos += 1

    return facteurConforme

def facteurConforme(x):
    pos = 0

    somme = 0

    for i in range(divX + 1):
        for j in range(divY + 1):
            somme += coeff[pos] * np.exp(- ((x[0] - 2 * i / divX + 1) ** 2 + (x[1] - 2 * j / divY + 1) ** 2) / 0.001)

            pos += 1

    return somme


def valeurPropre(params):
    """
    Fonction générant le fichier .edp puis l'exécutant.
    La fonction renvoie un tuple valeur propre normalisée / gradient.
    """

    # Une formule pour calculer le gradient est présente sur cette ligne.
    # Voir arXiv:2008.03385v2.
    programme = f'int divX = {divX};\nint divY = {divY};' + '\n\nreal angle1 = 0.4487989505128276;\n\nreal cx = 1.031969722;\nreal r = 0.254875473;\n\nreal s = 0.7770942488589633;\n\nreal x06 = s * cos(5 * pi / 7);\nreal y06 = s * sin(5 * pi / 7);\nreal x03 = s * cos(2 * pi / 7);\nreal y03 = s * sin(2 * pi / 7);\nreal x05 = s * cos(4 * pi / 7);\nreal y05 = s * sin(4 * pi / 7);\nreal x10 = s * cos(9 * pi / 7);\nreal y10 = s * sin(9 * pi / 7);\nreal x07 = s * cos(6 * pi / 7);\nreal y07 = s * sin(6 * pi / 7);\nreal x12 = s * cos(11 * pi / 7);\nreal y12 = s * sin(11 * pi / 7);\nreal x09 = s * cos(8 * pi / 7);\nreal y09 = s * sin(8 * pi / 7);\nreal x14 = s * cos(13 * pi / 7);\nreal y14 = s * sin(13 * pi / 7);\nreal x11 = s * cos(10 * pi / 7);\nreal y11 = s * sin(10 * pi / 7);\nreal x02 = s * cos(pi / 7);\nreal y02 = s * sin(pi / 7);\nreal x13 = s * cos(12 * pi / 7);\nreal y13 = s * sin(12 * pi / 7);\nreal x04 = s * cos(3 * pi / 7);\nreal y04 = s * sin(3 * pi / 7);\n\nborder G1 (t = 9 * pi / 7,5 * pi / 7){ x = cx + r * cos(t); y = r * sin(t); label = 1; };\nborder G2 (t = 9 * pi / 7,5 * pi / 7){ x = cos(pi / 7) * (cx + r * cos(t)) - sin(pi / 7) * (r * sin(t));\ny = sin(pi / 7) * (cx + r * cos(t)) + cos(pi / 7) * (r * sin(t)); label = 2; };\nborder G3 (t = 9 * pi / 7,5 * pi / 7){x = cos(2 * pi / 7) * (cx + r * cos(t)) - sin(2 * pi / 7) * r * sin(t);\ny = sin(2 * pi / 7) * (cx + r * cos(t)) + cos(2 * pi / 7) * r * sin(t); label = 3; };\nborder G4 (t = 9 * pi / 7,5 * pi / 7){x = cos(3 * pi / 7) * (cx + r * cos(t)) - sin(3 * pi / 7) * r * sin(t);\ny = sin(3 * pi / 7) * (cx + r * cos(t)) + cos(3 * pi / 7) * r * sin(t); label = 4; };\nborder G5 (t = 9 * pi / 7,5 * pi / 7){ x = cos(4 * pi / 7) * (cx + r * cos(t)) - sin(4 * pi / 7) * r * sin(t);\ny = sin(4 * pi / 7) * (cx + r * cos(t)) + cos(4 * pi / 7) * r * sin(t); label = 5; };\nborder G6(t = 9 * pi / 7,5 * pi / 7){x = cos(5 * pi / 7) * (cx + r * cos(t)) - sin(5 * pi / 7) * r * sin(t);\ny = sin(5 * pi / 7) * (cx + r * cos(t)) + cos(5 * pi / 7) * r * sin(t); label = 6; };\nborder G7(t = 9 * pi / 7,5 * pi / 7){x = cos(6 * pi / 7) * (cx + r * cos(t)) - sin(6 * pi / 7) * r * sin(t);\ny = sin(6 * pi / 7) * (cx + r * cos(t)) + cos(6 * pi / 7) * r * sin(t); label = 7; };\nborder G8(t = 9 * pi / 7,5 * pi / 7){x = cos(7 * pi / 7) * (cx + r * cos(t)) - sin(7 * pi / 7) * r * sin(t);\ny = sin(7 * pi / 7) * (cx + r * cos(t)) + cos(7 * pi / 7) * r * sin(t); label = 8; };\nborder G9(t = 9 * pi / 7,5 * pi / 7){x = cos(8 * pi / 7) * (cx + r * cos(t)) - sin(8 * pi / 7) * r * sin(t);\ny = sin(8 * pi / 7) * (cx + r * cos(t)) + cos(8 * pi / 7) * r * sin(t); label = 9; };\nborder G10(t = 9 * pi / 7,5 * pi / 7){x = cos(9 * pi / 7) * (cx + r * cos(t)) - sin(9 * pi / 7) * r * sin(t);\ny = sin(9 * pi / 7) * (cx + r * cos(t)) + cos(9 * pi / 7) * r * sin(t); label = 10; };\nborder G11(t = 9 * pi / 7,5 * pi / 7){x = cos(10 * pi / 7) * (cx + r * cos(t)) - sin(10 * pi / 7) * r * sin(t);\ny = sin(10 * pi / 7) * (cx + r * cos(t)) + cos(10 * pi / 7) * r * sin(t); label = 11; };\nborder G12(t = 9 * pi / 7,5 * pi / 7){x = cos(11 * pi / 7) * (cx + r * cos(t)) - sin(11 * pi / 7) * r * sin(t);\ny = sin(11 * pi / 7) * (cx + r * cos(t)) + cos(11 * pi / 7) * r * sin(t); label = 12; };\nborder G13(t = 9 * pi / 7,5 * pi / 7){x = cos(12 * pi / 7) * (cx + r * cos(t)) - sin(12 * pi / 7) * r * sin(t);\ny = sin(12 * pi / 7) * (cx + r * cos(t)) + cos(12 * pi / 7) * r * sin(t); label = 13; };\nborder G14(t = 9 * pi / 7,5 * pi / 7){x = cos(13 * pi / 7) * (cx + r * cos(t)) - sin(13 * pi / 7) * r * sin(t);\ny = sin(13 * pi / 7) * (cx + r * cos(t)) + cos(13 * pi / 7) * r * sin(t); label = 14; };\n\n' + f'int n = {precision}' + ';\n\n//plot(T1(precision) + T2(precision) + T3(precision) + T4(precision));\n\nmesh Th = buildmesh(G1(n) + G2(n) + G3(n) + G4(n) + G5(n) + G6(n) + G7(n) + G8(n) + G9(n) + G10(n) + G11(n) + G12(n) + G13(n) + G14(n), fixeborder = true);\n\nfespace Vh(Th, P2, periodic = [[1, y], [6, (cos(2 * pi/7)) * (y - y06) - (sin( - 2 * pi/7)) * (x - x06)],\n                                [3, (cos(2 * pi/7)) * (y - y03) - (sin(2 * pi/7)) * (x - x03)], [8, y], \n                                [5, (cos(4 * pi/7)) * (y - y05) - (sin(4 * pi/7)) * (x - x05)], \n                                [10, (cos(2 * pi/7)) * (y - y10) - (sin(2 * pi/7)) * (x - x10)], \n                                [7, (cos(6 * pi/7)) * (y - y07) - (sin(6 * pi/7)) * (x - x07)], \n                                [12, (cos(4 * pi/7)) * (y - y12) - (sin(4 * pi/7)) * (x - x12)], \n                                [9, (cos(6 * pi/7)) * (y - y09) - (sin( - 6 * pi/7)) * (x - x09)], \n                                [14, (cos(6 * pi/7)) * (y - y14) - (sin(6 * pi/7)) * (x - x14)], \n                                [11, (cos(4 * pi/7)) * (y - y11) - (sin( - 4 * pi/7)) * (x - x11)], \n                                [2, (cos(6 * pi/7)) * (y - y02) - (sin( - 6 * pi/7)) * (x - x02)], \n                                [13, (cos(2 * pi/7)) * (y - y13) - (sin( - 2 * pi/7)) * (x - x13)], \n                                [4, (cos(4 * pi/7)) * (y - y04) - (sin( - 4 * pi/7)) * (x - x04)]]);\n\nVh u1, u2;\nreal sigma = 0.00001;\n\nvarf op (u1, u2) = int2d(Th)(dx(u1) * dx(u2) + dy(u1) * dy(u2) - sigma * u1 * u2 * (' + f'{genererFacteurConforme(divX, divY, params)}));\nvarf bm ([u1], [u2]) = int2d(Th)(u1 * u2 * ({genererFacteurConforme(divX, divY, params)}' + '));\n\nmatrix OP = op(Vh, Vh, solver = Crout, factorize = 1);\nmatrix B = bm(Vh, Vh, solver = CG, eps = 1e-20);\n\nreal[int] ev(30);\nVh[int] eV(30);\n\nint k = EigenValue(OP, B, sym = true, sigma = sigma, value = ev, vector = eV, tol=1e-10, maxit = 0, ncv = 0);\n\nreal aire = int2d(Th)(' + f'{genererFacteurConforme(divX, divY, params)}' + ');\nreal vp1 = 0;\nint nbr = 0;\n\nfor (int i = 0; i < k; i++)\n{\n    cout << "lambda[" << i << "] = " << ev[i] << " | normalisation : " << ev[i] * aire << endl;\n\n    if (abs(ev[i] * aire) > 10e-2)\n    {\n        vp1 = ev[i] * aire;\n        nbr = i;\n        break;\n    }\n}\n\ncout << "Retour : " << vp1 << "|" << endl;\n\ncout << "GRADIENT [" << endl;\n\nint pos = 0;\n' + f'real denom = int2d(Th)(({genererFacteurConforme(divX, divY, params)}) * eV[nbr]^2);' + '\n\nfor (int i = 0; i < divX + 1; i++)\n{\n    for (int j = 0; j < divY + 1; j++)\n    {\n        cout << -' + f'ev[nbr] * int2d(Th)((aire * eV[nbr]^2 / (denom) - 1) * exp(- ( (x - 2 * i / divX + 1)^2 + (y - 2 * j / divY + 1)^2  ) / 0.001))' + ' << endl;\n        pos++;\n    }\n}\n\ncout << "]" << endl;\ncout << "Aire : " << aire << endl;'

    fichier = open(f"{chemin}{nom}.edp", "w")
    fichier.write(programme)
    fichier.close()

    texte = ""

    if cheminFF != "":
        texte = os.popen(f"{cheminFF} {chemin}{nom}.edp").read()
    else:
        texte = os.popen(f"{chemin}{nom}.edp").read()

    val = float(texte.split("Retour : ")[-1].split("|")[0])
    aire = float(texte.split("Aire : ")[-1].split("\n")[0])

    print(f"vpn = {val}, vp = {val / aire}, aire = {aire}")

    return (val, [float(texte.split("GRADIENT [")[-1].split("]")[0].split("\n")[i]) for i in range((divX + 1) * (divY + 1) + 1) if texte.split("GRADIENT [")[-1].split("]")[0].split("\n")[i] != ''])

def vp(params):
    programme = f'int divX = {divX};\nint divY = {divY};' + '\n\nreal angle1 = 0.4487989505128276;\n\nreal cx = 1.031969722;\nreal r = 0.254875473;\n\nreal s = 0.7770942488589633;\n\nreal x06 = s * cos(5 * pi / 7);\nreal y06 = s * sin(5 * pi / 7);\nreal x03 = s * cos(2 * pi / 7);\nreal y03 = s * sin(2 * pi / 7);\nreal x05 = s * cos(4 * pi / 7);\nreal y05 = s * sin(4 * pi / 7);\nreal x10 = s * cos(9 * pi / 7);\nreal y10 = s * sin(9 * pi / 7);\nreal x07 = s * cos(6 * pi / 7);\nreal y07 = s * sin(6 * pi / 7);\nreal x12 = s * cos(11 * pi / 7);\nreal y12 = s * sin(11 * pi / 7);\nreal x09 = s * cos(8 * pi / 7);\nreal y09 = s * sin(8 * pi / 7);\nreal x14 = s * cos(13 * pi / 7);\nreal y14 = s * sin(13 * pi / 7);\nreal x11 = s * cos(10 * pi / 7);\nreal y11 = s * sin(10 * pi / 7);\nreal x02 = s * cos(pi / 7);\nreal y02 = s * sin(pi / 7);\nreal x13 = s * cos(12 * pi / 7);\nreal y13 = s * sin(12 * pi / 7);\nreal x04 = s * cos(3 * pi / 7);\nreal y04 = s * sin(3 * pi / 7);\n\nborder G1 (t = 9 * pi / 7,5 * pi / 7){ x = cx + r * cos(t); y = r * sin(t); label = 1; };\nborder G2 (t = 9 * pi / 7,5 * pi / 7){ x = cos(pi / 7) * (cx + r * cos(t)) - sin(pi / 7) * (r * sin(t));\ny = sin(pi / 7) * (cx + r * cos(t)) + cos(pi / 7) * (r * sin(t)); label = 2; };\nborder G3 (t = 9 * pi / 7,5 * pi / 7){x = cos(2 * pi / 7) * (cx + r * cos(t)) - sin(2 * pi / 7) * r * sin(t);\ny = sin(2 * pi / 7) * (cx + r * cos(t)) + cos(2 * pi / 7) * r * sin(t); label = 3; };\nborder G4 (t = 9 * pi / 7,5 * pi / 7){x = cos(3 * pi / 7) * (cx + r * cos(t)) - sin(3 * pi / 7) * r * sin(t);\ny = sin(3 * pi / 7) * (cx + r * cos(t)) + cos(3 * pi / 7) * r * sin(t); label = 4; };\nborder G5 (t = 9 * pi / 7,5 * pi / 7){ x = cos(4 * pi / 7) * (cx + r * cos(t)) - sin(4 * pi / 7) * r * sin(t);\ny = sin(4 * pi / 7) * (cx + r * cos(t)) + cos(4 * pi / 7) * r * sin(t); label = 5; };\nborder G6(t = 9 * pi / 7,5 * pi / 7){x = cos(5 * pi / 7) * (cx + r * cos(t)) - sin(5 * pi / 7) * r * sin(t);\ny = sin(5 * pi / 7) * (cx + r * cos(t)) + cos(5 * pi / 7) * r * sin(t); label = 6; };\nborder G7(t = 9 * pi / 7,5 * pi / 7){x = cos(6 * pi / 7) * (cx + r * cos(t)) - sin(6 * pi / 7) * r * sin(t);\ny = sin(6 * pi / 7) * (cx + r * cos(t)) + cos(6 * pi / 7) * r * sin(t); label = 7; };\nborder G8(t = 9 * pi / 7,5 * pi / 7){x = cos(7 * pi / 7) * (cx + r * cos(t)) - sin(7 * pi / 7) * r * sin(t);\ny = sin(7 * pi / 7) * (cx + r * cos(t)) + cos(7 * pi / 7) * r * sin(t); label = 8; };\nborder G9(t = 9 * pi / 7,5 * pi / 7){x = cos(8 * pi / 7) * (cx + r * cos(t)) - sin(8 * pi / 7) * r * sin(t);\ny = sin(8 * pi / 7) * (cx + r * cos(t)) + cos(8 * pi / 7) * r * sin(t); label = 9; };\nborder G10(t = 9 * pi / 7,5 * pi / 7){x = cos(9 * pi / 7) * (cx + r * cos(t)) - sin(9 * pi / 7) * r * sin(t);\ny = sin(9 * pi / 7) * (cx + r * cos(t)) + cos(9 * pi / 7) * r * sin(t); label = 10; };\nborder G11(t = 9 * pi / 7,5 * pi / 7){x = cos(10 * pi / 7) * (cx + r * cos(t)) - sin(10 * pi / 7) * r * sin(t);\ny = sin(10 * pi / 7) * (cx + r * cos(t)) + cos(10 * pi / 7) * r * sin(t); label = 11; };\nborder G12(t = 9 * pi / 7,5 * pi / 7){x = cos(11 * pi / 7) * (cx + r * cos(t)) - sin(11 * pi / 7) * r * sin(t);\ny = sin(11 * pi / 7) * (cx + r * cos(t)) + cos(11 * pi / 7) * r * sin(t); label = 12; };\nborder G13(t = 9 * pi / 7,5 * pi / 7){x = cos(12 * pi / 7) * (cx + r * cos(t)) - sin(12 * pi / 7) * r * sin(t);\ny = sin(12 * pi / 7) * (cx + r * cos(t)) + cos(12 * pi / 7) * r * sin(t); label = 13; };\nborder G14(t = 9 * pi / 7,5 * pi / 7){x = cos(13 * pi / 7) * (cx + r * cos(t)) - sin(13 * pi / 7) * r * sin(t);\ny = sin(13 * pi / 7) * (cx + r * cos(t)) + cos(13 * pi / 7) * r * sin(t); label = 14; };\n\n' + f'int n = {precision}' + ';\n\n//plot(T1(precision) + T2(precision) + T3(precision) + T4(precision));\n\nmesh Th = buildmesh(G1(n) + G2(n) + G3(n) + G4(n) + G5(n) + G6(n) + G7(n) + G8(n) + G9(n) + G10(n) + G11(n) + G12(n) + G13(n) + G14(n), fixeborder = true);\n\nfespace Vh(Th, P2, periodic = [[1, y], [6, (cos(2 * pi/7)) * (y - y06) - (sin( - 2 * pi/7)) * (x - x06)],\n                                [3, (cos(2 * pi/7)) * (y - y03) - (sin(2 * pi/7)) * (x - x03)], [8, y], \n                                [5, (cos(4 * pi/7)) * (y - y05) - (sin(4 * pi/7)) * (x - x05)], \n                                [10, (cos(2 * pi/7)) * (y - y10) - (sin(2 * pi/7)) * (x - x10)], \n                                [7, (cos(6 * pi/7)) * (y - y07) - (sin(6 * pi/7)) * (x - x07)], \n                                [12, (cos(4 * pi/7)) * (y - y12) - (sin(4 * pi/7)) * (x - x12)], \n                                [9, (cos(6 * pi/7)) * (y - y09) - (sin( - 6 * pi/7)) * (x - x09)], \n                                [14, (cos(6 * pi/7)) * (y - y14) - (sin(6 * pi/7)) * (x - x14)], \n                                [11, (cos(4 * pi/7)) * (y - y11) - (sin( - 4 * pi/7)) * (x - x11)], \n                                [2, (cos(6 * pi/7)) * (y - y02) - (sin( - 6 * pi/7)) * (x - x02)], \n                                [13, (cos(2 * pi/7)) * (y - y13) - (sin( - 2 * pi/7)) * (x - x13)], \n                                [4, (cos(4 * pi/7)) * (y - y04) - (sin( - 4 * pi/7)) * (x - x04)]]);\n\nVh u1, u2;\nreal sigma = 0.00001;\n\nvarf op (u1, u2) = int2d(Th)(dx(u1) * dx(u2) + dy(u1) * dy(u2) - sigma * u1 * u2 * (' + f'{genererFacteurConforme(divX, divY, params)}));\nvarf bm ([u1], [u2]) = int2d(Th)(u1 * u2 * ({genererFacteurConforme(divX, divY, params)}' + '));\n\nmatrix OP = op(Vh, Vh, solver = Crout, factorize = 1);\nmatrix B = bm(Vh, Vh, solver = CG, eps = 1e-20);\n\nreal[int] ev(30);\nVh[int] eV(30);\n\nint k = EigenValue(OP, B, sym = true, sigma = sigma, value = ev, vector = eV, tol=1e-10, maxit = 0, ncv = 0);\n\nreal aire = int2d(Th)(' + f'{genererFacteurConforme(divX, divY, params)}' + ');\nreal vp1 = 0;\nint nbr = 0;\n\nfor (int i = 0; i < k; i++)\n{\n    cout << "lambda[" << i << "] = " << ev[i] << " | normalisation : " << ev[i] * aire << endl;\n\n    if (abs(ev[i] * aire) > 10e-2)\n    {\n        vp1 = ev[i] * aire;\n        nbr = i;\n        break;\n    }\n}\n\ncout << "Retour : " << vp1 << "|" << endl;\ncout << "Aire : " << aire << endl;'

    fichier = open(f"{chemin}{nom}.edp", "w")
    fichier.write(programme)
    fichier.close()

    texte = ""

    if cheminFF != "":
        texte = os.popen(f"{cheminFF} {chemin}{nom}.edp").read()
    else:
        texte = os.popen(f"{chemin}{nom}.edp").read()

    return float(texte.split("Retour : ")[-1].split("|")[0])


x = np.linspace(-0.6, 0.6, 500)
y = np.linspace(-0.6, 0.6, 500)

x, y = np.meshgrid(x, y)

z = [facteurConforme([x[i], y[i]]) for i in range(len(x))]

z = np.array(z)

fig = plt.figure(figsize = plt.figaspect(0.5))
ax = fig.add_subplot(1, 2, 1, projection='3d')

# Plot the surface.
surf = ax.plot_surface(x, y, z, cmap = cm.coolwarm, linewidth = 0, antialiased = True)

plt.show()

cible = (16 * (4 - np.sqrt(7)) * np.pi)

def rechercheLineaire(x, gradient):
    valPermises = [-1000, -500, -200, -100, -50, -10, -5, -2, -1, -0.1, -0.01, -0.001, 0.001, 0.01, 0.1, 1, 2, 5, 10, 50, 100, 200, 500, 1000]

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
    
    return valPermises[iMax]

base = (0, [0] * len(coeff))

while abs(cible - base[0]) > 0.01:
    base = valeurPropre(coeff)

    eta = rechercheLineaire(coeff, base[1])

    for i in range(len(coeff)):
        if base[1][i] != 0 and approx(abs(base[1][i] * 100), 1) != 0:
            coeff[i] = coeff[i] - eta * base[1][i] / approx(abs(base[1][i] * 100), 1)

    fichier = open("path", "a")
    fichier.write(str([base[0], coeff]))
    fichier.write("\n")
    fichier.close()