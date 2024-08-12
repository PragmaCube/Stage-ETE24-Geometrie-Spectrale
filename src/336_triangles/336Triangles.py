import os
import numpy as np
import random
import math

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
nom = "336Triangles"

# Chemin du fichier .edp
chemin = "path"

# Chemin de l'exécutable modifié, si besoin, de FreeFem++
cheminFF = "path"

angle = 1.102

longueur = 0.5

divX = 10
divY = 10

fichier = open(chemin + "Isom(K).txt", "r")

associations = []

for elem in fichier:
    associations.append(eval(elem.split(",\n")[0]))

fichier.close()

for elem in associations:
    elem[0].sort()

fini = False

while not fini:
    recommencer = False

    for i in range(len(associations)):
        for j in range(len(associations)):
            if i != j and associations[i] == associations[j]:
                del associations[j]
                recommencer = True
                break

        if recommencer:
            break

    else:
        fini = True

def genererFacteurConforme(divX, divY, coeff):
    pos = 0

    facteurConforme = ""

    for i in range(divX + 1):
        for j in range(divY + 1):
            facteurConforme += f"{coeff[pos]} * exp(- ( (x - {i / divX} - floor(x))^2 + (y - {j / divY})^2  ) / 0.1)"
            
            if pos != (divX + 1) * (divY + 1) - 1:
                facteurConforme += " + "

            pos += 1

    return facteurConforme

def facteurConforme(x, coeff):
    pos = 0

    somme = 0

    for i in range(divX + 1):
        for j in range(divY + 1):
            somme += coeff[pos] * np.exp(- ((x[0] - i / divX + 1) ** 2 + (x[1] - j / divY + 1) ** 2) / 0.1)

            pos += 1

    return somme

coeff = [random.random() for i in range((divX + 1) * (divY + 1))]

def valeurPropre(params):
    texte = f"int divX = {divX};\nint divY = {divY};\nint nbrVP = 10;\nint n = 10;\n\n"

    for i in range(336):
        texte += f"border T{i + 1}C1 (t = {2 * i}, {2 * i + 1})" + " {x = t; y = 0; label = " + f"{3 * i + 1}" + ";};\n" + f"border T{i + 1}C2 (t = 0, tan(pi / 2 - {angle}))" + " {" + f"x = {2 * i}; y = t; label = {3 * i + 2};" + "};\n" + f"border T{i + 1}C3 (t = {2 * i}, {2 * i + 1})" + " {x = t; " + f"y = - tan(pi / 2 - {angle}) * (t - {2 * i}) + tan(pi / 2 - {angle}); label = {3 * (i + 1)};" + "};\n"

    texte += "\n\nmesh Th = buildmesh("

    for i in range(336):
        texte += f"T{i + 1}C1(n) + T{i + 1}C2(-n) + T{i + 1}C3(-n)"

        if i != 335:
            texte += " + "

    texte += ");"

    fichier = open(chemin + nom + ".edp", "w")
    fichier.write(texte)
    fichier.close()

    conditionsPeriodiques = "["

    for i in range(len(associations)):
        if associations[i][1] == 0:
            conditionsPeriodiques += f"[{3 * associations[i][0][0] + 1}, x], [{3 * associations[i][0][1] + 1}, x - 2 * {associations[i][0][1] - associations[i][0][0]}]"

        elif associations[i][1] == 1:
            conditionsPeriodiques += f"[{3 * associations[i][0][0] + 2}, y], [{3 * associations[i][0][1] + 2}, y]"

        else:
            conditionsPeriodiques += f"[{3 * (associations[i][0][0] + 1)}, y], [{3 * (associations[i][0][1] + 1)}, y]"

        if i != len(associations) - 1:
            conditionsPeriodiques += ",\n"

    conditionsPeriodiques += "]"
    
    texte = '\nfespace Vh (Th, P2, periodic = ' + conditionsPeriodiques + ');\n\n' + f'Vh u1, u2;\nreal sigma = 0.00001;\n\nvarf op (u1, u2) = int2d(Th)(dx(u1) * dx(u2) + dy(u1) * dy(u2) - sigma * u1 * u2 * ({genererFacteurConforme(divX, divY, params)}));\nvarf bm ([u1], [u2]) = int2d(Th)(u1 * u2 * ({genererFacteurConforme(divX, divY, params)}));' + '\n\nmatrix OP = op(Vh, Vh, solver = Crout, factorize = 1);\nmatrix B = bm(Vh, Vh, solver = CG, eps = 1e-20);\n\nreal[int] ev(nbrVP);\nVh[int] eV(nbrVP);\n\nint k = EigenValue(OP, B, sym = true, sigma = sigma, value = ev, vector = eV, tol=1e-10, maxit = 0, ncv = 0);\nreal val = 0;\n' + f'real aire = int2d(Th)(({genererFacteurConforme(divX, divY, params)}));' + '\nreal vp1 = 0;\nint nbr = 0;\n\nfor (int i = 0; i < k; i++)\n{\n    cout << "lambda[" << i << "] = " << ev[i] << " | normalisation : " << ev[i] * aire << endl;\n\n    if (abs(ev[i] * aire) > 10e-2)\n    {\n        vp1 = ev[i] * aire;\n        nbr = i;\n        break;\n    }\n}\n\ncout << "Retour : " << vp1 << "|" << endl;\n\ncout << "GRADIENT [" << endl;\n\nint pos = 0;\n' + f'real denom = int2d(Th)(({genererFacteurConforme(divX, divY, params)}) * eV[nbr]^2);' + '\n\nfor (int i = 0; i < divX + 1; i++)\n{\n    for (int j = 0; j < divY + 1; j++)\n    {\n        cout << -' + f'ev[nbr] * int2d(Th)((aire * eV[nbr]^2 / (denom) - 1) * exp(- ((x - i / divX + 1 - floor(x))^2 + (y - j / divY + 1)^2 ) / 0.1))' + ' << endl;\n        pos++;\n    }\n}\n\ncout << "]" << endl;\ncout << "Aire : " << aire << endl;'

    fichier = open(chemin + nom + ".edp", "a")
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

    return 0

def vp(params):
    texte = f"int divX = {divX};\nint divY = {divY};\nint nbrVP = 10;\nint n = 10;\n\n"

    for i in range(336):
        texte += f"border T{i + 1}C1 (t = {2 * i}, {2 * i + 1})" + " {x = t; y = 0; label = " + f"{3 * i + 1}" + ";};\n" + f"border T{i + 1}C2 (t = 0, tan(pi / 2 - {angle}))" + " {" + f"x = {2 * i}; y = t; label = {3 * i + 2};" + "};\n" + f"border T{i + 1}C3 (t = {2 * i}, {2 * i + 1})" + " {x = t; " + f"y = - tan(pi / 2 - {angle}) * (t - {2 * i}) + tan(pi / 2 - {angle}); label = {3 * (i + 1)};" + "};\n"

    texte += "\n\nmesh Th = buildmesh("

    for i in range(336):
        texte += f"T{i + 1}C1(n) + T{i + 1}C2(-n) + T{i + 1}C3(-n)"

        if i != 335:
            texte += " + "

    texte += ");"

    fichier = open(chemin + nom + ".edp", "w")
    fichier.write(texte)
    fichier.close()

    conditionsPeriodiques = "["

    for i in range(len(associations)):
        if associations[i][1] == 0:
            conditionsPeriodiques += f"[{3 * associations[i][0][0] + 1}, x], [{3 * associations[i][0][1] + 1}, x - 2 * {associations[i][0][1] - associations[i][0][0]}]"

        elif associations[i][1] == 1:
            conditionsPeriodiques += f"[{3 * associations[i][0][0] + 2}, y], [{3 * associations[i][0][1] + 2}, y]"

        else:
            conditionsPeriodiques += f"[{3 * (associations[i][0][0] + 1)}, y], [{3 * (associations[i][0][1] + 1)}, y]"

        if i != len(associations) - 1:
            conditionsPeriodiques += ",\n"

    conditionsPeriodiques += "]"
    
    texte = '\nfespace Vh (Th, P2, periodic = ' + conditionsPeriodiques + ');\n\n' + f'Vh u1, u2;\nreal sigma = 0.00001;\n\nvarf op (u1, u2) = int2d(Th)(dx(u1) * dx(u2) + dy(u1) * dy(u2) - sigma * u1 * u2 * ({genererFacteurConforme(divX, divY, params)}));\nvarf bm ([u1], [u2]) = int2d(Th)(u1 * u2 * ({genererFacteurConforme(divX, divY, params)}));' + '\n\nmatrix OP = op(Vh, Vh, solver = Crout, factorize = 1);\nmatrix B = bm(Vh, Vh, solver = CG, eps = 1e-20);\n\nreal[int] ev(nbrVP);\nVh[int] eV(nbrVP);\n\nint k = EigenValue(OP, B, sym = true, sigma = sigma, value = ev, vector = eV, tol=1e-10, maxit = 0, ncv = 0);\nreal val = 0;\n' + f'real aire = int2d(Th)(({genererFacteurConforme(divX, divY, params)}));' + '\nreal vp1 = 0;\nint nbr = 0;\n\nfor (int i = 0; i < k; i++)\n{\n    cout << "lambda[" << i << "] = " << ev[i] << " | normalisation : " << ev[i] * aire << endl;\n\n    if (abs(ev[i] * aire) > 10e-2)\n    {\n        vp1 = ev[i] * aire;\n        nbr = i;\n        break;\n    }\n}\n\ncout << "Retour : " << vp1 << "|" << endl;\ncout << "Aire : " << aire << endl;'

    fichier = open(chemin + nom + ".edp", "a")
    fichier.write(texte)
    fichier.close()

    if cheminFF != "":
        texte = os.popen(f"{cheminFF} {chemin}{nom}.edp").read()
    else:
        texte = os.popen(f"{chemin}{nom}.edp").read()

    if len(texte.split("Exec error ")) == 1:

        return float(texte.split("Retour : ")[-1].split("|")[0])

    return 0


cible = (16 * (4 - np.sqrt(7)) * np.pi)

temp = []
def rechercheLineaire(x, gradient):
    valPermises = [-1000, -500, -200, -100, -50, -10, -5, -2, -1, -0.1, -0.01, -0.001, 0.001, 0.01, 0.1, 1, 2, 5, 10, 50, 100, 200, 500, 1000]

    resultats = []

    for i in range(len(valPermises)):
        temp = list(x)

        for j in range(len(c2)):
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