import os

# Nom du fichier .edp
nom = "trianglesSpheriques"

# Chemin du fichier .edp
chemin = "path"

# Chemin de l'exécutable modifié, si besoin, de FreeFem++
cheminFF = "path"

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

texte = f"\nint nbrVP = 10;\nint n = 10;\n\n" + "real a = acos(cos(4 * pi / 7) / sin(pi / 3));\nreal b = acos(cos(pi / 3) / sin(4 * pi / 7));\nreal c = acos(1 / tan(4 * pi / 7) / tan(pi / 3));\n\nfunc complex projection(real[int] x)\n{\n    if (x[2] == 1)\n    {\n        return 0;\n    }\n\n    return (x[0] + x[1] * 1i) / (1 - x[2]);\n}\n\nfunc real norme(real[int] v)\n{\n    real somme = 0;\n\n    for (int i = 0; i < v.n; i++)\n    {\n        somme += v[i]^2;\n    }\n\n    return sqrt(somme);\n}\n\nfunc real[int] produitVectoriel(real[int] v1, real[int] v2)\n{\n    real[int] w(3);\n\n    w[0] = v1[1] * v2[2] - v1[2] * v2[1];\n    w[1] = v1[2] * v2[0] - v1[0] * v2[2];\n    w[2] = v1[0] * v2[1] - v1[1] * v2[0];\n\n    return w;\n} \n\nfunc real produitScalaire(real[int] v1, real[int] v2)\n{\n    if (v1.n != v2.n)\n    {\n        return 0;\n    }\n\n    real somme = 0;\n\n    for (int i = 0; i < v1.n; i++)\n    {\n        somme += v1[i] * v2[i];\n    }\n\n    return somme;\n}\n\nfunc complex projectionComplete(real angle, real[int] k, real[int] A)\n{\n    real[int] xphi(3);\n    real[int] produit(3);\n\n    xphi += cos(angle) * A;\n    xphi += (1 - cos(angle)) * produitScalaire(k, A) * k;\n    xphi += sin(angle) * produitVectoriel(k, A);\n\n    return projection(xphi);\n\n    return 1;\n}\n\nreal[int] A(3);\nreal[int] B(3);\nreal[int] k(3);\n\nA = [cos(b - pi / 2), 0, sin(b - pi / 2)];\nB = [0, cos(a - pi / 2), sin(a - pi / 2)];\n\nk = produitVectoriel(A, B) / norme(produitVectoriel(A, B));\n\nreal theta = acos(produitScalaire(A, B) / (norme(A) * norme(B)));\n\n"

for i in range(336):
    texte += f"border T{i + 1}C1 (t = 0, imag(projection(B)))" + " {" + f"x = {i} * 10; y = t; label = " + f"{3 * i + 1}" + ";};\n" + f"border T{i + 1}C2 (t = 0, real(projection(A)))" + " {" + f"x = {i} * 10 + t; y = 0; label = {3 * i + 2};" + "};\n" + f"border T{i + 1}C3 (t = 0, theta)" + " {" + f"x = real(projectionComplete(t, k, A)) + {i} * 10; " + f"y = imag(projectionComplete(t, k, A)); label = {3 * (i + 1)};" + "};\n"

texte += "\n\nmesh Th = buildmesh("

for i in range(336):
    texte += f"T{i + 1}C1(-n) + T{i + 1}C2(n) + T{i + 1}C3(n)"

    if i != 335:
        texte += " + "

texte += ");"

fichier = open(chemin + nom + ".edp", "w")
fichier.write(texte)
fichier.close()

conditionsPeriodiques = "["

for i in range(len(associations)):
    if associations[i][1] == 0:
        conditionsPeriodiques += f"[{3 * associations[i][0][0] + 1}, y], [{3 * associations[i][0][1] + 1}, y]"

    elif associations[i][1] == 1:
        conditionsPeriodiques += f"[{3 * associations[i][0][0] + 2}, x], [{3 * associations[i][0][1] + 2}, x - {associations[i][0][1] - associations[i][0][0]} * 10]"

    else:
        conditionsPeriodiques += f"[{3 * (associations[i][0][0] + 1)}, y], [{3 * (associations[i][0][1] + 1)}, y]"

    if i != len(associations) - 1:
        conditionsPeriodiques += ",\n"

conditionsPeriodiques += "]"

texte = '\nfespace Vh (Th, P2, periodic = ' + conditionsPeriodiques + ');\n\n' + f'Vh u1, u2;\nreal sigma = 0.00001;\n\nvarf op (u1, u2) = int2d(Th)(dx(u1) * dx(u2) + dy(u1) * dy(u2) - sigma * u1 * u2 * (1 / (1 + (x - floor(x))^2 + y^2)^2));\nvarf bm ([u1], [u2]) = int2d(Th)(u1 * u2 * (1 / (1 + (x - floor(x))^2 + y^2)^2));' + '\n\nmatrix OP = op(Vh, Vh, solver = Crout, factorize = 1);\nmatrix BM = bm(Vh, Vh, solver = CG, eps = 1e-20);\n\nreal[int] ev(nbrVP);\nVh[int] eV(nbrVP);\n\nint K = EigenValue(OP, BM, sym = true, sigma = sigma, value = ev, vector = eV, tol=1e-10, maxit = 0, ncv = 0);\nreal val = 0;\n' + f'real aire = int2d(Th)((1 / (1 + (x - floor(x))^2 + y^2)^2));' + '\nreal vp1 = 0;\nint nbr = 0;\n\nfor (int i = 0; i < K; i++)\n{\n    cout << "lambda[" << i << "] = " << ev[i] << " | normalisation : " << ev[i] * aire << endl;\n\n    if (abs(ev[i] * aire) > 10e-2)\n    {\n        vp1 = ev[i] * aire;\n        break;\n        nbr = i;\n    }\n}\n\ncout << "Retour : " << vp1 << "|" << endl;\ncout << "Aire : " << aire << endl;'

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
