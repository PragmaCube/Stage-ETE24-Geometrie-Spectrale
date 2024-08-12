import os

from scipy.optimize import minimize

###########################################################
# Nombre de valeurs propres à calculer
nbrVP = 10

# Nombre de divisions par côté
precision = 25

signeBords = ["-", "-", "-", "-", "", "-", "-", "", "", "-", "", "-", "-", ""]

# Nom du fichier .edp
nom = "14gone"

# Chemin du fichier .edp
chemin = "path"

# Chemin de l'exécutable modifié, si besoin, de FreeFem++
cheminFF = "path"

# Chemin du fichier de sauvegarde des données
sauvegarde = "path"

def presqueNul(nbr, tol):
    if abs(nbr) < tol:
        return 0
    
    return nbr

def valeurPropre(params):
    texte = ""

    while True: # v[0] = 1 - coeffb * exp(1i * pi / 7) + {params[0]};\nv[1] = exp(2i * pi / 7) - coeffb * exp(3i * pi / 7) + {params[1]};\nv[2] = exp(4i * pi / 7) - coeffb * exp(5i * pi / 7) + {params[2]};\nv[3] = exp(6i * pi / 7) - coeffb * exp(7i * pi / 7) + {params[3]};\nv[4] = exp(8i * pi / 7) - coeffb * exp(9i * pi / 7) + {params[4]};\nv[5] = exp(10i * pi / 7) - coeffb * exp(11i * pi / 7) + {params[5]};\nv[6] = exp(12i * pi / 7) - coeffb * exp(13i * pi / 7) + {params[6]};
        programme = f'int n = {precision};\nint nbrVP = {nbrVP};\n\nreal coeffb = sin(2 * pi / 7) / sin(4 * pi / 7);\n\ncomplex[int] v(7);\nv[0] = 1 - coeffb * exp(1i * pi / 7) + {params[7]}i + {params[0]};\nv[1] = exp(2i * pi / 7) - coeffb * exp(3i * pi / 7) + {params[8]}i + {params[1]};\nv[2] = exp(4i * pi / 7) - coeffb * exp(5i * pi / 7) + {params[9]}i + {params[2]};\nv[3] = exp(6i * pi / 7) - coeffb * exp(7i * pi / 7) + {params[10]}i + {params[3]};\nv[4] = exp(8i * pi / 7) - coeffb * exp(9i * pi / 7) + {params[11]}i + {params[4]};\nv[5] = exp(10i * pi / 7) - coeffb * exp(11i * pi / 7) + {params[12]}i + {params[5]};\nv[6] = exp(12i * pi / 7) - coeffb * exp(13i * pi / 7) + {params[13]}i + {params[6]};' + '\n\ncomplex[int] points(14);\npoints[0] = 0;\npoints[1] = v[0];\npoints[2] = v[0] - v[5];\npoints[3] = v[0] - v[5] + v[1];\npoints[4] = v[0] - v[5] + v[1] - v[6];\npoints[5] = v[0] - v[5] + v[1] - v[6] + v[2];\npoints[6] = - v[5] + v[1] - v[6] + v[2];\npoints[7] = -v[5] + v[1] - v[6] + v[2] + v[3];\npoints[8] = -v[5] - v[6] + v[2] + v[3];\npoints[9] = -v[5] - v[6] + v[2] + v[3] + v[4];\npoints[10] = -v[5] - v[6] + v[3] + v[4];\npoints[11] = -v[6] + v[3] + v[4];\npoints[12] = -v[6] + v[4];\npoints[13] = v[4];\n\nreal[int] pentes(14);\nint[int] obs(14);\nint[int] cond(7);\n\nfor (int i = 0; i < 14; i++)\n{\n    pentes[i] = 1;\n    obs[i] = 1;\n\n    if (real(points[(i + 1) % 14] - points[i]) != 0.0)\n    {\n        pentes[i] = imag(points[(i + 1) % 14] - points[i]) / real(points[(i + 1) % 14] - points[i]);\n        obs[i] = 0;\n    }\n}\n\nfor (int i = 0; i < 7; i++)\n{\n    cond[i] = 1;\n\n    if (real(points[(2 * i + 5) % 14] - points[(2 * i + 1) % 14] != 0.0))\n    {\n        cond[i] = 0;\n    }\n}\n\nborder C1 (t = 0, real(points[1]) + (obs[0]) * imag(points[1])) {x = !obs[0] * t; y = pentes[0] * t; label = 1;};\nborder C2 (t = !obs[1] * real(points[1]) + obs[1] * imag(points[1]), !obs[1] * real(points[2]) + obs[1] * imag(points[2])) {x = !obs[1] * t + obs[1] * real(points[1]); y = pentes[1] * (t - !obs[1] * real(points[1])) + !obs[1] * imag(points[1]); label = 2;};\nborder C3 (t = !obs[2] * real(points[2]) + obs[2] * imag(points[2]), !obs[2] * real(points[3]) + obs[2] * imag(points[3])) {x = !obs[2] * t + obs[2] * real(points[2]); y = pentes[2] * (t - !obs[2] * real(points[2])) + !obs[2] * imag(points[2]); label = 3;};\nborder C4 (t = !obs[3] * real(points[3]) + obs[3] * imag(points[3]), !obs[3] * real(points[4]) + obs[3] * imag(points[4])) {x = !obs[3] * t + obs[3] * real(points[3]); y = pentes[3] * (t - !obs[3] * real(points[3])) + !obs[3] * imag(points[3]); label = 4;};\nborder C5 (t = !obs[4] * real(points[4]) + obs[4] * imag(points[4]), !obs[4] * real(points[5]) + obs[4] * imag(points[5])) {x = !obs[4] * t + obs[4] * real(points[4]); y = pentes[4] * (t - !obs[4] * real(points[4])) + !obs[4] * imag(points[4]); label = 5;};\nborder C6 (t = !obs[5] * real(points[5]) + obs[5] * imag(points[5]), !obs[5] * real(points[6]) + obs[5] * imag(points[6])) {x = !obs[5] * t + obs[5] * real(points[5]); y = pentes[5] * (t - !obs[5] * real(points[5])) + !obs[5] * imag(points[5]); label = 6;};\nborder C7 (t = !obs[6] * real(points[6]) + obs[6] * imag(points[6]), !obs[6] * real(points[7]) + obs[6] * imag(points[7])) {x = !obs[6] * t + obs[6] * real(points[6]); y = pentes[6] * (t - !obs[6] * real(points[6])) + !obs[6] * imag(points[6]); label = 7;};\nborder C8 (t = !obs[7] * real(points[7]) + obs[7] * imag(points[7]), !obs[7] * real(points[8]) + obs[7] * imag(points[8])) {x = !obs[7] * t + obs[7] * real(points[7]); y = pentes[7] * (t - !obs[7] * real(points[7])) + !obs[7] * imag(points[7]); label = 8;};\nborder C9 (t = !obs[8] * real(points[8]) + obs[8] * imag(points[8]), !obs[8] * real(points[9]) + obs[8] * imag(points[9])) {x = !obs[8] * t + obs[8] * real(points[8]); y = pentes[8] * (t - !obs[8] * real(points[8])) + !obs[8] * imag(points[8]); label = 9;};\nborder C10 (t = !obs[9] * real(points[9]) + obs[9] * imag(points[9]), !obs[9] * real(points[10]) + obs[9] * imag(points[9])) {x = !obs[9] * t + obs[9] * real(points[9]); y = pentes[9] * (t - !obs[9] * real(points[9])) + !obs[9] * imag(points[9]); label = 10;};\nborder C11 (t = !obs[10] * real(points[10]) + obs[10] * imag(points[10]), !obs[10] * real(points[11]) + obs[10] * imag(points[10])) {x = !obs[10] * t + obs[10] * real(points[10]); y = pentes[10] * (t - !obs[10] * real(points[10])) + !obs[10] * imag(points[10]); label = 11;};\nborder C12 (t = !obs[11] * real(points[11]) + obs[11] * imag(points[11]), !obs[11] * real(points[12]) + obs[11] * imag(points[11])) {x = !obs[11] * t + obs[11] * real(points[11]); y = pentes[11] * (t - !obs[11] * real(points[11])) + !obs[11] * imag(points[11]); label = 12;};\nborder C13 (t = !obs[12] * real(points[12]) + obs[12] * imag(points[12]), !obs[12] * real(points[13]) + obs[12] * imag(points[12])) {x = !obs[12] * t + obs[12] * real(points[12]); y = pentes[12] * (t - !obs[12] * real(points[12])) + !obs[12] * imag(points[12]); label = 13;};\nborder C14 (t = !obs[13] * real(points[13]) + obs[13] * imag(points[13]), !obs[13] * real(points[0]) + obs[13] * imag(points[13])) {x = !obs[13] * t + obs[13] * real(points[13]); y = pentes[13] * (t - !obs[13] * real(points[13])) + !obs[13] * imag(points[13]); label = 14;};\n\n//plot(C9(-n) + C10(-n) + C11(-n) + C12(-n) + C13(-n) + C14(n));\n\nmesh Th = buildmesh(C1(n) + C2(n) + C3(n) + C4(n) + C5(n) + C6(n) + C7(n) + C8(n) + C9(n) + C10(n) + C11(n) + C12(n) + C13(n) + C14(n), fixeborder = true);\n\n//plot(Th);\n\nfespace Vh(Th, P2, periodic = [[1, !cond[0] * (x + real(points[5] - points[1])) + cond[0] * (imag(points[5] - points[1]))], [6, !cond[0] * x + cond[0] * y],\n                                [3, !cond[1] * (x + real(points[7] - points[3])) + cond[1] * (imag(points[7] - points[3]))], [8, !cond[1] * x + cond[1] * y],\n                                [5, !cond[2] * (x + real(points[9] - points[5])) + cond[2] * (imag(points[9] - points[5]))], [10, !cond[2] * x + cond[2] * y],\n                                [7, !cond[3] * (x + real(points[11] - points[7])) + cond[3] * (imag(points[11] - points[7]))], [12, !cond[3] * x + cond[3] * y],\n                                [9, !cond[4] * (x + real(points[13] - points[9])) + cond[4] * (imag(points[13] - points[9]))], [14, !cond[4] * x + cond[4] * y],\n                                [11, !cond[5] * (x + real(points[1] - points[11])) + cond[5] * (imag(points[1] - points[11]))], [2, !cond[5] * x + cond[5] * y],\n                                [13, !cond[6] * (x + real(points[3] - points[13])) + cond[6] * (imag(points[3] - points[13]))], [4, !cond[6] * x + cond[6] * y]]);\n\n\nVh u1, u2;\nreal sigma = 0.00001;\nvarf op(u1, u2) = int2d(Th)(dx(u1) * dx(u2) + dy(u1) * dy(u2) - sigma * u1 * u2);\nvarf b([u1], [u2]) = int2d(Th)(u1 * u2);\nmatrix OP = op (Vh, Vh, solver = Crout, factorize = 1);\nmatrix B = b(Vh, Vh, solver = CG, eps = 1e-20);\nreal[int] ev(nbrVP);\nVh[int] eV(nbrVP);\nint k = EigenValue(OP, B, sym = true, sigma = sigma, value = ev, vector = eV, tol = 1e-10, maxit = 0, ncv = 0);\n\nreal aire = int2d(Th)(1);\n\nfor(int i = 0; i < k; i++)\n{\n    cout << ev[i] << " normalisation : " << ev[i] * aire << endl;\n\n    if (ev[i] * aire > 10e-2)\n    {\n        break;\n    }\n}'

        fichier = open(f"{chemin}{nom}.edp", "w")
        fichier.write(programme)
        fichier.close()

        if cheminFF != "":
            texte = os.popen(f"{cheminFF} {chemin}{nom}.edp").read()
        else:
            texte = os.popen(f"{chemin}{nom}.edp").read()

        if len(texte.split("Error in the definition of subdomain ")) == 1:
            break

        else:
            positionsChangements = []

            for pos in range(1, len(texte.split("Error in the definition of subdomain "))):
                positionsChangements.append(int(texte.split("Error in the definition of subdomain ")[pos].split("/")[0][-1]) - 1)

            for pos in positionsChangements:
                signeBords[pos] = (signeBords[pos] == "-") * "" + (signeBords[pos] == "") * "-"

    if len(texte.split("Missing edge")) > 1 or len(texte.split("Some giving point are outside the domain")) != 1 or len(texte.split("Exec error")) > 1:
        return 0
    
    return - float(texte.split("normalisation : ")[-1].split("\n")[0])

# Valeurs obtenues après déformations du 14gone de base.
val = [-0.05051446018243234, 0.5214320634823795, 0.539792418683307, 0.007147218278215177, 0.029083601131721447, 0.2479515729877344, 0.30184834621683265, 0.27430284766283064, -0.4771859414119586, 0.8978340422615234, 0.38720600446859565, 0.040375669637777424, 0.37132267833862775, 0.24960243022677797]

res = minimize(valeurPropre, val, method = 'nelder-mead', options={'xatol': 1e-6, 'disp': True})
print(res.x)