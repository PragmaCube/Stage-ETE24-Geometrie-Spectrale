import os
import random

from scipy.optimize import minimize

###########################################################
# Nombre de valeurs propres à calculer
nbrVP = 10

# Nombre de divisions par côté
precision = 45

parametres = []
valeursPropresNormalisees = []

# signeBords = ["", "", "", "-", "", "", "-", "-", "", "-", "", "-", "", "-", "", "", "-", "-", "", "-", "", "-", "", ""]
signeBords = ["-", "-", "-", "-", "", "-", "-", "", "", "-", "", "-", "-", "", "", "-", "", "-", "", ""]

# Nom du fichier .edp
nom = "tore3Fentes"

# Chemin du fichier .edp
chemin = "path"

# Chemin de l'exécutable modifié, si besoin, de FreeFem++
cheminFF = "path"

def presqueNul(nbr, tol):
    if abs(nbr) < tol:
        return 0
    
    return nbr

def valeurPropre(params):
    i, j, k, l = params

    if ((((i < 0 and l < (j / i) * (k - 1 / 3) and l > (j / i) * k) or (i > 0 and l > (j / i) * (k - 1 / 3)) and l < (j / i) * k) and ((k != 0 and abs((j / i) - (l / k)) > 10e-6 and (l / k) < (j / i)))) or (presqueNul(i, 10e-6) == 0)) and abs(j - l) > 10e-6:
        texte = ""

        while True:
            programme = f'int nbrValeursPropres = 10;\n\nreal tau1 = {presqueNul(i, 10e-6)};\nreal tau2 = {j};\n\nreal sigma1 = {presqueNul(k, 10e-6)};\nreal sigma2 = {l};' + '\n\nif (abs(sigma1 - tau1) > 10e-6)\n{\n    real pente1 = 1;\n\n    if (sigma1 != 0)\n    {\n        pente1 = sigma2 / sigma1;\n    }\n\n    real pente2 = (sigma2 - tau2) / (sigma1 - tau1);\n\n    real pente3 = 1;\n\n    if (tau1 != 0)\n    {\n        pente3 = tau2 / tau1;\n    }\n\n    ' + f'int n = {precision};' + '\n\n    border TG (t = 0, tau1 + tau2 * (tau1 == 0)) {x = t * (tau1 != 0); y = pente3 * t; label = 100;};\n    border TD (t = 2.5 * (tau1 != 0), (2.5 + tau1) * (tau1 != 0) + tau2 * (tau1 == 0)) {x = t * (tau1 != 0) + 2.5 * (tau1 == 0); y = pente3 * t + (tau2 - pente3 * (2.5 + tau1)) * (tau1 != 0); label = 200;};\n\n    border TB2 (t = 0.5, real(1) / real(3) + 0.5) {x = t; y = 0; label = 1;};\n    border TB3 (t = real(1) / real(3) + 1, real(2) / real(3) + 1) {x = t; y = 0; label = 2;};\n    border TB4 (t = real(2) / real(3) + 1.5, 2.5) {x = t; y = 0; label = 3;};\n\n    border TH2 (t = tau1 + 0.5, tau1 + real(1) / real(3) + 0.5) {x = t; y = tau2; label = 4;};\n    border TH3 (t = tau1 + 1 + real(1) / real(3), tau1 + real(2) / real(3) + 1) {x = t; y = tau2; label = 5;};\n    border TH4 (t = real(2) / real(3) + 1.5 + tau1, 2.5 + tau1) {x = t; y = tau2; label = 6;};\n\n    border TC11 (t = 0, sigma1 + sigma2 * (sigma1 == 0)) {x = t * (sigma1 != 0); y = pente1 * t; label = 7;};\n    border TC12 (t = min(tau1, sigma1), max(tau1, sigma1)) {x = t; y = pente2 * t + tau2 - pente2 * tau1; label = 8;};\n\n    border TC21 (t = 0.5 * (sigma1 != 0), (sigma1 + 0.5) * (sigma1 != 0) + sigma2 * (sigma1 == 0)) {x = t * (sigma1 != 0) + 0.5 * (sigma1 == 0); y = pente1 * t + (sigma2 - pente1 * (sigma1 + 0.5)) * (sigma1 != 0); label = 9;};\n    border TC22 (t = min(0.5 + tau1, sigma1 + 0.5), max(0.5 + tau1, sigma1 + 0.5)) {x = t; y = pente2 * t + tau2 - pente2 * (tau1 + 0.5); label = 10;};\n    border TC23 (t = (0.5 + real(1) / real(3)) * (sigma1 != 0), (sigma1 + 0.5 + real(1) / real(3)) * (sigma1 != 0) + sigma2 * (sigma1 == 0)) {x = t * (sigma1 != 0) + (0.5 + real(1) / real(3)) * (sigma1 == 0); y = pente1 * t + (sigma2 - pente1 * (sigma1 + 0.5 + real(1) / real(3))) * (sigma1 != 0); label = 11;};\n    border TC24 (t = min(0.5 + real(1) / real(3) + tau1, sigma1 + 0.5 + real(1) / real(3)), max(0.5 + real(1) / real(3) + tau1, sigma1 + 0.5 + real(1) / real(3))) {x = t; y = pente2 * t + tau2 - pente2 * (tau1 + 0.5 + real(1) / real(3)); label = 12;};\n\n    border TC31 (t = (1 + real(1) / real(3)) * (sigma1 != 0), (1 + real(1) / real(3) + sigma1) * (sigma1 != 0) + sigma2 * (sigma1 == 0)) {x = t * (sigma1 != 0) + (1 + real(1) / real(3)) * (sigma1 == 0); y = pente1 * t + (sigma2 - pente1 * (sigma1 + 1 + real(1) / real(3))) * (sigma1 != 0); label = 13;};\n    border TC32 (t = 1 + real(1) / real(3) + tau1, 1 + real(1) / real(3) + sigma1) {x = t; y = pente2 * t + tau2 - pente2 * (tau1 + 1 + real(1) / real(3)); label = 14;};\n    border TC33 (t = (1 + real(2) / real(3)) * (sigma1 != 0), (1 + real(2) / real(3) + sigma1) * (sigma1 != 0) + sigma2 * (sigma1 == 0)) {x = t * (sigma1 != 0) + (1 + real(2) / real(3)) * (sigma1 == 0); y = pente1 * t + (sigma2 - pente1 * (sigma1 + 1 + real(2) / real(3))) * (sigma1 != 0); label = 15;};\n    border TC34 (t = 1 + real(2) / real(3) + tau1, 1 + real(2) / real(3) + sigma1) {x = t; y = pente2 * t + tau2 - pente2 * (tau1 + 1 + real(2) / real(3)); label = 16;};\n\n    border TC41 (t = (1.5 + real(2) / real(3)) * (sigma1 != 0), (1.5 + real(2) / real(3) + sigma1) * (sigma1 != 0) + sigma2 * (sigma1 == 0)) {x = t * (sigma1 != 0) + (1.5 + real(2) / real(3)) * (sigma1 == 0); y = pente1 * t + (sigma2 - pente1 * (sigma1 + 1.5 + real(2) / real(3))) * (sigma1 != 0); label = 17;};\n    border TC42 (t = 1.5 + real(2) / real(3) + tau1, 1.5 + real(2) / real(3) + sigma1) {x = t; y = pente2 * t + tau2 - pente2 * (tau1 + 1.5 + real(2) / real(3)); label = 18;};\n\n    //plot(TH1(n));\n    //plot(TC11(n) + TC12(n) + TG(n));// + TC21(n) + TC22(n) + TH2(n) + TC24(n) + TC23(n) + TB2(n) + TC31(n) + TC32(n) + TH3(n) + TC34(n) + TC33(n) + TB3(n) + TC41(n) + TC42(n) + TH4(n) + TD(n) + TB4(n));\n\n   ' + f'mesh Th = buildmesh(TC11({signeBords[0]}n) + TC12({signeBords[1]}n) + TG({signeBords[2]}n) + TC21({signeBords[3]}n) + TC22({signeBords[4]}n) + TH2({signeBords[5]}n) + TC24({signeBords[6]}n) + TC23({signeBords[7]}n) + TB2({signeBords[8]}n) + TC31({signeBords[9]}n) + TC32({signeBords[10]}n) + TH3({signeBords[11]}n) + TC34({signeBords[12]}n) + TC33({signeBords[13]}n) + TB3({signeBords[14]}n) + TC41({signeBords[15]}n) + TC42({signeBords[16]}n) + TH4({signeBords[17]}n) + TD({signeBords[18]}n) + TB4({signeBords[19]}n), fixedborder = 1);' + '\n\n    fespace Vh(Th, P1, periodic = [[100, y], [200, y], [1, x + tau1], [4, x], [2, x + tau1], [5, x], [3, x + tau1], [6, x],\n                                    [8, y], [10, y], [12, y], [14, y], [16, y], [18, y], [7, y], [13, y], [11, y], [17, y],\n                                    [15, y], [9, y]]);\n\n    Vh u1, u2;\n    real sigma = 0.00001;\n\n    varf op (u1, u2) = int2d(Th)(dx(u1) * dx(u2) + dy(u1) * dy(u2) - sigma * u1 * u2);\n    varf bm ([u1], [u2]) = int2d(Th)(u1 * u2);\n\n    matrix OP = op(Vh, Vh, solver = Crout, factorize = 1);\n    matrix B = bm(Vh, Vh, solver = CG, eps = 1e-20);\n\n    real[int] ev(nbrValeursPropres);\n    Vh[int] eV(nbrValeursPropres);\n\n    int k = EigenValue(OP, B, sym = true, sigma = sigma, value = ev, vector = eV, tol=1e-10, maxit = 0, ncv = 0);\n\n    for (int i = 0; i < k; i++)\n    {\n        u1 = eV[i];\n\n        cout << "lambda[" << i << "] = " << ev[i] << ", normalisation : " << ev[i] * (abs(tau2)) << ", err = " << int2d(Th)(dx(u1) * dx(u1) + dy(u1) * dy(u1) - (ev[i]) * u1 * u1) << endl;\n    }\n}\n\nelse\n{\n    cout << "tau1 trop proche de sigma1";\n}'

            fichier = open(f"{chemin}{nom}.edp", "w")
            fichier.write(programme)
            fichier.close()

            if cheminFF != "":
                texte = os.popen(f"{cheminFF} {chemin}{nom}.edp").read()
            else:
                texte = os.popen(f"{chemin}{nom}.edp").read()

            if len(texte.split("Error in the definition of subdomain ")) == 1:
                break

            # Il est fréquent que l'orientation des bords soient mal choisies.
            # Dans ce cas, une erreur de compilation est levée et on peut
            # récupérer de l'information pour déterminer quels bords sont
            # mal orientés.
            else:
                positionsChangements = []

                for pos in range(1, len(texte.split("Error in the definition of subdomain "))):
                    positionsChangements.append(int(texte.split("Error in the definition of subdomain ")[pos].split("/")[0][-1]) - 1)

                for pos in positionsChangements:
                    signeBords[pos] = (signeBords[pos] == "-") * "" + (signeBords[pos] == "") * "-"
                            
        if len(texte.split("tau1 trop proche de sigma1")) == 3 or len(texte.split("Some giving point are outside the domain")) != 1  or len(texte.split("normalisation : ")) < 4:
            return 0

        return - float(texte.split("normalisation : ")[3].split(",")[0])
    
    return 0

for i in range(100):
    res = minimize(valeurPropre, [random.random(), random.random(), random.random(), random.random()], method = 'nelder-mead', options={'xatol': 1e-6, 'disp': True})

    val = abs(valeurPropre(list(res.x)))

    if (val > 0):
        print(f"Test #{i} :")
        print(f"{val} : {list(res.x)}")
###########################################################