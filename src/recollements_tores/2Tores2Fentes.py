import os
import random

from scipy import optimize
from scipy.optimize import minimize

###########################################################
# Nombre de valeurs propres à calculer
nbrVP = 10

# Nombre de divisions par côté
precision = 45

signeBords = ["", "", "", "-", "", "", "-", "-", "", "-", "", "-", "", "-", "", "", "-", "-", "", "-", "", "-", "", "", "", "", "", ""]

# Nom du fichier .edp
nom = "2_tores_2_fentes"

# Chemin du fichier .edp
chemin = "path"

# Chemin de l'exécutable modifié, si besoin, de FreeFem++
cheminFF = "path"

def presqueNul(nbr, tol):
    if abs(nbr) < tol:
        return 0
    
    return nbr

def trouveIntersection(k, l, m, n, o, p):
    a = 0
    b = 0

    if m > max(k, 0):
            return False

    if k > 0:
        # Intervalle de superposition : [max(m, 0), min(k, o)]
        a = max(m, 0)
        b = min(k, o)

    elif k < 0:
        # Intervalle de superposition : [max(k, m), min(o, 0)]
        a = max(k, m)
        b = min(o, 0)

    else:
        return (m > 0) ^ (p < 0)

    return (round((l / k) * a, 6) > round(((p - n) / (o - m)) * a + n - ((p - n) / (o - m)) * m, 6)) ^ (round((l / k) * b, 6) > round(((p - n) / (o - m)) * b + n - ((p - n) / (o - m)) * m, 6)) or (round((l / k) * a, 6) == round(((p - n) / (o - m)) * a + n - ((p - n) / (o - m)) * m, 6)) ^ (round((l / k) * b, 6) == round(((p - n) / (o - m)) * b + n - ((p - n) / (o - m)) * m, 6))

def valeurPropre(params):
    bugDeSigne = 0

    i, j, k, l, m, n, o, p = params

    if (((i < 0 and l > (j / i) * k and l < (j / i) * (k - 1)) or (i > 0 and l < (j / i) * k and l > (j / i) * (k - 1))) and (k != 0 and abs((j / i) - (l / k)) > 10e-3)):
        if ((i < 0 and n > (j / i) * m and p < round((j / i) * (o - 1), 2)) or (i > 0 and n < (j / i) * m and p > (j / i) * (o - 1))) and (presqueNul(m - o, 10e-6) != 0 and not trouveIntersection(k, l, m, n, o, p)) and (presqueNul(m - k, 10e-6) != 0 and presqueNul(n - l, 10e-6) != 0 and presqueNul(1 + i - o, 10e-6) != 0 and presqueNul(k - o, 10e-6)):
            texte = ""

            while True:
                programme = f'int nbrValeursPropres = 10;\n\nreal tau1 = {presqueNul(i, 10e-6)};\nreal tau2 = {presqueNul(j, 10e-6)};\n\nreal sigma1 = {presqueNul(k, 10e-6)};\nreal sigma2 = {presqueNul(l, 10e-6)};\n\n// Point de base\nreal zeta1 = {presqueNul(m, 10e-6)};\nreal zeta2 = {presqueNul(n, 10e-6)};\n\n// Tête\nreal gamma1 = {presqueNul(o, 10e-6)};\nreal gamma2 = {presqueNul(p, 10e-6)};\n\nint n = {precision};' + '\n\nreal pente1 = 1;\n\nif (sigma1 != 0)\n{\n    pente1 = sigma2 / sigma1;\n}\n\nreal pente2 = 1;\n\nif (gamma1 - zeta1 != 0)\n{\n    pente2 = (gamma2 - zeta2) / (gamma1 - zeta1);\n}\n\nreal pente3 = 1;\n\nif (zeta1 - sigma1 != 0)\n{\n    pente3 = (zeta2 - sigma2) / (zeta1 - sigma1);\n}\n\nreal pente4 = 1;\n\nif (1 + tau1 - gamma1)\n{\n    pente4 = (tau2 - gamma2) / (1 + tau1 - gamma1);\n}\n\nreal pente5 = 1;\n\nif (tau1 != 0)\n{\n    pente5 = tau2 / tau1;\n}\n\nreal pente6 = (tau2 - gamma2 / 5) / (1 + tau1 - gamma1);\n\nreal tX = 1;\nreal tY = 0.2;\nreal pente = ((gamma2 - zeta2 * tY) / (gamma1 - zeta1 * tX));\n\nreal p = ((tau2 - zeta2 * tY) / (1 + tau1 - zeta1 * tX));\n\nborder TG1 (t = 0, tau1 + tau2 * (tau1 == 0)) {x = t * (tau1 != 0); y = pente5 * t; label = 1;};\nborder TD1 (t = (4 + tau1) * (tau1 != 0), (4 + 2 * tau1) * (tau1 != 0) + tau2 * (tau1 == 0)) {x = t * (tau1 != 0) + 2 * (tau1 == 0); y = pente5 * t + (tau2 - pente5 * (4 + 2 * tau1)) * (tau1 != 0); label = 2;};\n\nborder TC11 (t = 0, sigma1 + sigma2 * (sigma1 == 0)) {x = t * (sigma1 != 0); y = t * pente1; label = 3;};\nborder TC12 (t = sigma1 + (sigma2 - sigma1) * (zeta1 - sigma1 == 0), zeta1 + (zeta2 - zeta1) * (zeta1 - sigma1 == 0)) {x = t * (zeta1 - sigma1 != 0) + sigma1 * (zeta1 - sigma1 == 0); y = (t * pente3 + sigma2 - pente3 * sigma1) * (zeta1 - sigma1 != 0) + t * (zeta1 - sigma1 == 0); label = 4;};\nborder TC13 (t = zeta1 + (zeta2 - zeta1) * (gamma1 - zeta1 == 0), gamma1 + gamma2 * (gamma1 - zeta1 == 0)) {x = t * (gamma1 - zeta1 != 0) + zeta1 * (gamma1 - zeta1 == 0); y = (t * pente2 + zeta2 - pente2 * zeta1) * (gamma1 - zeta1 != 0) + t * (gamma1 - zeta1 == 0); label = 5;};\nborder TC14 (t = gamma1 + (gamma2 - gamma1) * (1 + tau1 - gamma1 == 0), (1 + tau1) + tau2 * (1 + tau1 - gamma1 == 0)) {x = t * (1 + tau1 - gamma1 != 0) + gamma1 * (1 + tau1 - gamma1 == 0); y = (t * pente4 + gamma2 - gamma1 * pente4) * (1 + tau1 - gamma1 != 0) + t * (1 + tau1 - gamma1 == 0); label = 6;};\n\nborder TC15 (t = gamma1, zeta1 * tX) {x = t; y = t * pente + gamma2 - pente * gamma1; label = 100;};\nborder TC16 (t = zeta1 * tX, (1 + tau1) + tau2 * (1 + tau1 - gamma1 == 0)) {x = t * (1 + tau1 - gamma1 != 0) + gamma1 * (1 + tau1 - gamma1 == 0); y = (t * p + zeta2 * tY - zeta1 * tX * p) * (1 + tau1 - gamma1 != 0) + t * (1 + tau1 - gamma1 == 0); label = 200;};\n\nborder TH1 (t = tau1, 1 + tau1) {x = t; y = tau2; label = 7;};\n\n\nborder TC21 (t = (3 + tau1) * (sigma1 != 0), 3 + tau1 + sigma1 + (sigma2 - 3 - tau1) * (sigma1 == 0)) {x = t * (sigma1 != 0) + (3 + tau1) * (sigma1 == 0); y = t * pente1 - (3 + tau1) * pente1 * (sigma1 != 0); label = 8;};\nborder TC22 (t = 3 + tau1 + sigma1 + (sigma2 - sigma1 - 3 + tau1) * (zeta1 - sigma1 == 0), 3 + tau1 + zeta1 + (zeta2 - zeta1 - 3 + tau1) * (zeta1 - sigma1 == 0)) {x = t * (zeta1 - sigma1 != 0) + (sigma1 + 3 + tau1) * (zeta1 - sigma1 == 0); y = (t * pente3 + sigma2 - pente3 * (sigma1 + 3 + tau1)) * (zeta1 - sigma1 != 0) + t * (zeta1 - sigma1 == 0); label = 9;};\nborder TC23 (t = 3 + tau1 + zeta1 + (zeta2 - zeta1 - 3 + tau1) * (gamma1 - zeta1 == 0), 3 + gamma1 + tau1 + (gamma2 - 3 + tau1) * (gamma1 - zeta1 == 0)) {x = t * (gamma1 - zeta1 != 0) + (zeta1 + 3 + tau1) * (gamma1 - zeta1 == 0); y = (t * pente2 + zeta2 - pente2 * (zeta1 + 3 + tau1)) * (gamma1 - zeta1 != 0) + t * (gamma1 - zeta1 == 0); label = 10;};\nborder TC24 (t = 3 + tau1 + gamma1 + (gamma2 - gamma1 - 3 + tau1) * (1 + tau1 - gamma1 == 0), (4 + 2 * tau1) + tau2 * (1 + tau1 - gamma1 == 0)) {x = t * (1 + tau1 - gamma1 != 0) + (gamma1 + 3 + tau1) * (1 + tau1 - gamma1 == 0); y = (t * pente4 + gamma2 - (gamma1 + 3 + tau1) * pente4) * (1 + tau1 - gamma1 != 0) + t * (1 + tau1 - gamma1 == 0); label = 11;};\n\nborder TC25 (t = 3 + tau1 + gamma1, 3 + tau1 + zeta1 * tX) {x = t; y = t * pente + gamma2 - pente * (gamma1 + 3 + tau1); label = 300;};\nborder TC26 (t = 3 + tau1 + zeta1 * tX, (4 + 2 * tau1) + tau2 * (1 + tau1 - gamma1 == 0)) {x = t * (1 + tau1 - gamma1 != 0) + gamma1 * (1 + tau1 - gamma1 == 0); y = (t * p + zeta2 * tY - (zeta1 * tX + 3 + tau1) * p) * (1 + tau1 - gamma1 != 0) + t * (1 + tau1 - gamma1 == 0); label = 400;};\n\nborder TB1 (t = 3 + tau1, 4 + tau1) {x = t; y = 0; label = 12;};\n\n\nborder TG2 (t = 6 + 3 * tau1, 6 + 4 * tau1 + tau2 * (tau1 == 0)) {x = t * (tau1 != 0); y = pente5 * t - (6 + 3 * tau1) * pente5; label = 13;};\nborder TD2 (t = (10 + 4 * tau1) * (tau1 != 0), (10 + 5 * tau1) * (tau1 != 0) + tau2 * (tau1 == 0)) {x = t * (tau1 != 0) + 3 * (tau1 == 0); y = pente5 * t + (tau2 - pente5 * (10 + 5 * tau1)) * (tau1 != 0); label = 14;};\n\nborder TC31 (t = (6 + 3 * tau1) * (sigma1 != 0), 6 + 3 * tau1 + sigma1 + (sigma2 - 6 - 3 * tau1) * (sigma1 == 0)) {x = t * (sigma1 != 0) + (6 + 3 * tau1) * (sigma1 == 0); y = t * pente1 - (6 + 3 * tau1) * pente1 * (sigma1 != 0); label = 15;};\nborder TC32 (t = 6 + 3 * tau1 + sigma1 + (sigma2 - sigma1 - 6 - 3 * tau1) * (zeta1 - sigma1 == 0), 6 + 3 * tau1 + zeta1 + (zeta2 - zeta1 - 6 - 3 * tau1) * (zeta1 - sigma1 == 0)) {x = t * (zeta1 - sigma1 != 0) + (sigma1 + 6 + 3 * tau1) * (zeta1 - sigma1 == 16); y = (t * pente3 + sigma2 - pente3 * (sigma1 + 6 + 3 * tau1)) * (zeta1 - sigma1 != 0) + t * (zeta1 - sigma1 == 0); label = 16;};\nborder TC33 (t = 6 + 3 * tau1 + zeta1 + (zeta2 - zeta1 - 6 - 3 * tau1) * (gamma1 - zeta1 == 0), 6 + 3 * tau1 + gamma1 + (gamma2 - 6 - 3 * tau1) * (gamma1 - zeta1 == 0)) {x = t * (gamma1 - zeta1 != 0) + (zeta1 + 6 + 3 * tau1) * (gamma1 - zeta1 == 0); y = (t * pente2 + zeta2 - pente2 * (zeta1 + 6 + 3 * tau1)) * (gamma1 - zeta1 != 0) + t * (gamma1 - zeta1 == 0); label = 17;};\nborder TC34 (t = 6 + 3 * tau1 + gamma1 + (gamma2 - gamma1 - 6 - 3 * tau1) * (1 + tau1 - gamma1 == 0), (7 + 4 * tau1) + tau2 * (1 + tau1 - gamma1 == 0)) {x = t * (1 + tau1 - gamma1 != 0) + (gamma1 + 6 + 3 * tau1) * (1 + tau1 - gamma1 == 0); y = (t * pente4 + gamma2 - (gamma1 + 6 + 3 * tau1) * pente4) * (1 + tau1 - gamma1 != 0) + t * (1 + tau1 - gamma1 == 0); label = 18;};\n\nborder TC35 (t = 6 + 3 * tau1 + gamma1, 6 + 3 * tau1 + zeta1 * tX) {x = t; y = t * pente + gamma2 - pente * (gamma1 + 6 + 3 * tau1); label = 500;};\nborder TC36 (t = 6 + 3 * tau1 + zeta1 * tX, (7 + 4 * tau1) + tau2 * (1 + tau1 - gamma1 == 0)) {x = t * (1 + tau1 - gamma1 != 0) + gamma1 * (1 + tau1 - gamma1 == 0); y = (t * p + zeta2 * tY - (zeta1 * tX + 6 + 3 * tau1) * p) * (1 + tau1 - gamma1 != 0) + t * (1 + tau1 - gamma1 == 0); label = 600;};\n\nborder TH2 (t = 6 + 4 * tau1, 7 + 4 * tau1) {x = t; y = tau2; label = 19;};\n\nborder TC41 (t = (9 + 4 * tau1) * (sigma1 != 0), 9 + 4 * tau1 + sigma1 + (sigma2 - 9 - 4 * tau1) * (sigma1 == 0)) {x = t * (sigma1 != 0) + (9 + 4 * tau1) * (sigma1 == 0); y = t * pente1 - (9 + 4 * tau1) * pente1 * (sigma1 != 0); label = 20;};\nborder TC42 (t = 9 + 4 * tau1 + sigma1 + (sigma2 - sigma1 - 9 - 4 * tau1) * (zeta1 - sigma1 == 0), 9 + 4 * tau1 + zeta1 + (zeta2 - zeta1 - 9 - 4 * tau1) * (zeta1 - sigma1 == 0)) {x = t * (zeta1 - sigma1 != 0) + (sigma1 + 9 + 4 * tau1) * (zeta1 - sigma1 == 0); y = (t * pente3 + sigma2 - pente3 * (sigma1 + 9 + 4 * tau1)) * (zeta1 - sigma1 != 0) + t * (zeta1 - sigma1 == 0); label = 21;};\nborder TC43 (t = 9 + 4 * tau1 + zeta1 + (zeta2 - zeta1 - 9 - 4 * tau1) * (gamma1 - zeta1 == 0), 9 + gamma1 + 4 * tau1 + (gamma2 - 9 - 4 * tau1) * (gamma1 - zeta1 == 0)) {x = t * (gamma1 - zeta1 != 0) + (zeta1 + 9 + 4 * tau1) * (gamma1 - zeta1 == 0); y = (t * pente2 + zeta2 - pente2 * (zeta1 + 9 + 4 * tau1)) * (gamma1 - zeta1 != 0) + t * (gamma1 - zeta1 == 0); label = 22;};\nborder TC44 (t = 9 + 4 * tau1 + gamma1 + (gamma2 - gamma1 - 9 - 4 * tau1) * (1 + tau1 - gamma1 == 0), (10 + 5 * tau1) + tau2 * (1 + tau1 - gamma1 == 0)) {x = t * (1 + tau1 - gamma1 != 0) + (gamma1 + 9 + 4 * tau1) * (1 + tau1 - gamma1 == 0); y = (t * pente4 + gamma2 - (gamma1 + 9 + 4 * tau1) * pente4) * (1 + tau1 - gamma1 != 0) + t * (1 + tau1 - gamma1 == 0); label = 23;};\n\nborder TC45 (t = 9 + 4 * tau1 + gamma1, 6 + 4 * tau1 + zeta1 * tX) {x = t; y = t * pente + gamma2 - pente * (gamma1 + 9 + 4 * tau1); label = 700;};\nborder TC46 (t = 9 + 4 * tau1 + zeta1 * tX, (10 + 5 * tau1) + tau2 * (1 + tau1 - gamma1 == 0)) {x = t * (1 + tau1 - gamma1 != 0) + gamma1 * (1 + tau1 - gamma1 == 0); y = (t * p + zeta2 * tY - (zeta1 * tX + 9 + 4 * tau1) * p) * (1 + tau1 - gamma1 != 0) + t * (1 + tau1 - gamma1 == 0); label = 800;};\n\nborder TB2 (t = 9 + 4 * tau1, 10 + 4 * tau1) {x = t; y = 0; label = 24;};\n\nif (gamma2 > gamma1 * pente3 + sigma2 - pente3 * sigma1 || gamma1 > zeta1)\n{\n    //plot(TD1(n) + TC21(n) + TC22(n) + TC23(n) + TC24(n) + TB1(n) + TG2(n) + TC31(n) + TC32(n) + TC33(n) + TC34(n) + TH2(n));// + TD2(n) + TC41(n) + TC42(n) + TC43(n) + TC44(n) + TB2(n));\n    //plot(TG1(n) + TC11(n) + TC12(n) + TC13(n) + TC14(n) + TH1(n));// + TD1(n) + TC21(n) + TC22(n) + TC23(n) + TC24(n) + TB1(n));\n    ' + f'mesh Th = buildmesh(TG1({signeBords[0]}n) + TC11({signeBords[1]}n) + TC12({signeBords[2]}n) + TC13({signeBords[3]}n) + TC14({signeBords[4]}n) + TH1({signeBords[5]}n) + TD1({signeBords[6]}n) + TC21({signeBords[7]}n) + TC22({signeBords[8]}n) + TC23({signeBords[9]}n) + TC24({signeBords[10]}n) + TB1({signeBords[11]}n) + TG2({signeBords[12]}n) + TC31({signeBords[13]}n) + TC32({signeBords[14]}n) + TC33({signeBords[15]}n) + TC34({signeBords[16]}n) + TH2({signeBords[17]}n) + TD2({signeBords[18]}n) + TC41({signeBords[19]}n) + TC42({signeBords[20]}n) + TC43({signeBords[21]}n) + TC44({signeBords[22]}n) + TB2({signeBords[23]}n), fixedborder = 1);' + '\n\n    fespace Vh(Th, P1, periodic = [[1, y], [2, y], [3, y], [20, y], [4, y], [9, y], [5, y], [22, y], [6, y], [11, y],\n                                    [7, x + 3], [12, x], [13, y], [14, y], [15, y], [8, y], [16, y], [21, y],\n                                    [17, y], [10, y], [18, y], [23, y], [19, x + 3], [24, x]]);\n    \n    Vh u1, u2;\n    real sigma = 0.00001;\n\n    varf op (u1, u2) = int2d(Th)(dx(u1) * dx(u2) + dy(u1) * dy(u2) - sigma * u1 * u2);\n    varf bm ([u1], [u2]) = int2d(Th)(u1 * u2);\n\n    matrix OP = op(Vh, Vh, solver = Crout, factorize = 1);\n    matrix B = bm(Vh, Vh, solver = CG, eps = 1e-20);\n\n    real[int] ev(nbrValeursPropres);\n    Vh[int] eV(nbrValeursPropres);\n\n    int k = EigenValue(OP, B, sym = true, sigma = sigma, value = ev, vector = eV, tol=1e-10, maxit = 0, ncv = 0);\n\n    for (int i = 0; i < k; i++)\n    {\n        u1 = eV[i];\n\n        cout << "lambda[" << i << "] = " << ev[i] << ", normalisation : " << ev[i] * (2 * abs(tau2)) << ", err = " << int2d(Th)(dx(u1) * dx(u1) + dy(u1) * dy(u1) - (ev[i]) * u1 * u1) << endl;\n    }\n}\n\nelse\n{\n    //plot(TG1(n) + TC11(n) + TC12(n) + TC13(n) + TC15(n) + TC16(n) + TH1(n) + TD1(n) + TC21(n) + TC22(n) + TC23(n) + TC25(n) + TC26(n) + TB1(n) + TG2(n) + TC31(n) + TC32(n) + TC33(n) + TC35(n) + TC36(n) + TH2(n) + TD2(n) + TC41(n) + TC42(n) + TC43(n) + TC45(n) + TC46(n) + TB2(n));\n    //plot(TG1(n) + TC11(n) + TC12(n) + TC13(n) + TC14(n) + TH1(n) + TC31(n) + TC32(n) + TC33(n) + TC35(n) + TC36(n) + TH2(n));// + TD1(n) + TC21(n) + TC22(n) + TC23(n) + TC25(n) + TC26(n) + TB1(n))\n    ' + f'mesh Th = buildmesh(TG1({signeBords[0]}n) + TC11({signeBords[1]}n) + TC12({signeBords[2]}n) + TC13({signeBords[3]}n) + TC15({signeBords[4]}n) + TC16({signeBords[5]}n) + TH1({signeBords[6]}n) + TD1({signeBords[7]}n) + TC21({signeBords[8]}n) + TC22({signeBords[9]}n) + TC23({signeBords[10]}n) + TC25({signeBords[11]}n) + TC26({signeBords[12]}n) + TB1({signeBords[13]}n) + TG2({signeBords[14]}n) + TC31({signeBords[15]}n) + TC32({signeBords[16]}n) + TC33({signeBords[17]}n) + TC35({signeBords[18]}n) + TC36({signeBords[19]}n) + TH2({signeBords[20]}n) + TD2({signeBords[21]}n) + TC41({signeBords[22]}n) + TC42({signeBords[23]}n) + TC43({signeBords[24]}n) + TC45({signeBords[25]}n) + TC46({signeBords[26]}n) + TB2({signeBords[27]}n), fixedborder = 1);' + '\n\n    fespace Vh(Th, P1, periodic = [[1, y], [2, y], [3, y], [20, y], [4, y], [9, y], [5, y], [22, y], [100, y], [300, y],\n                                    [200, y], [400, y], [7, x + 3], [12, x], [13, y], [14, y], [15, y], [8, y], [16, y], [21, y],\n                                    [17, y], [10, y], [500, y], [700, y], [19, x + 3], [24, x], [600, y], [800, y]]);\n\n    Vh u1, u2;\n    real sigma = 0.00001;\n\n    varf op (u1, u2) = int2d(Th)(dx(u1) * dx(u2) + dy(u1) * dy(u2) - sigma * u1 * u2);\n    varf bm ([u1], [u2]) = int2d(Th)(u1 * u2);\n\n    matrix OP = op(Vh, Vh, solver = Crout, factorize = 1);\n    matrix B = bm(Vh, Vh, solver = CG, eps = 1e-20);\n\n    real[int] ev(nbrValeursPropres);\n    Vh[int] eV(nbrValeursPropres);\n\n    int k = EigenValue(OP, B, sym = true, sigma = sigma, value = ev, vector = eV, tol=1e-10, maxit = 0, ncv = 0);\n\n    for (int i = 0; i < k; i++)\n    {\n        u1 = eV[i];\n\n        cout << "lambda[" << i << "] = " << ev[i] << ", normalisation : " << ev[i] * (2 * abs(tau2)) << ", err = " << int2d(Th)(dx(u1) * dx(u1) + dy(u1) * dy(u1) - (ev[i]) * u1 * u1) << endl;\n    }\n}'

                fichier = open(f"{chemin}{nom}.edp", "w")
                fichier.write(programme)
                fichier.close()

                if cheminFF != "":
                    texte = os.popen(f"{cheminFF} {chemin}{nom}.edp").read()
                else:
                    texte = os.popen(f"{chemin}{nom}.edp").read()

                if len(texte.split("Error in the definition of subdomain ")) == 1 or bugDeSigne == 2:
                    break

                # Il est fréquent que l'orientation des bords soient mal choisies.
                # Dans ce cas, une erreur de compilation est levée et on peut
                # récupérer de l'information pour déterminer quels bords sont
                # mal orientés.
                else:
                    positionsChangements = []
                    bugDeSigne += 1

                    for pos in range(1, len(texte.split("Error in the definition of subdomain "))):
                        positionsChangements.append(int(texte.split("Error in the definition of subdomain ")[pos].split("/")[0].split("border ")[1]) - 1)

                    for pos in positionsChangements:
                        signeBords[pos] = (signeBords[pos] == "-") * "" + (signeBords[pos] == "") * "-"
                                
            if len(texte.split("tau1 trop proche de sigma1")) == 3 or len(texte.split("Some giving point are outside the domain")) != 1 or len(texte.split("The boundary is crossing maybe!")) != 1 or bugDeSigne == 2 or len(texte.split("Number of Negative triangles")) != 1 or len(texte.split("normalisation : ")) < 5:
                return 0

            return - float(texte.split("normalisation : ")[4].split(",")[0])
        
    return 0

for i in range(50):
    res = minimize(valeurPropre, [random.random(), random.random(), random.random(), random.random(), random.random(), random.random(), random.random(), random.random()], method = 'nelder-mead', options={'xatol': 1e-6, 'disp': False})

    val = abs(valeurPropre(list(res.x)))

    if (val > 0):
        print(f"Test #{i} :")
        print(f"{val} : {list(res.x)}")
###########################################################