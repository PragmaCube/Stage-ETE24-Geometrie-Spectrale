import numpy as np
import math as mt
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

import os

###########################################################
# Calculs pour déterminer les points du triangles et l'équation
# du cercle (côté courbé).

class SpectreTrianglesHyperboliques:
    def __init__(self):
        self.angleTriangle = 0
        self.p1 = 0
        self.p2 = 0
        self.p3 = 0
        self.p4 = 0
        self.p5 = 0
        self.pTriangle1 = [0, 0]
        self.pTriangle2 = [0, 0]
        self.r1 = 0
        self.donneesCercle = [0, 0, 0]
        pass
    
    def setAngle(self, angle):
        self.angleTriangle = angle

    def equationsP1(self, x):
        return [(x[0] - np.cos(self.angleTriangle)) ** 2 + x[1] ** 2 - np.cos(self.angleTriangle) ** 2, x[0] - 1]

    def equationP3(self, x):
        return [(x[0] - self.p2[0]) * self.p2[0] + (x[1] - self.p2[1]) * self.p2[1], x[1]]

    def equationPT1(self, x):
        return [(x[0] - self.p3[0]) ** 2 + self.p3[1] ** 2 - self.r1 ** 2, x[1]]

    def equationP4(self, x):
         return [x[0] ** 2 + x[1] ** 2 - 1, x[0] - self.pTriangle1[0]]

    def equationP5(self, x):
        return [(x[0] - self.p4[0]) * self.p4[0] + (x[1] - self.p4[1]) * self.p4[1], x[1]]

    def cercle(self, x):
        return [(self.pTriangle1[0] - x[0]) ** 2 + (self.pTriangle1[1] - x[1]) ** 2 - x[2] ** 2, (self.pTriangle2[0] - x[0]) ** 2 + (self.pTriangle2[1] - x[1]) ** 2 - x[2] ** 2, (self.p5[0] - x[0]) ** 2 + (self.p5[1] - x[1]) ** 2 - x[2] ** 2]

    def update(self):
        self.p1 = fsolve(self.equationsP1, [1, 1])

        if self.p1[1] > 0:
            self.p1[1] *= - 1

        self.angle2 = np.arcsin((1 - np.cos(self.angleTriangle)) / np.cos(self.angleTriangle))
        self.p2 = [np.cos(self.angle2), - np.sin(self.angle2)]
        self.p3 = fsolve(self.equationP3, [1, 1])
        self.r1 = mt.sqrt((self.p3[0] - self.p2[0]) ** 2 + (self.p3[1] - self.p2[1]) ** 2)
        self.pTriangle1 = fsolve(self.equationPT1, [1, 0])
        self.pTriangle2 = [np.cos(self.angleTriangle) * self.pTriangle1[0] - np.sin(self.angleTriangle) * self.pTriangle1[1], np.sin(self.angleTriangle) * self.pTriangle1[0] + np.cos(self.angleTriangle) * self.pTriangle1[1]]
        self.p4 = fsolve(self.equationP4, [1, 1])

        if self.p4[1] > 0:
            self.p4[1] *= - 1

        self.p5 = fsolve(self.equationP5, [1, 1])
        self.donneesCercle = fsolve(self.cercle, [1, 1, 1])

        if self.donneesCercle[2] < 0:
            self.donneesCercle[2] *= - 1

# donneesCercle contient le centre et le rayon du cercle qui sert de 3 côté
###########################################################

###########################################################
# Nombre de valeurs propres à calculer
nbrVP = 10

# Nombre de division par côté
precision = 40

angleMin = mt.pi / 100
angleMax = mt.pi / 3 - mt.pi / 100
pas = mt.pi / 200

instance = SpectreTrianglesHyperboliques()

angles = []
valeursPropresNormalisees = []

for i in np.arange(angleMin, angleMax, pas):
    instance.setAngle(i)
    instance.update()

    # Programme FreeFem++
    programme = f'int nbrValeursPropres = {nbrVP};\n\nint n = {precision};\n\nreal angle1 = {instance.angleTriangle};\n\nreal[int] pTriangle1(2);\npTriangle1[0] = {instance.pTriangle1[0]};\npTriangle1[1] = {instance.pTriangle1[1]};\n\nreal[int] pTriangle2(2);\npTriangle2[0] = {instance.pTriangle2[0]};\npTriangle2[1] = {instance.pTriangle2[1]};' + '\n\nreal[int] ordonneesOrigines(14);\n\nfor (int i = 0; i < 14; i++)\n{\n    ordonneesOrigines[i] = pTriangle2[1] - (pTriangle2[1] / pTriangle2[0]) * (2 * i + pTriangle2[0]);\n} \n\nreal p = pTriangle2[1] / pTriangle2[0];\n\n' + f'real a = {instance.donneesCercle[0]};\nreal b = {instance.donneesCercle[1]};\nreal r = {instance.donneesCercle[2]}' + ';\n\nborder B1 (t = 0, pTriangle1[0]) {x = t; y = 0; label = 1;};\nborder B2 (t = 2, 2 + pTriangle1[0]) {x = t; y = 0; label = 2;};\nborder B3 (t = 4, 4 + pTriangle1[0]) {x = t; y = 0; label = 3;};\nborder B4 (t = 6, 6 + pTriangle1[0]) {x = t; y = 0; label = 4;};\nborder B5 (t = 8, 8 + pTriangle1[0]) {x = t; y = 0; label = 5;};\nborder B6 (t = 10, 10 + pTriangle1[0]) {x = t; y = 0; label = 6;};\nborder B7 (t = 12, 12 + pTriangle1[0]) {x = t; y = 0; label = 7;};\nborder B8 (t = 14, 14 + pTriangle1[0]) {x = t; y = 0; label = 8;};\nborder B9 (t = 16, 16 + pTriangle1[0]) {x = t; y = 0; label = 9;};\nborder B10 (t = 18, 18 + pTriangle1[0]) {x = t; y = 0; label = 10;};\nborder B11 (t = 20, 20 + pTriangle1[0]) {x = t; y = 0; label = 11;};\nborder B12 (t = 22, 22 + pTriangle1[0]) {x = t; y = 0; label = 12;};\nborder B13 (t = 24, 24 + pTriangle1[0]) {x = t; y = 0; label = 13;};\nborder B14 (t = 26, 26 + pTriangle1[0]) {x = t; y = 0; label = 14;};\n\nborder D1 (t = 0, pTriangle2[0]) {x = t; y = p * t + ordonneesOrigines[0]; label = 15;};\nborder D2 (t = 2, 2 + pTriangle2[0]) {x = t; y = p * t + ordonneesOrigines[1]; label = 16;};\nborder D3 (t = 4, 4 + pTriangle2[0]) {x = t; y = p * t + ordonneesOrigines[2]; label = 17;};\nborder D4 (t = 6, 6 + pTriangle2[0]) {x = t; y = p * t + ordonneesOrigines[3]; label = 18;};\nborder D5 (t = 8, 8 + pTriangle2[0]) {x = t; y = p * t + ordonneesOrigines[4]; label = 19;};\nborder D6 (t = 10, 10 + pTriangle2[0]) {x = t; y = p * t + ordonneesOrigines[5]; label = 20;};\nborder D7 (t = 12, 12 + pTriangle2[0]) {x = t; y = p * t + ordonneesOrigines[6]; label = 21;};\nborder D8 (t = 14, 14 + pTriangle2[0]) {x = t; y = p * t + ordonneesOrigines[7]; label = 22;};\nborder D9 (t = 16, 16 + pTriangle2[0]) {x = t; y = p * t + ordonneesOrigines[8]; label = 23;};\nborder D10 (t = 18, 18 + pTriangle2[0]) {x = t; y = p * t + ordonneesOrigines[9]; label = 24;};\nborder D11 (t = 20, 20 + pTriangle2[0]) {x = t; y = p * t + ordonneesOrigines[10]; label = 25;};\nborder D12 (t = 22, 22 + pTriangle2[0]) {x = t; y = p * t + ordonneesOrigines[11]; label = 26;};\nborder D13 (t = 24, 24 + pTriangle2[0]) {x = t; y = p * t + ordonneesOrigines[12]; label = 27;};\nborder D14 (t = 26, 26 + pTriangle2[0]) {x = t; y = p * t + ordonneesOrigines[13]; label = 28;};\n\nreal angle2 = atan((b - pTriangle1[1]) / (a - pTriangle1[0])) + pi;\nreal angle3 = atan((b - pTriangle2[1]) / (a - pTriangle2[0])) + pi;\n\nborder C1 (t = angle2, angle3) {x = a + r * cos(t); y = b + r * sin(t); label = 29;}\nborder C2 (t = angle2, angle3) {x = 2 + a + r * cos(t); y = b + r * sin(t); label = 30;}\nborder C3 (t = angle2, angle3) {x = 4 + a + r * cos(t); y = b + r * sin(t); label = 31;}\nborder C4 (t = angle2, angle3) {x = 6 + a + r * cos(t); y = b + r * sin(t); label = 32;}\nborder C5 (t = angle2, angle3) {x = 8 + a + r * cos(t); y = b + r * sin(t); label = 33;}\nborder C6 (t = angle2, angle3) {x = 10 + a + r * cos(t); y = b + r * sin(t); label = 34;}\nborder C7 (t = angle2, angle3) {x = 12 + a + r * cos(t); y = b + r * sin(t); label = 35;}\nborder C8 (t = angle2, angle3) {x = 14 + a + r * cos(t); y = b + r * sin(t); label = 36;}\nborder C9 (t = angle2, angle3) {x = 16 + a + r * cos(t); y = b + r * sin(t); label = 37;}\nborder C10 (t = angle2, angle3) {x = 18 + a + r * cos(t); y = b + r * sin(t); label = 38;}\nborder C11 (t = angle2, angle3) {x = 20 + a + r * cos(t); y = b + r * sin(t); label = 39;}\nborder C12 (t = angle2, angle3) {x = 22 + a + r * cos(t); y = b + r * sin(t); label = 40;}\nborder C13 (t = angle2, angle3) {x = 24 + a + r * cos(t); y = b + r * sin(t); label = 41;}\nborder C14 (t = angle2, angle3) {x = 26 + a + r * cos(t); y = b + r * sin(t); label = 42;}\n\nmesh Th = buildmesh(B1(n) + C1(n) + D1(-n) + B2(n) + C2(n) + D2(-n) + B3(n) + C3(n) + D3(-n) + B4(n) + C4(n) + D4(-n) + B5(n) + C5(n) + D5(-n) + B6(n) + C6(n) + D6(-n) + B7(n) + C7(n) + D7(-n) + B8(n) + C8(n) + D8(-n) + B9(n) + C9(n) + D9(-n) + B10(n) + C10(n) + D10(-n) + B11(n) + C11(n) + D11(-n) + B12(n) + C12(n) + D12(-n) + B13(n) + C13(n) + D13(-n) + B14(n) + C14(n) + D14(-n));\n\nfespace Vh(Th, P2, periodic = [[2, x + 2], [3, x], [4, x + 2], [5, x], [6, x + 2], [7, x], [8, x + 2], [9, x], [10, x + 2], [11, x],\n                             [12, x + 2], [13, x], [14, x - 26], [1, x] , [15, y], [16, y], [17, y], [18, y], [19, y], [20, y],\n                             [21, y], [22, y], [23, y], [24, y], [25, y], [26, y], [27, y], [28, y], [29, y], [34, y],\n                             [31, y], [36, y], [33, y], [38, y], [35, y], [40, y], [37, y], [42, y], [39, y],\n                             [30, y], [41, y], [32, y]]);\n\nVh u1, u2;\nreal sigma = 0.00001;\n\nvarf op (u1, u2) = int2d(Th)(dx(u1) * dx(u2) + dy(u1) * dy(u2) - 4 * sigma * u1 * u2 / (1 - (x - 1 * floor(x))^2 - y^2)^2);\nvarf bm ([u1], [u2]) = int2d(Th)(4 * u1 * u2 / (1 -  (x - 1 * floor(x))^2 - y^2)^2);\n\nmatrix OP = op(Vh, Vh, solver = Crout, factorize = 1);\nmatrix B = bm(Vh, Vh, solver = CG, eps = 1e-20);\n\nreal[int] ev(nbrValeursPropres);\nVh[int] eV(nbrValeursPropres);\n\nint k = EigenValue(OP, B, sym = true, sigma = sigma, value = ev, vector = eV, tol=1e-10, maxit = 0, ncv = 0);\n\ncout << "Angle :' + f'{instance.angleTriangle}' + '" << endl;\n\nfor (int i = 0; i < k; i++)\n{\n    u1 = eV[i];\n\n    cout << "lambda[" << i << "] = " << ev[i] << ", normalisation : " << ev[i] * (pi - 3 * ' + f'{instance.angleTriangle}' + ') * 14 << ", err = " << int2d(Th)(dx(u1) * dx(u1) + dy(u1) * dy(u1) - 4 * (ev[i]) * u1 * u1 / (1 - (x - 1 * floor(x))^2 - y^2)^2) << endl;\n}'

    # Nom du fichier .edp
    nom = "trianglesHyperboliques"

    # Chemin du fichier .edp
    chemin = "path"

    # Chemin de l'exécutable modifié, si besoin, de FreeFem++
    cheminFF = "path"

    fichier = open(f"{chemin}{nom}.edp", "w")
    fichier.write(programme)
    fichier.close()

    texte = ""

    if cheminFF != "":
        texte = os.popen(f"{cheminFF} {chemin}{nom}.edp").read()
    else:
        texte = os.popen(f"{chemin}{nom}.edp").read()

    angles.append(i)
    valeursPropresNormalisees.append(float(texte.split("normalisation : ")[3].split(",")[0]))

plt.plot(angles, valeursPropresNormalisees)
plt.xlabel("Angles en radian")
plt.ylabel("1ère valeur propre normalisée")
plt.title("Valeur propre normalisée en fonction de l'angle")
plt.show()
###########################################################