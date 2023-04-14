#!/usr/bin/env python
# coding: utf-8

# # Université de Mons
# ## Projet d'Analyse Numérique
# ## Thème: Sociobiologie
# ### Auteur: Sinclair TSANA


import numpy as np
from scipy.integrate import odeint
from matplotlib import pyplot as plt


def gain_moyen_individu(A, i, x):
    """Gain moyen pour de toutes les rencontres d'un individu du groupe i"""
    return sum([A[i][j] * x[j] for j in range(len(x))])

def gain_moyen_population(x):
    """Gain moyen pour toute la poppulation"""
    return sum([x[j] * gain_moyen_individu(A, j, x) for j in range(len(x))])

def dxdt(A, x, t):
    """Système dynamique"""
    return [x[i] * (gain_moyen_individu(A, i, x) - gain_moyen_population(x))
            for i in range(len(x))]


# Question 8
# on peut supposer que la matrice est sous la forme indiquée à la question 7

class Point:
    def __init__(self, x):
            self.x1 = x[0]
            self.x2 = x[1]
            self.x3 = x[2]
            
class Segment:
    def __init__(self, x, xtild):
        self.x = x
        self.xtild = xtild

# Routine renvoyant la liste des equilibres de S^N
e1 = [1,0,0]
e2 = [0,1,0],
e3 = [0,0,1]
e4 = [1/3,1/3,1/3]
def equilibres(A):
    """prend en entrée une matrice 3x3 et renvoie la liste des equilibres de appartenant au simplexe S^N"""
    if not np.array_equal(A, np.zeros([3, 3], dtype = float)):
        # On met la matrice sous la forme indiquée à la question 7
        # a chaque colonne j, on retranche l'element  a_1j de la premiere ligne
        for j in range(len(A)):
            A[:,j]=A[:,j]-A[0][j]
        
        #on supprime la première ligne de la nouvelle matrice
        A = np.delete(A,0,0).reshape(2, 3)
        
        #on calcule le determinant des sous-matrices
        d1= np.linalg.det(np.delete(A, 0, 1))
        d2= np.linalg.det(np.delete(A, 1, 1))
        d3= np.linalg.det(np.delete(A, 2, 1))
        equi = [e1,e2,e3,e4]  #Liste qui va contenir tous les equilibres
        res = [e1,e2,e3,e4]   #Liste contenant les equilibres de S^N
        if (d1 !=0 and d2!=0 and d3 !=0):
            x = [d1,d2,d3]
            equi.append(x)
            if sum(x)==1:
                p = Point(x)
                res.append(p)
        elif d1 == 0:
            x = [0, d2, d3]
            equi.append(x)
            if sum(x)==1:
                p = Segment(d2,d3)
                res.append()
        elif d2 == 0:
            x = [d1,0,d3]
            equi.append(x)
            if sum(x)==1:
                p = Segment(d1,d3)
                res.append(p)
        elif d3 == 0:
            x = [d1,d2,0]
            equi.append(x)
            if sum(x)==1:
                p = Segment(d1,d2)
                res.append(p)
        else:
            print("On a une infinité de solution")
    else:
        exit('Cas dégénéré: la matrice est nulle')
        
    return res

#Question 9
def Graphe(A):
    def f(x,t):
        return [x[i] * (gain_moyen_individu(A, i, x) - gain_moyen_population(x)) for i in range(len(x))]
    
    t = np.arange(0, 100, .1)
    
    for s0 in np.arange(0.05, 0.95, .1):
        #1-s0 to ensure that population vector lies in simplex
        for r0 in np.arange(0.05, 1-s0, .1):
            y0 = [s0, r0]
            y = odeint(f, y0, t)
            X, Y = y[:,1], y[:,0]
            plt.plot(X, Y)
            
            # Ajout des directions
            adx0, adx1 = 0, len(t) // 3
            adx2 = adx1 * 2
            arrow0 = X[adx0+1], Y[adx0+1], X[adx0+1]-X[adx0], Y[adx0+1]-Y[adx0]
            arrow1 = X[adx1+1], Y[adx1+1], X[adx1+1]-X[adx1], Y[adx1+1]-Y[adx1]
            arrow2 = X[adx2+1], Y[adx2+1], X[adx2+1]-X[adx2], Y[adx2+1]-Y[adx2]
            plt.arrow(*arrow0, shape='full', lw=0, length_includes_head=True, head_width=0.02) 
            plt.arrow(*arrow1, shape='full', lw=0, length_includes_head=True, head_width=0.02)  
            plt.arrow(*arrow2, shape='full', lw=0, length_includes_head=True, head_width=0.02)
    plt.show()

# Question 10
# pour V=2 et C= 4 on obtient la matrice suivante 
M = np.array([[-1,2,-1],[0,1,1],[-1,1,1]])
Graphe(M)


#Question 14
A1 = np.array([[0,0,1],[0,0,1],[1,1,0]])
A2 = np.array([[0,0,1],[0,0,1],[1,-1,0]])
A3 = np.array([[0,1,2],[1,0,1],[2,-1,0]])
A4 = np.array([[0,1,1],[1,0,1],[1,1,0]])
A5 = np.array([[0,-1,3],[-1,0,3],[1,1,0]])
A6 = np.array([[0,1,1],[-1,0,3],[1,1,0]])
A7 = np.array([[0,1,-1],[1,0,-3],[-1,3,0]])
A8 = np.array([[0,0,1],[0,0,0],[1,-1,0]])
A9 = np.array([[0,0,1],[0,0,0],[1,1,0]])
A10 = np.array([[0,0,-1],[0,0,1],[-1,0,0]])
A11 = np.array([[0,3,1],[3,0,1],[1,1,0]])
A12 = np.array([[0,1,1],[-1,0,3],[1,3,0]])

liste_matrices = [A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12]
L = ['A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12']

    
for i in range(len(L)):
    print('Trajectoires pour la matrice {}'.format(L[i]))
    Graphe(liste_matrices[i])
