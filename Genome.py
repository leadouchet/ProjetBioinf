#! /usr/bin/env python
# -*- coding:utf-8 -*-
import random
import string
import numpy as np
from Gene import Gene

class Genome() : 
	
	def __init__(self, nbr_gene, mean_size, genome = None) : 
		self.nbr_gene = nbr_gene
		self.mean_size = mean_size
		self.nbr_inversion = 0
		self.nbr_mutation = 0
		if genome == None : 
			self.genome = Create_genome()
		else : 
			self.genome = genome
		self.genome_ancestral = self.genome
			
	def Create_genome(self) : 
		positions = 0
		Genome = []
		for i in range(0,self.nbr_gene) : 
			taille_gene = int(random.gauss(self.mean_size,1))
			Genome += [Gene(taille_gene, position)]
			positions += taille_gene
		return(Genome)
	
	def Invertion_aleatoire(self) : 
		pos_1 = random.randint(0, self.nbr_gene-1)
		pos_2 = random.randint(0, self.nbr_gene-1)
		if pos_1 != pos_2 :
			self.nbr_inversion += 1 
			debut = min([pos_1, pos_2])
			fin = max([pos_1, pos_2])-1
			for g in range(debut, fin) : 
				self.Genome[g].Invertion_sequence()
		# D = DistanceGenome(Ancestrale, NewGenome)
		# Distance  = le nombre minimum de réarrangements qui transforme un génome en un autre,
		# Comparez le nombre d'inversions effectivement appliquées et la distance calculée
		# Nombre d'inversions (fonction du nombre de gènes) à partir duquel la distance estimée ne capture plus la distance évolutive. 
	
	def MutateGenome(self,m) : 
		for gene in self.genome : 
			self.nbr_mutation += gene.MutationPonctuelles(m)
	
	def Aligment(self) :
		nbr_mutation = 0 
		for gene in self.genome : 
			nbr_mutation += gene.NeedlemanWunsch()
		return(nbr_mutation)

"""	
	def NeedlemanWunsch(self,GenomeB) :
		# Creation de la matrice F 
		Sb = len(B)
		F = np.ones((Sa, Sb))
		d = -1 # Pénalité de trou
		S = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]) #Matrice de similarité
		for i in range(Sa) : 
			F[i, 0] = d*i
		for j in range(Sb) : 
			F[0,j] = d*j
		for i in range(1,Sa) : 
			for j in range(1,Sb) :
				pA = seq.index(A[i])
				pB = seq.index(B[j])
				Choice1 = F[i-1,j-1] + S[pA, pB]
				Choice2 = F[i-1, j] + d
				Choice3 = F[i, j-1] + d
				F[i, j] = max(Choice1, Choice2, Choice3)
	# Alignment
		AlignmentA = ""
		AlignmentB = ""
		i = Sa - 1
		j = Sb - 1
		Sc = Score = F[i, j]
		while (i > 0 and j > 0) : 
			Score = F[i, j]
			ScoreDiag = F[i - 1, j - 1]
			ScoreUp = F[i, j - 1]
			ScoreLeft = F[i - 1, j]
			pA = seq.index(A[i])
			pB = seq.index(B[j])
			if (Score == ScoreDiag + S[pA, pB]) : 
				AlignmentA = A[i] + AlignmentA
				AlignmentB = B[j] + AlignmentB
				i = i - 1
				j = j - 1
			elif (Score == ScoreLeft + d) : 
				AlignmentA = A[i] + AlignmentA
				AlignmentB = "-" + AlignmentB
				i = i - 1
			elif Score == ScoreUp + d : 
				AlignmentA = "-" + AlignmentA
				AlignmentB = B[j] + AlignmentB
				j = j - 1
		while (i >= 0) : 
			AlignmentA = A[i] + AlignmentA
			AlignmentB = "-" + AlignmentB
			i = i - 1
		while (j >= 0) : 
			AlignmentA = "-" + AlignmentA
			AlignmentB = B[j] + AlignmentB
			j = j - 1
		return((AlignmentA, AlignmentB, Sc))
"""
"""
19 inversions en 25 ans sur une bactérie comportant approximativement 4200 gènes, 
4Mb de séquence codante. 
environ 120 mutations ponctuelles sur la même durée
"""
