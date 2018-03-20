#! /usr/bin/env python
# -*- coding:utf-8 -*-
from Gene import Gene
import random
from matplotlib import use
use('qt4agg')
import matplotlib.pyplot as plt

class Genome() : 
	
	def __init__(self, nbr_gene, mean_size, genome = None) : 
		self.nbr_gene = nbr_gene
		self.mean_size = mean_size
		self.nbr_inversion = 0
		self.nbr_mutation = 0
		self.generation = 0
		if genome == None : 
			self.genome = self.Create_genome()
		else : 
			self.genome = genome
		self.genome_ancestral = self.genome
			
	def Create_genome(self) : 
		positions = 0
		Genome = []
		for i in range(0,self.nbr_gene) : 
			taille_gene = int(random.gauss(self.mean_size,1))
			Genome += [Gene(taille_gene, positions)]
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
			gene.MutationPonctuelles(m)
			self.nbr_mutation += gene.mutation
	
	def Aligment(self) :
		nbr_mutation = 0 
		for gene in self.genome : 
			nbr_mutation += gene.NeedlemanWunsch()
		return(nbr_mutation)
	
	def Estimated_mutation(self) :
		ObservedMutation = 0 
		for gene in self.genome : 
			ObservedMutation += gene.NeedlemanWunsch()
		return(ObservedMutation)
			
	def Evolution_mutation(self, nbr_generation, m) :
		TrueMutation = [0]
		ObservedMutation = [0] 
		while self.generation < nbr_generation : 
			self.generation += 1
			self.MutateGenome(m)
			print(self.generation)
			#print(TrueMutation)
			#print(ObservedMutation)
			if self.generation % 100 == 0 : 
				print("yess")
				TrueMutation += [self.nbr_mutation]
				ObservedMutation += [self.Estimated_mutation()]
		return([TrueMutation, ObservedMutation])

		
"""
19 inversions en 25 ans sur une bactérie comportant approximativement 4200 gènes, 
4Mb de séquence codante. 
environ 120 mutations ponctuelles sur la même durée
"""
G1 = Genome(70, 3000)

res = G1.Evolution_mutation(1000,0.1)  #Drosophile 3,4*10-10
print(res)
plt.plot(res[0], label = "True")
plt.plot(res[1],label = "False")
plt.show()
