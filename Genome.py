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
			Genome += [Gene(taille_gene, positions, i)]
			positions += taille_gene
		return(Genome)
	
	def Inversion_aleatoire(self) : 
		pos_1 = random.randint(0, self.nbr_gene-1)
		pos_2 = random.randint(0, self.nbr_gene-1)
		if pos_1 != pos_2 :
			self.nbr_inversion += 1 
			debut = min([pos_1, pos_2])
			fin = max([pos_1, pos_2])-1
			invertseq = self.genome[debut:fin]
			self.genome = self.genome[0:debut] + invertseq[::-1] + self.genome[fin : self.nbr_gene]  
			for g in range(debut, fin) : 
				self.genome[g].Inversion_sequence()

	def MutateGenome(self,m) : 
		for gene in self.genome : 
			gene.MutationPonctuelles(m)
			self.nbr_mutation += gene.mutation

	
	def Estimated_mutation(self) :
		ObservedMutation = 0 
		for gene in self.genome : 
			#ObservedMutation += gene.NeedlemanWunsch()
			ObservedMutation += gene.nbMut()
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
		
	def Evolution_inversion(self, nbr_inversion) : 
		inversion_observee = [0]
		while self.nbr_inversion < nbr_inversion : 
			self.Inversion_aleatoire()
			if self.nbr_inversion % 100 == 0 : 
				graph = self.Comp_graph()
				print(graph)
				nbr_cycle = self.number_of_cycle(graph)
				print(nbr_cycle)
				inversion_observee += [self.nbr_gene - nbr_cycle]
		return(inversion_observee)

	def Comp_graph(self):
		comp_graphe = dict()
		for gene in self.genome_ancestral :
			g = gene.posInit
			gd = str(g)+"_d"
			gf = str(g)+"_f"
			if (g == self.genome_ancestral[0].posInit):
				comp_graphe["D"] = [gd]
				comp_graphe[gd] = ["D"]
				comp_graphe[gf] = [str(g+1)+"_d"]
			elif (g == self.genome_ancestral[-1].posInit):
				comp_graphe[gd] = [str(g-1)+"_f"]
				comp_graphe[gf] = ["F"]
				comp_graphe["F"] = [gf]
			else:
				comp_graphe[gd] = [str(g-1)+"_f"]
				comp_graphe[gf] = [str(g+1)+"_d"]
		
		for i in range(len(self.genome)):
			g = self.genome[i].posInit
			gav = self.genome[i-1].posInit
			gap = self.genome[(i+1)%len(self.genome)].posInit
			if self.genome[i].sens == 1:
				gd = str(g)+"_d"
				gf = str(g)+"_f"
			else:
				gd = str(g)+"_f"
				gf = str(g)+"_d"

			if self.genome[i-1].sens==1:
				gf_av = str(gav)+"_f"
			else:
				gf_av = str(gav)+"_d"

			if self.genome[(i+1)%len(self.genome)].sens ==1:
				gd_ap = str(gap)+"_d"
			else:
				gd_ap = str(gap)+"_f"

			if (i==0):
				comp_graphe["D"].append(gd)
				comp_graphe[gd].append("D")
				comp_graphe[gf].append(gd_ap)
			elif (i==len(self.genome)-1):
				comp_graphe[gd].append(gf_av)
				comp_graphe[gf].append("F")
				comp_graphe["F"].append(gf)
			else:
				comp_graphe[gd].append(gf_av)
				comp_graphe[gf].append(gd_ap)
		return(comp_graphe)

	def number_of_cycle(self,Dico):
		c = 0
		c2 = 0
		queue = list(Dico.keys())
		while len(queue) != 0:
			myKey = queue[0]
			queue.remove(myKey)
			visited = [myKey]
			start = Dico[myKey][0]
			end = Dico[myKey][1]
			visited.append(start)
			visited.append(end)
			queue.remove(start)
			if (start != end):
				queue.remove(end)
			else:
				c2 += 1
			while 1:
				if (start == end):
					break
	#			print("myKey : ",myKey, "start : ", start, "end : ",end)
	#			print("---------------------------------")
				voisins = Dico[start]
				if voisins[0] == myKey:
					myKey = start
					start = voisins[1]
				else:
					myKey = start
					start = voisins[0]
				visited.append(start)
				if (start != end):
					queue.remove(start)
			c += 1
		return (c)
		
"""
19 inversions en 25 ans sur une bactérie comportant approximativement 4200 gènes, 
4Mb de séquence codante. 
environ 120 mutations ponctuelles sur la même durée
"""
G1 = Genome(10, 3000)
I = G1.Evolution_inversion(1)
plt.plot(range(0,len(I)), I , label = "Inversion observées")
plt.plot(range(0,len(I)), range(0,len(I)), label = "Inversion effectuée")
plt.show()
plt.legend()


#G1.Evolution_mutation(100,0.000001)  #Drosophile 3,4*10-10
