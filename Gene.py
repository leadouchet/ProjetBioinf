#! /usr/bin/env python
# -*- coding:utf-8 -*-
import random
import string
import numpy as np
from Bio import pairwise2

class Gene() : 
	
	def __init__(self, size, PositionDebut, posInit, sequence = None) : 
		self.posInit = posInit
		self.size = size
		self.debut = PositionDebut
		self.end = self.debut + self.size
		self.sens = 1
		self.bases = ["A", "T", "G" , "C"]
		self.mutation = 0
		self.inversions = 0
		if sequence == None : 
			self.seq = self.CreateGene(self.size)
		else : 
			self.seq = sequence
		self.seq_ancestral = str(self.seq)

	def CreateGene(self,taille) : 
		k = 0
		choix = []
		while k < taille : 
			k += 1
			choix += [random.choice(self.bases)]
		gene = string.join(choix,"")
		return gene
		
	def Invertion_sequence(self, m) :
		self.seq = self.seq[::-1] 
		self.sens *= -1 
		deb = self.end
		self.debut = self.end
		self.end = deb
		self.inversions +=1

	def MutationPonctuelles(self,m) :
		n = 0
		for base in range(self.size) : 
			mute = random.random()
			if mute < m : 
				seq = list(self.bases) # On ne peut pas muter avec la meme base
				n += 1
				seq.remove(self.seq[base])
				self.seq = self.seq[0:base] + random.choice(seq)+ self.seq[base+1 : self.size]
				self.mutation += 1
	
	def NeedlemanWunsch(self) :
		# Creation de la matrice F 
		F = np.ones((self.size, self.size))
		d = -1 # Pénalité de trou
		S = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]) #Matrice de similarité
		for i in range(self.size) : 
			F[i, 0] = d*i
		for j in range(self.size) : 
			F[0,j] = d*j
		for i in range(1,self.size) : 
			for j in range(1,self.size) :
				pA = self.bases.index(self.seq[i])
				pB = self.bases.index(self.seq_ancestral[j])
				Choice1 = F[i-1,j-1] + S[pA, pB]
				Choice2 = F[i-1, j] + d
				Choice3 = F[i, j-1] + d
				F[i, j] = max(Choice1, Choice2, Choice3)
	# Alignment
		#AlignmentA = ""
		#AlignmentB = ""
		i = self.size - 1
		j = self.size - 1
		Sc = Score = F[i, j]
		nbr_mutation = 0
		while (i > 0 and j > 0) : 
			Score = F[i, j]
			ScoreDiag = F[i - 1, j - 1]
			ScoreUp = F[i, j - 1]
			ScoreLeft = F[i - 1, j]
			pA = self.bases.index(self.seq[i])
			pB = self.bases.index(self.seq_ancestral[j])
			if (Score == ScoreDiag + S[pA, pB]) : 
				#AlignmentA = self.seq[i] + AlignmentA
				#AlignmentB = self.seq_ancestral[j] + AlignmentB
				i = i - 1
				j = j - 1
				if S[pA, pB] == 0 : 
					nbr_mutation += 1
			elif (Score == ScoreLeft + d) : 
				#AlignmentA = self.seq[i] + AlignmentA
				#AlignmentB = "-" + AlignmentB
				i = i - 1
				nbr_mutation += 1
			elif Score == ScoreUp + d : 
				#AlignmentA = "-" + AlignmentA
				#AlignmentB = self.seq_ancestral[j] + AlignmentB
				j = j - 1
				nbr_mutation += 1
		while (i >= 0) : 
			#AlignmentA = self.seq[i] + AlignmentA
			#AlignmentB = "-" + AlignmentB
			i = i - 1
		while (j >= 0) : 
			#AlignmentA = "-" + AlignmentA
			#AlignmentB = self.seq_ancestral[j] + AlignmentB
			j = j - 1
		return(nbr_mutation)

