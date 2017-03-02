#!/usr/bin/python

from collections import defaultdict

_end = '_end'

class TrieSearch():
	def __init__(self, words):
		self.words = words

	def make_trie(self):
		"""Build trie tree index for str match."""
		self.root = dict()
		for word in self.words:
			current_dict = self.root
			for letter in word:
				current_dict = current_dict.setdefault(letter, {})
			current_dict[_end] = word
	
	def in_trie(self, word):
		"""Search 'word' against the trie tree."""
		if len(self.root) == 0: return []
		current_dict = self.root
		mapping = []
		for letter in word:
			if _end in current_dict:
				mapping.append(current_dict[_end])
			if letter in current_dict:
				current_dict = current_dict[letter]
			else:
				return mapping
		else:
			if _end in current_dict:
				mapping.append(current_dict[_end])
				return mapping
			else:
				return mapping
	
	def doSearch(self, data, first=['K','R']):
		"""Do batch search.
		Input: 
			1. seqs dict with fasta format
			2. only search peps started with first AA
		"""
		result = []
		peps = defaultdict(list)
		for head in data:
			seq = str.upper(data[head])
			firstindex = [i for i, ltr in enumerate(seq) if ltr in first]
			for i in firstindex:
				peps[seq[i:]].append(head)
		for p in peps:
			search = self.in_trie(p)
			if len(search) == 0: continue
			id = [str(i) for i in peps[p]]
			for s in search:
				result.append(''.join(first)+'\t'+ s +'\t'+ ';'.join(id))
		return result

def custom_formatwarning(msg, *a):
	return str(msg) + '\n'

if __name__ == '__main__':
	pass
