#!/usr/bin/env python

__author__ = "Gaja Dreszler, Bioinformatyka"

"""Backtracking version of script coRNA."""

import networkx as nx
import warnings
import argparse
import matplotlib.pyplot as plt
import random
import time
import re
import itertools

class BACKTRACKING():

    def __init__(self):

        self.filename = ""
        self.structures = []
        self.constraints = ""
        self.constraints_nodes = ""
        self.G = nx.Graph()
        self.bases = ['A', 'C', 'G', 'U']
        self.all_possible_characters = ['A', 'C', 'G', 'U', 'W', 'S', 'M', 'K', 'R', 'Y', 'B', 'D', 'H', 'V', 'N', '-']
        self.basepairs = {
            'G': ['C', 'U'],
            'C': ['G'],
            'A': ['U'],
            'U': ['A', 'G'],
            'W': ['G', 'U'],
            'S': ['G', 'C', 'U'],
            'M': ['U', 'G'],
            'K': ['A', 'C', 'G', 'U'],
            'R': ['U', 'C'],
            'Y': ['G', 'A'],
            'B': ['A', 'C', 'G', 'U'],
            'D': ['A', 'C', 'G', 'U'],
            'H': ['A', 'C', 'G', 'U'],
            'V': ['A', 'C', 'G', 'U'],
            'N': ['A', 'C', 'G', 'U']
            }
        self.constraints_dict = {
                    'N': ['A', 'C', 'G', 'U'],
                    'W': ['A', 'U'],
                    'S': ['C', 'G'],
                    'M': ['A', 'C'],
                    'K': ['G', 'U'],
                    'R': ['A', 'G'],
                    'Y': ['C', 'U'],
                    'B': ['C', 'G', 'U'],
                    'D': ['A', 'G', 'U'],
                    'H': ['A', 'C', 'U'],
                    'V': ['A', 'C', 'G']
                    }
        self.characters = [".", "(", ")", "[", "]"]
        self.empty_nodes = []
        self.seq_options = []
        self.universal_sequence = []
        self.random_seq = ""
        self.mutation_seq = ""
        self.structures_copy = []
        self.list_of_bps = []
        self.set_filename = ""
        self.most_stable_motifs = []
        self.motifs_found = []
        self.prevented_patterns = []

    def file_parser(self):
        """Parses filename and optional information from a user."""
      
        parser = argparse.ArgumentParser()
        parser.add_argument("-f", "--filename", required = True, dest = "name", help = "Name of a file that contains secondary structures and constraints.")
        args = parser.parse_args()
        self.filename = args.name

    def get_data_from_file(self, filename):
        """Gets data from a given file."""

        file1 = open(filename, "r")
        try: 
            data = file1.readlines()
        except IOError:
            print "Error occured during reading a file.\n"
        finally:
            file1.close()
        return data
    
    def get_structures_from_data(self):
        """Saves file information in proper variables."""
        
        data = self.get_data_from_file(self.filename)
        for i in data:
            if i[0] in self.characters:
                self.structures.append(i.strip())
            elif i.startswith(">prevented"):
                self.prevented_patterns.append(i[10:].strip())
            elif (i[0].isalpha() or i.startswith("-"))and self.constraints == "":
                self.constraints = i.strip()
                
        self.structures_copy = self.structures[:]
        if self.constraints == "":
            self.constraints = "N" * len(self.structures[0])

    def check_input(self):
        """Checks if input file is correct."""

        for i in range(1, len(self.structures)):
            if len(self.structures[i]) != len(self.structures[i-1]):
                raise Exception("Various lengths of structures. Check input file.")

        if self.constraints != "" and len(self.constraints) != len(self.structures[0]):
            raise Exception("Various lengths of structures and constraints. Check input file.")
        
        for j in self.constraints:
            if j not in self.all_possible_characters:
                raise Exception("Wrong character ", j, " in constraints. Check input file.")
        
        for k in self.structures:
            for l in k:
                if l not in self.characters:
                    raise Exception("Wrong character ", l, " in structure. Check input file.")

        if "-" in self.constraints:
            self.constraints = list(self.constraints)
            for l in range(len(self.constraints)):
                if self.constraints[l] == "-":
                    self.constraints[l] = "N"
            "".join(self.constraints)
    
    def check_brackets(self, strs):
        """Checks if all secondary structures are correct."""
        
        open_brackets = []
        open_brackets_pseudoknots = []
        for structure in strs:
            for nt in range(len(structure)):
                if structure[nt] == "(":
                    open_brackets.append(nt)
                elif structure[nt] == "[":
                    open_brackets_pseudoknots.append(nt)
                elif structure[nt] == ")":
                    if all(x in strs for x in self.structures):
                        self.G.add_edge(open_brackets[-1]+1, nt+1)
                        self.list_of_bps.append((open_brackets[-1], nt))
                    try:
                        open_brackets.pop()
                    except IndexError:
                        pass
                elif structure[nt] == "]":
                    if all(x in strs for x in self.structures):
                        self.G.add_edge(open_brackets_pseudoknots[-1]+1, nt+1)
                        self.list_of_bps.append((open_brackets[-1], nt))
                    try:
                        open_brackets_pseudoknots.pop()
                    except IndexError:
                        pass
                elif structure[nt] != ".":
                    raise Exception("Wrong character in structures. Check input file.")
        
        if open_brackets == [] and open_brackets_pseudoknots == []:
            return True
        else:
            return False

    def draw_graph(self):
        """Draws a single graph for all structures."""

        for i in range(len(self.structures[0])):
            self.G.add_node(i+1)
            
        result = self.check_brackets(self.structures)
        if result is False:
            raise Exception("Wrong base pairings in structures. Check input file.")

        degrees = self.G.degree()
        degreeslist = degrees.values()
        for i in range(len(degreeslist)):
            if degreeslist[i] == 0:
                self.G.remove_node(i+1)
                # print nx.nodes(self.G)

        pos=nx.spring_layout(self.G)#draw graph
        nx.draw_networkx_nodes(self.G,pos,
                       node_color='r',
                       node_size=500,
                   alpha=0.8)
        nx.draw_networkx_edges(self.G,pos,width=1.0,alpha=0.5)
        labels={}
        for i in nx.nodes(self.G):
            labels[i]=r'$'+str(i)+'$'
        nx.draw_networkx_labels(self.G,pos,labels,font_size=16)
        plt.axis('off')
        # plt.show()

    def check_if_bipartite(self):
        """Checks if graph is bipartite."""

        if nx.is_bipartite(self.G):
            pass
        else:
            raise_exc = 1
            for i in range(len(self.structures_copy)):
                self.structures = self.structures_copy[:]
                self.structures.remove(self.structures_copy[i])
                self.G.clear()
                self.draw_graph()
                if nx.is_bipartite(self.G):
                    print "Structure number ", i+1 , " is inconsistent."
                    raise_exc = 0
                    break
            if raise_exc == 1:
                raise Exception("Inconsistent structures. Unsolvable case.")

    def get_universal_sequence(self):
        """Creates universal sequence for given structures."""

        self.constraints_nodes = ""
        self.random_seq = ""
        self.universal_sequence = []
        for i in self.G.nodes_iter():
            if self.constraints != "":
                self.constraints_nodes += self.constraints[i-1] # only nodes constraints
                self.random_seq += self.constraints[i-1]
            else:
                self.constraints_nodes += "N"
                self.random_seq += "N"

        for i in range(len(nx.nodes(self.G))):
            self.G.node[nx.nodes(self.G)[i]]['nt'] = self.constraints_nodes[i]
            self.universal_sequence.append(self.constraints_nodes[i])
            if self.constraints_nodes[i] == "N":
                self.empty_nodes.append(nx.nodes(self.G)[i])

        # for i in self.G.nodes_iter(): #print values stored in nodes
            # print self.G.node[i]['nt']

        stored_values_before = []
        stored_values_after = []
        while(True):
            for i in range(len(self.universal_sequence)):
                if self.constraints_nodes[i] == "N":
                    nt = random.choice(self.bases)
                    result = self.assignable_universal_nt(nx.nodes(self.G)[i], i, nt)
                    if result == False:
                        raise Exception("Constraints are inconsistent. Sequence assignment is impossible.")
                    self.bases = ['A', 'C', 'G', 'U']

            for i in self.universal_sequence: 
                stored_values_after.append(i)

            if stored_values_after == stored_values_before: # no changes -> break
                break
            else:
                stored_values_before = stored_values_after[:]
                stored_values_after = []

        for node in range(len(self.universal_sequence)):
            for key, value in self.constraints_dict.items():
                contained = 0
                if len(self.universal_sequence[node]) > 1:
                    for letter in self.universal_sequence[node]:
                        if letter in value:
                            contained += 1
                    if contained == len(value) and len(self.universal_sequence[node]) == contained:
                        self.universal_sequence[node] = key
                        break
        return self.universal_sequence

    def assignable_universal_nt(self, node, place_in_seq, nt):
        """Tries to assign a random nucleotide to given position in sequence."""

        neighbors = self.G.neighbors(node)
        right_neighbors = []
        for j in neighbors:
            for key, value in self.basepairs.items():
                if self.G.node[j]['nt'] == "N" or (nt == key and self.G.node[j]['nt'] in value) or (self.G.node[j]['nt'] == key and nt in value):
                    if self.universal_sequence[place_in_seq].startswith("N"):
                        self.universal_sequence[place_in_seq] = self.universal_sequence[place_in_seq][1:]
                    if nt not in self.universal_sequence[place_in_seq]:
                        right_neighbors.append(j)
                        break

        if len(right_neighbors) == len(neighbors):
            self.universal_sequence[place_in_seq] += nt
            if "N" in self.random_seq:
                self.random_seq = list(self.random_seq)
                self.random_seq[place_in_seq] = nt
                self.random_seq = "".join(self.random_seq)

        self.bases.remove(nt) # removes used nt from list self.bases
        if self.bases == []:
            return False
        nt = random.choice(self.bases)
        result = self.assignable_universal_nt(node, place_in_seq, nt)
    
    def check_forbidden_nts(self, seq):
        """Checks if forbidden patterns are present in a given sequence."""

        all_patterns = []
        for pt in self.prevented_patterns:
            temp = []
            in_bases = 0
            for i in pt:
                if i not in self.bases:
                    [temp.append(x) for x in pt]
                    temp = self.convert_shortcuts(temp)
                    for j in range(len(temp)):
                        if len(temp[j]) > 1:
                            temp[j] = list(temp[j])
                            
                    all_pt = []
                    all_pt = list(itertools.product(*temp)) # get all combinations of a pattern
                    all_pt = ["".join(x) for x in all_pt]
                    all_patterns.extend(all_pt)
                    break
                else:
                    in_bases += 1
            if in_bases == len(pt):
                all_patterns.append(pt)
        
        sequence = []
        [sequence.append(x) for x in seq]
        sequence = self.convert_shortcuts(seq)
        
        for i in range(len(sequence)):
            if len(sequence[i]) > 1:
                sequence[i] = list(sequence[i])
        
        all_sequences = []
        all_sequences = list(itertools.product(*sequence)) # get all combinations of a sequence
        all_sequences = ["".join(x) for x in all_sequences]
        
        forbidden_seqs = []
        for pattern in all_patterns:
            for s in all_sequences:
                if pattern in s and seq not in forbidden_seqs:
                    forbidden_seqs.append(s)
        return forbidden_seqs
    
    def convert_shortcuts(self, sequence):
        """Converts shortcuts into all possible nucleotides in given positions."""

        # change more complicated letters into all possible letters in this position
        for i in range(len(sequence)):
            if sequence[i][0] not in self.bases:
                sequence[i] = "".join(self.constraints_dict.get(sequence[i]))
        return sequence
    
    def get_random_mutation(self, seq_from_program):
        """Incorporates random mutations into a sequence."""
        
        constr = []
        [constr.append(x) for x in self.constraints]
        self.seq_options = self.convert_shortcuts(constr)
        used_bases = 0
        timeout = time.time() + 10
        mutations = 0
        while(True):
            if time.time() > timeout:
                print "Adding random mutations failed."
                break

            self.mutation_seq = seq_from_program
            while(True):
                random_position = random.randint(1, len(seq_from_program))
                if len(self.seq_options[random_position-1]) > 1:
                    break
            bases = self.seq_options[random_position-1]
            random_mutation = random.choice(bases)
            
            bps = []
            if self.mutation_seq[random_position-1] != random_mutation:
                [bps.append(item) for item in self.list_of_bps if random_position-1 in item]
                if bps != []:
                    for pair in bps:
                        for number in pair:
                            if number != random_position-1:
                                options = []
                                options = self.basepairs.get(random_mutation)[:]
                                while(True):
                                    if options == []:
                                        break
                                    bps_mutation = random.choice(options)
                                    if (bps_mutation in self.seq_options[number]) and (bps_mutation in self.mutation_seq[number] or len(self.seq_options[number]) > 1):
                                        if bps_mutation not in self.mutation_seq[number]:
                                            self.mutation_seq = list(self.mutation_seq)
                                            self.mutation_seq[number] = bps_mutation
                                            self.mutation_seq = "".join(self.mutation_seq)

                                        self.mutation_seq = list(self.mutation_seq)
                                        self.mutation_seq[random_position-1] = random_mutation
                                        self.mutation_seq = "".join(self.mutation_seq)
                                        mutations += 1
                                        break         
                                    else:
                                        options.remove(bps_mutation)
                else:
                    self.mutation_seq = list(self.mutation_seq)
                    self.mutation_seq[random_position-1] = random_mutation
                    self.mutation_seq = "".join(self.mutation_seq)
                    mutations += 1
            else:
                continue
        return self.mutation_seq

    def calling_function(self):
        """Calls all functions in the right order."""

        self.file_parser()
        self.get_structures_from_data()
        self.check_input()
        self.draw_graph()
        self.check_if_bipartite()
        result = self.get_universal_sequence()
        self.random_seq = list(self.random_seq)

        if self.prevented_patterns != []:
            forbidden = self.check_forbidden_nts(self.random_seq)
            if forbidden != []:
                print "Forbidden pattern found."
            else:
                print "No forbidden patterns found."

def main():

    inverse_folding = BACKTRACKING()
    inverse_folding.calling_function()


if __name__ == "__main__":  
    main()
      