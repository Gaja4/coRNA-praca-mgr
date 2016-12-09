#!/usr/bin/env python

"""
Class coRNA serves as an initial checking of data given by user to DesiRNA design program.
It's task is to detect and raise proper errors if a given mistake occurs.
It also checks if initial secondary structures and constraints are suitable and if any
sequence matching those structures and constraints, can be created.
Besides that, this class incorporates random mutations to sequence and returns a new alternative sequence.

Input: file in fasta format with comments begining with ">" character, containing secondary
structures in dot-bracket notation and constraints with all possible shortcuts
optional: - lines starting with ">prevented " and a sequence that can not show up in a solution sequence
optional: - file with a library of the most stable motifs, that can be present in initial secondary structures

Usage: "python [name of a program].py -f [main filename].fasta -s [library filename].fasta -v [0 or 1]
-m [number of mutation cycles] -p [probability]"
"""

__author__ = "Gaja Klaudel, Bioinformatyka"

import argparse
import re
import itertools
import pycosat
import random
from numpy.random import choice


class coRNA:

    def __init__(self):

        self.filename = ""
        self.structures = []
        self.constraints = ""
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
        self.dict_logic_nt_values = {
            'A': (0, 0),
            'U': (0, 1),
            'G': (1, 0),
            'C': (1, 1)
        }
        self.dict_additional_expressions = {
                    'W': lambda xi, yi: [[-xi]],
                    'S': lambda xi, yi: [[xi]],
                    'Y': lambda xi, yi: [[yi]],
                    'R': lambda xi, yi: [[-yi]],
                    'M': lambda xi, yi: [[-xi, yi], [xi, -yi]],
                    'K': lambda xi, yi: [[-xi, -yi], [xi, yi]],
                    'B': lambda xi, yi: [[xi, yi]],
                    'D': lambda xi, yi: [[-xi, -yi]],
                    'H': lambda xi, yi: [[-xi, yi]],
                    'V': lambda xi, yi: [[xi, -yi]],
                    'A': lambda xi, yi: [[-xi], [-yi]],
                    'U': lambda xi, yi: [[-xi], [yi]],
                    'G': lambda xi, yi: [[xi], [-yi]],
                    'C': lambda xi, yi: [[xi], [yi]]
        }
        self.seq_options = []
        self.random_seq = ""
        self.mutation_seq = ""
        self.list_of_bps = []
        self.set_filename = ""
        self.most_stable_motifs = []
        self.motifs_found = []
        self.prevented_patterns = []
        self.verbose = ""
        self.cnf = []
        self.converted_sequence = []
        self.mutation_cycles = 0
        self.solutions = []
        self.mutated_seqs_used = []
        self.cnf_mutations = []
        self.mutation_history = []
        self.cnf_forbidden_patterns = []
        self.motifs_locations = []
        self.sequences_to_check = []
        self.cnf_motifs = []
        self.sequences_stable_motifs = []
        self.probability = 0

    def file_parser(self):
        """Parses filename and optional information from a user."""

        parser = argparse.ArgumentParser(prog='inverse_folding_cnf', usage='%(prog)s [options]')
        parser.add_argument("-f", "--filename", required=True, dest="name",
                            help="Name of a file that contains secondary structures and constraints.")
        parser.add_argument("-s", "--set", required=False, dest="set",
                            help="Name of a file containing set of the most stable motifs.")
        parser.add_argument("-m", "--mutations", required=False, default=100, dest="number_of_mutations", type=int,
                            help="Number of mutation cycles.")
        parser.add_argument("-v", "--verbose", required=False, dest="verbose", default=0, type=int, choices=[0, 1],
                            help="Displaying solution from PYCOSAT solver.")
        parser.add_argument("-p", "--probability", required=False, dest="prob", default=0.5, type=float,
                            help="Probability of mutating the whole fragment from library of motifs.")
        args = parser.parse_args()
        self.filename = args.name
        self.set_filename = args.set
        self.verbose = args.verbose
        self.probability = args.prob
        self.mutation_cycles = args.number_of_mutations

    def get_data_from_file(self, filename):
        """Gets data from a given file. Takes a filename and returns data from it."""

        data = []
        file1 = open(filename, "r")
        try:
            data = file1.readlines()
        except IOError:
            print "Error occured during reading a file.\n"
        finally:
            file1.close()
        return data

    def get_structures_from_data(self, data):
        """Gets information about structures, prevented patterns and constraints from file data
         and saves them in proper variables."""

        for i in data:
            if i[0] in self.characters:
                self.structures.append(i.strip())
            elif i.startswith(">prevented"):
                self.prevented_patterns.append(i[10:].strip())
            elif (i[0].isalpha() or i.startswith("-")) and self.constraints == "":
                self.constraints = i.strip()

        if self.constraints == "":
            self.constraints = "N" * len(self.structures[0])

    def check_input(self):
        """Checks if input file is correct(the same lengths of all structures, structures and constraints
         and if constraints and structures consist only of proper characters).
         If any error detected, it will raise a suitable exception message."""

        for i in range(1, len(self.structures)):
            if len(self.structures[i]) != len(self.structures[i-1]):
                raise ValueError("Various lengths of structures. Check input file.")

        if self.constraints != "" and len(self.constraints) != len(self.structures[0]):
            raise ValueError("Various lengths of structures and constraints. Check input file.")

        for j in self.constraints:
            if j not in self.all_possible_characters:
                raise ValueError("Wrong character ", j, " in constraints. Check input file.")

        for k in self.structures:
            for l in k:
                if l not in self.characters:
                    raise ValueError("Wrong character ", l, " in structure."
                                                            " Check input file.")

        if "-" in self.constraints:
            self.constraints = list(self.constraints)
            for l in range(len(self.constraints)):
                if self.constraints[l] == "-":
                    self.constraints[l] = "N"
            "".join(self.constraints)

        if self.prevented_patterns:
            for l in self.prevented_patterns:
                if l.startswith("N") or l.endswith("N") or all(x == "N" for x in l):
                    raise ValueError("Wrong forbidden patterns. Check input file.")

    def check_brackets(self, strs):
        """Checks if initial secondary structures are created due to dot-bracket notation rules
        and creates a list with numbers of nucleotides that basepair with each other in these structures.
        Returns True if structures are correct, and False - otherwise."""

        open_brackets = []
        open_brackets_pseudoknots = []

        for structure in strs:
            for nt in range(len(structure)):
                if structure[nt] == "(":
                    open_brackets.append(nt)
                elif structure[nt] == "[":
                    open_brackets_pseudoknots.append(nt)
                elif structure[nt] == ")":
                    if all(x for x in strs):
                        if strs == self.structures:
                            self.list_of_bps.append((open_brackets[-1], nt))
                    try:
                        open_brackets.pop()
                    except IndexError:
                        pass
                elif structure[nt] == "]":
                    if all(x for x in strs):
                        if strs == self.structures:
                            self.list_of_bps.append((open_brackets[-1], nt))
                    try:
                        open_brackets_pseudoknots.pop()
                    except IndexError:
                        pass
                elif structure[nt] != ".":
                    raise ValueError("Wrong character in structures. Check input file.")
        if open_brackets == [] and open_brackets_pseudoknots == []:
            return True
        else:
            return False

    def check_constraints_logic(self):
        """Checks if given constraints are suitable for initial structures
        and if it is possible to create a sequence for given constraints."""

        constr = []
        [constr.append(x) for x in self.constraints]

        for pair in self.list_of_bps:
            i = constr[pair[0]]
            j = constr[pair[1]]
            xi = 2 * pair[0] + 1
            yi = 2 * pair[0] + 2
            xj = 2 * pair[1] + 1
            yj = 2 * pair[1] + 2
            # append main formula to cnf
            self.cnf.append([-xi, xj, -yi])
            self.cnf.append([xi, -xj, -yj])
            self.cnf.append([-yi, -yj])
            self.cnf.append([yi, yj])
            # add i formulas
            if i in self.dict_additional_expressions.keys() and i not in self.bases:
                for x in self.dict_additional_expressions[i](xi, yi):
                    self.cnf.append(x)
            elif i in self.dict_logic_nt_values.keys():
                if self.dict_logic_nt_values[i][0] == 0:
                    xi = -xi
                if self.dict_logic_nt_values[i][1] == 0:
                    yi = -yi
                if [xi, yi] not in self.cnf:
                    self.cnf.append([xi, yi])
            # add j formulas
            if j in self.dict_additional_expressions.keys() and j not in self.bases:
                for x in self.dict_additional_expressions[j](xj, yj):
                    self.cnf.append(x)
            elif j in self.dict_logic_nt_values.keys():
                if self.dict_logic_nt_values[j][0] == 0:
                    xj = -xj
                if self.dict_logic_nt_values[j][1] == 0:
                    yj = -yj
                if [xj, yj] not in self.cnf:
                    self.cnf.append([xj, yj])

        formula = self.create_formula_whole_sequence(self.constraints, 0)
        self.cnf += formula
        self.cnf += self.cnf_forbidden_patterns
        solution = pycosat.solve(self.cnf)
        return solution

    def check_forbidden_nts(self):
        """Creates a list of all possible forbidden patterns and converts it into a logic formula
         that forbidds their presence in every position in the sequence."""

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
                    all_pt = list(itertools.product(*temp))  # get all combinations of a forbidden pattern
                    all_pt = ["".join(x) for x in all_pt]
                    all_patterns.extend(all_pt)
                    break
                else:
                    in_bases += 1
            if in_bases == len(pt):
                all_patterns.append(pt)
        all_patterns = list(set(all_patterns))  # unique list of patterns

        # create phrases of forbidden patterns
        for pattern in all_patterns:
            x = 0
            for i in range(len(self.constraints) - (len(pattern) - 1)):
                k = 0
                phrase = []
                for j in range(x, len(pattern) + x):
                    xi = 2 * j + 1
                    yi = 2 * j + 2
                    if self.dict_logic_nt_values[pattern[k]][0] == 1:
                        xi = -xi
                    if self.dict_logic_nt_values[pattern[k]][1] == 1:
                        yi = -yi
                    phrase.append(xi)
                    phrase.append(yi)
                    k += 1
                x += 1
                self.cnf_forbidden_patterns.append(phrase)

        # unique list of forbidden patterns
        self.cnf_forbidden_patterns = [list(x) for x in set(tuple(x) for x in self.cnf_forbidden_patterns)]

    def convert_shortcuts(self, sequence):
        """Converts nucleotides shortcuts into a list of all possible nucleotides in a given position.
        Takes a list with latters from a sequence. Returns list with nucleotides without shortcuts."""

        # change more complicated letters into all possible nts in a given position
        for i in range(len(sequence)):
            if sequence[i][0] not in self.bases:
                sequence[i] = "".join(self.constraints_dict.get(sequence[i]))
        return sequence

    def create_formula_whole_sequence(self, sequence, con):
        """Returns a logic formula that describes a given sequence.
        Con = 1 means that a whole sequence should be treated like constraints."""

        phrase = [[]]
        for i in range(len(sequence)):
            xi = 2 * i + 1
            yi = 2 * i + 2
            # not in constraints
            if (self.constraints[i] == 'N' or self.constraints == "" or self.constraints[i] == '-') \
                    and sequence[i] != 'N' and con != 1 and sequence[i] != '-':
                if self.dict_logic_nt_values[sequence[i]][0] == 0:
                    xi = -xi
                if self.dict_logic_nt_values[sequence[i]][1] == 0:
                    yi = -yi
                phrase[0].append(xi)
                phrase[0].append(yi)
            # in constraints
            elif (sequence[i] == self.constraints[i] and sequence[i] != 'N' and sequence[i] != '-') or con == 1:
                for x in self.dict_additional_expressions[sequence[i]](xi, yi):
                    phrase.append(x)

        if not phrase[0]:
            phrase.pop(0)
        return phrase

    def unfreeze_basepairs(self, position, seq, try_cnf):
        """Unlocks given positions in sequence, so that they can be changed by SATsolver."""

        for pair in self.list_of_bps:
            if position in pair:
                bps = [x for x in pair if x != position][0]
                phrase = self.try_motif(seq[bps], [bps])
                for i in phrase:
                    if i in try_cnf:
                        try_cnf.remove(i)
        return try_cnf

    def get_random_mutation(self, seq):
        """Incorporates mutations into the sequence. The kind of mutation depends on
        a random choice and a probability value specified by a user.
        Function can incorporate single point mutations or mutations of a sequence segment,
        based on sequences from a library, also specified by a user."""

        seq = ''.join(seq)
        seq_first = seq

        # add condition -> seq has to be different from initial sequence
        phrase = self.create_formula_whole_sequence(seq, 0)
        phrase[0] = [-x for x in phrase[0]]
        self.cnf_mutations += phrase
        self.cnf_mutations = [list(x) for x in set(tuple(x) for x in self.cnf_mutations)]  # unique self.cnf_mutations

        # draw position and mutation nt
        constr = []
        [constr.append(x) for x in self.constraints]
        self.seq_options = self.convert_shortcuts(constr)

        while True:
            while True:
                random_position = random.randint(0, len(seq) - 1)
                if len(self.seq_options[random_position]) > 1:
                    break
            # check if drawn position in any motif
            in_motifs = []
            if self.motifs_locations:
                for motif_found in self.motifs_locations:
                    for i in range(1, len(motif_found[1:]), 2):
                        if (random_position >= motif_found[i]) and (random_position <= motif_found[i + 1]):
                            in_motifs.append(motif_found)
                            # draw what to do next with given probability
                            draw = choice([0, 1], 1, p=[self.probability, 1 - self.probability])
            # if position not found in motifs -> single mutation
            if not in_motifs:
                draw = 1
            # change whole fragment
            if draw == 0:
                # find sequences for the selected motif
                index = 0
                for element in self.sequences_to_check:
                    try:
                        if element[0] == in_motifs[0][0]:
                            index = self.sequences_to_check.index(element)
                            break
                    except ValueError:
                        pass
                # get random sequence for selected motif
                seqs_to_check = self.sequences_to_check[index][1:]
                while True:
                    if not seqs_to_check:
                        draw = 1  # make single mutation if cant incorporate any motif
                        break
                    random_seq = random.choice(seqs_to_check)
                    seqs_to_check.remove(random_seq)
                    if "&" in random_seq:
                        random_seq = random_seq.split("&")
                    # get formula for this sequence
                    try_cnf = self.cnf_mutations[:]
                    if type(random_seq) == list:
                        k = 1
                        for s in random_seq:
                            formula = self.try_motif(s, [in_motifs[0][k], in_motifs[0][k + 1]])
                            k += 2
                            try_cnf += formula
                    else:
                        formula = self.try_motif(random_seq, in_motifs[0][1:])
                        try_cnf += formula
                    solution = pycosat.solve(try_cnf)
                    if solution != "UNSAT":
                        if self.verbose == 1:
                            print "solution: ", solution
                        converted = self.solution_to_sequence_converter(solution)
                        converted = ''.join(converted)
                        number_of_mutations = sum(1 for a, b in zip(seq_first, converted) if a != b)
                        return converted
            # make single/double mutation
            if draw == 1:
                # freeze the whole sequence
                try_cnf = self.cnf_mutations[:]
                while True:
                    mutation = random.choice(self.bases)
                    if mutation != seq[random_position]:
                        break
                seq = list(seq)
                seq[random_position] = mutation
                "".join(seq)
                phrase = self.create_formula_whole_sequence(seq, 1)
                try_cnf += phrase
                solution = pycosat.solve(try_cnf)
                if solution == "UNSAT":
                    # check if nt at random position pairs with any other nts
                    cnf_unfreeze = self.unfreeze_basepairs(random_position, seq, try_cnf)
                    try_cnf = cnf_unfreeze[:]
                    try_cnf += self.cnf
                    solution = pycosat.solve(try_cnf)
                    if solution != "UNSAT":
                        if self.verbose == 1:
                            print "solution: ", solution
                        converted = self.solution_to_sequence_converter(solution)
                        converted = ''.join(converted)
                        number_of_mutations = sum(1 for a, b in zip(seq_first, converted) if a != b)
                        return converted
                    # if still UNSAT - unfreeze a random naighbor
                    neighbors = []
                    # append to list mutation's neighbors
                    try:
                        neighbors.append((seq[random_position - 1], random_position - 1))
                    except IndexError:
                        pass
                    try:
                        neighbors.append((seq[random_position + 1], random_position + 1))
                    except IndexError:
                        pass
                    while True:
                        # can't make any mutations
                        if not neighbors:
                            raise Exception("Cannot incorporate any single mutation.")
                        neighbors = list(set(neighbors))
                        random_neighbor = random.choice(neighbors)
                        neighbors.remove(random_neighbor)
                        phrase = self.try_motif(random_neighbor[0], [random_neighbor[1]])
                        for i in phrase:
                            if i in try_cnf:
                                try_cnf.remove(i)
                        cnf_unfreeze = self.unfreeze_basepairs(random_neighbor[1], seq, try_cnf)
                        try_cnf = cnf_unfreeze[:]
                        try:
                            neighbors.append((seq[random_neighbor[1] - 1], random_neighbor[1] - 1))
                        except IndexError:
                            pass
                        try:
                            neighbors.append((seq[random_neighbor[1] + 1], random_neighbor[1] + 1))
                        except IndexError:
                            pass

                        solution = pycosat.solve(try_cnf)
                        if solution != "UNSAT":
                            if self.verbose == 1:
                                print "solution: ", solution
                            converted = self.solution_to_sequence_converter(solution)
                            converted = ''.join(converted)
                            number_of_mutations = sum(1 for a, b in zip(seq_first, converted) if a != b)
                            return converted
                else:
                    if self.verbose == 1:
                        print "solution: ", solution
                    converted = self.solution_to_sequence_converter(solution)
                    converted = ''.join(converted)
                    number_of_mutations = sum(1 for a, b in zip(seq_first, converted) if a != b)
                    return converted

    def check_patterns(self):
        """Checks if most stable motifs and patterns specified by a user are present in initial structures."""

        data = self.get_data_from_file(self.set_filename)

        for i in data:
            if i.startswith(">"):
                self.most_stable_motifs.append(i[1:].strip())
                self.sequences_to_check.append([i[1:].strip()])
            else:
                self.sequences_to_check[-1].append(i.strip())

        self.motifs_found = []
        for j in self.most_stable_motifs:
            if "&" not in j:
                for k in self.structures:
                    if j in k:
                        pattern_ready = re.compile(re.escape(j))
                        result_start = [n.start() for n in re.finditer(pattern_ready, k)]
                        result_end = [n.end() - 1 for n in re.finditer(pattern_ready, k)]
                        for x in range(len(result_start)):
                            self.motifs_locations.append([j, result_start[x], result_end[x]])
                        self.motifs_found.append(j)
            elif "&" in j:
                patterns = j.split("&")
                patterns_locations = []
                result = self.check_brackets(''.join(patterns))
                if j.count("(") == j.count(")") and j.count("[") == j.count("]") and result:
                    for l in patterns:
                        for m in self.structures:
                            pattern = l
                            pattern_ready = re.compile(re.escape(pattern))
                            result_start = [n.start() for n in re.finditer(pattern_ready, m)]
                            result_end = [n.end() - 1 for n in re.finditer(pattern_ready, m)]
                            for x in range(len(result_start)):
                                patterns_locations.append((result_start[x], result_end[x]))
                    patterns_locations = set(patterns_locations)
                    patterns_locations = list(patterns_locations)

                    results = []
                    for structure in self.structures:
                        for pair in range(len(patterns_locations) - 1):
                            interior = structure[patterns_locations[pair][1] + 1:patterns_locations[pair + 1][0]]
                            result = self.check_brackets(interior)
                            results.append(result)
                    if all(results):
                        for i in range(0, len(patterns_locations), 2):
                            self.motifs_locations.append([j])
                            self.motifs_locations[-1].append(patterns_locations[i][0])
                            self.motifs_locations[-1].append(patterns_locations[i][1])
                            self.motifs_locations[-1].append(patterns_locations[i + 1][0])
                            self.motifs_locations[-1].append(patterns_locations[i + 1][1])
                        self.motifs_found.append(j)
                else:
                    raise ValueError("Wrong pattern. Check file again.")

    def try_motif(self, motif_seq, locations):
        """Returns a logic formula for a given motif.
        In this logic formula motif is described as constraints."""

        phrase = self.create_formula_whole_sequence(motif_seq, 1)
        temp = []
        for k in phrase:
            for j in k:
                if j < 0:
                    temp.append([-1 * (abs(j) + (locations[0] * 2))])
                else:
                    temp.append([j + locations[0] * 2])
        phrase = temp[:]
        return phrase

    def solution_to_sequence_converter(self, logic_list):
        """Converts a given solution to nucleotide sequence."""

        assert map(abs, logic_list) == range(1, 1 + len(logic_list))
        logic_values = map(lambda x: (0, 1)[x > 0], logic_list)
        dictionary = dict((v, k) for k, v in self.dict_logic_nt_values.iteritems())
        return map(lambda nt: dictionary[nt], zip(logic_values[0::2], logic_values[1::2]))

    def calling_function(self):
        """Calls all functions in the right order."""

        self.file_parser()
        data = self.get_data_from_file(self.filename)
        self.get_structures_from_data(data)
        self.check_input()
        result = self.check_brackets(self.structures)
        if result is False:
            raise ValueError("Wrong base pairings in structures. Check input file.")
        if self.set_filename is not None:
            self.check_patterns()
        self.check_forbidden_nts()
        solution = self.check_constraints_logic()
        if solution != "UNSAT":
            self.solutions.append(solution)
            if self.verbose == 1:
                print "solution: ", solution
            self.cnf_mutations = self.cnf[:]
        else:
            raise Exception("No solution found. Check input file again.")
        converted_sequence = self.solution_to_sequence_converter(solution)
        print converted_sequence
        seq = self.get_random_mutation(converted_sequence)
        print seq


def main():

    inverse_folding = coRNA()
    inverse_folding.calling_function()


if __name__ == "__main__":
    main()
