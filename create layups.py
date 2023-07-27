import itertools
import numpy as np
from collections import defaultdict
import csv
orientation_list = [0, 45, 90, -45]
N = 10
test_layups = []
layer_counts = defaultdict(int)
check_symmetrical = True  # Set to True to enforce symmetrical layups
check_balanced = True  # Set to True to enforce balanced layups
check_outermost = True  # Set to True to enforce that outermost ply condition cannot be 0 or 90 deg
check_adjacent = True  # Set to True to enforce that adjacent plies cannot have an angle difference greater or equal to 90 deg
for nr_of_layers in range(2, N+1):
    for combination in itertools.product(orientation_list, repeat=nr_of_layers):
        valid_layup = True
        if check_symmetrical and combination != combination[::-1]:
            valid_layup = False
        if check_balanced and combination.count(45) != combination.count(-45):
            valid_layup = False
        if check_outermost and (combination[0] == 0 or combination[0] == 90 or combination[-1] == 0 or combination[-1] == 90):
            valid_layup = False
        if check_adjacent:
            for i in range(1, len(combination)):
                if (combination[i] == 90 and combination[i-1] == 0) or (combination[i] == 0 and combination[i-1] == 90):
                    valid_layup = False
                    break
                if (combination[i] == 45 and combination[i-1] == -45) or (combination[i] == -45 and combination[i-1] == 45):
                    valid_layup = False
                    break
        if valid_layup:
            test_layups.append(combination)
print(test_layups)