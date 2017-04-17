#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np

# generate A from the textr
# A should be 2D list

def _parse_input_text(path):
    ret_list = []
    with open(path, 'r') as inputfile:
        for line in inputfile:
            inner_list = []
            text_line = line.strip().split(',')
            for val in text_line:
                inner_list.append(int(val))
            ret_list.append(inner_list)
            #print(line.strip().split(','))
        return ret_list

def main():
    # parse input text values
    file_path = 'path.txt'
    parsed_values = []
    parsed_values = _parse_input_text(file_path)
    #print('Parsed Values: '+str(parsed_values))
    A = np.asarray(parsed_values)
    plt.figure(1)
    plt.imshow(A, interpolation='nearest', cmap='gray')
    plt.grid(False)
    plt.show()

if __name__ == '__main__':
    main()