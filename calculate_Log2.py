import os
import numpy as np
import math

def read_matrix_from_file(matrix_file_path):
    matrix = []
    with open(matrix_file_path, 'r') as file:
        lines = file.readlines()

        for line in lines:
            row = line.strip().split()
            name, values = row[0], row[1:]
            row = [float(value) if float(value) >= 1 else 1 for value in values]
            matrix.append((name, row))
    return matrix

def calculate_log_of_matrix(matrix):
    log_matrix = []
    for item in matrix:
        name, row = item
        log_row = [math.log2(value) for value in row]  # Calculate log base 2
        log_matrix.append((name, log_row))
    return log_matrix

def main(matrix_file_path, output_file_path):
    # Read the matrix from the file
    input_matrix = read_matrix_from_file(matrix_file_path)

    # Calculate log base 2 for the numeric values in the matrix
    log_matrix = calculate_log_of_matrix(input_matrix)

    # Write the results to the output file
    with open(output_file_path, 'w') as output_file:
        for item in log_matrix:
            name, row = item
            non_numeric_chars = name.replace(' ', '\t')  # Convert spaces to tabs
            log_row_str = "\t".join(str(value) for value in row)
            output_file.write(f"{non_numeric_chars}\t{log_row_str}\n")

if __name__ == "__main__":
    # Specify the input and output file paths
    matrix_file_path = "/home/rohit/rohit/grossart_things/hisat2/FUNEPA_new_raw/RPMM_
output_files/RPMM_matrix"
    output_file_path = "/home/rohit/rohit/grossart_things/hisat2/FUNEPA_new_raw/RPMM_
output_files/log2_rpmm"

    main(matrix_file_path, output_file_path)
