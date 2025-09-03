import os

def parse_alignment_file(filename, stats_filename):
    var1 = None
    var2 = None
    var3 = None
    var4 = None
    var5 = None
    var6 = None
    var7 = None
    var8 = None
    N_reads = None
    total_residues = None

    with open(filename, 'r') as file:
        lines = file.readlines()

        for line in lines:
            if 'reads' in line:
                try:
                    N_reads = int(line.split()[0])
                except ValueError:
                    continue

            if 'aligned concordantly exactly 1 time' in line:
                var1 = int(line.split()[0])
                var3 = var1 * 2
            elif 'aligned concordantly >1 times' in line:
                var2 = int(line.split()[0])
                var4 = var2 * 2
            elif 'aligned exactly 1 time' in line:
                var5 = int(line.split()[0])
            elif 'aligned >1 times' in line:
                var6 = int(line.split()[0])
            elif 'aligned discordantly 1 time' in line:
                var7 = int(line.split()[0])
                var8 = var7 * 2

    alignment_name = filename.split('--VS--')[0].strip()
    alignment_name = alignment_name.split('/')[-1].split('--VS--')[0]

    with open(stats_filename, 'r') as stats_file:
        stats_lines = stats_file.readlines()

        for line in stats_lines:
            if alignment_name in line:
                total_residues = int(line.split("Total # residues:")[1])
                break

    return var1, var2, var3, var4, var5, var6, var7, var8, N_reads, total_residues

# Specify the directory path
directory = '/home/rohit/rohit/grossart_things/hisat2/FUNEPA_new_raw'

# Specify the path to the stats file
stats_file = "/home/rohit/rohit/grossart_things/hisat2/FUNEPA_new_raw/mag_size_and
_names.txt"

# Specify the path for the output folder
output_folder = "/home/rohit/rohit/grossart_things/hisat2/FUNEPA_new_raw/RPMM_outp
ut_files"

# Create the output folder if it doesn't exist
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Iterate over all files in the directory
for filename in os.listdir(directory):
    if filename.endswith('.fasta.gz.summary.txt'):
        # Generate the full file paths
        file_path = os.path.join(directory, filename)
        output_file = os.path.join(output_folder, f"{filename}_parsed.txt")

        # Parse the alignment file and get the results
        var1, var2, var3, var4, var5, var6, var7, var8, N_reads, total_residues =
parse_alignment_file(file_path, stats_file)

        # Perform additional calculations
        num_reads = N_reads * 2
        total_reads_in_millions = num_reads / 1000000
        genome_length = total_residues
        genome_length_in_ml = total_residues / 1000000
        numhits = var3 + var4 + var5 + var6 + var8
        numhits_per_million = numhits / total_reads_in_millions
        RPKM = numhits_per_million / genome_length_in_ml

        # Write the results to the output file
        with open(output_file, 'w') as output:
            #output.write(f"Number of reads: {num_reads}\n")
            #output.write(f"Total reads in millions: {total_reads_in_millions}\n")
            #output.write(f"Genome length: {genome_length}\n")
            #output.write(f"Genome length in kb: {genome_length_in_kb}\n")
            #output.write(f"Numhits: {numhits}\n")
            #output.write(f"Numhits per million: {numhits_per_million}\n")
            output.write(f"RPKM: {RPKM}\n")

        print(f"Output saved to: {output_file}")
