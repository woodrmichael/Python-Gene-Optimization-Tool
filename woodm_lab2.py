"""
CSC 2611 131
Lab 2 - Gene Optimization
Michael Wood
"""


def optimize_gene(gene_id, fasta_file_path, codon_freq_table_file_path):
    """
    optimize_gene takes in a file path to a fasta file with gene sequences, a file path
    to a codon frequency table, and a specific gene id. It will find that gene id and optimize
    it by replacing its codons with the most common codon resulting in the same amino acid based
    on the codon frequency table. This function will return the new gene sequence,
    if the gene id wasn't found, it will return -1.

    :param gene_id: the gene id to search for.
    :param fasta_file_path: a file path to a fasta file which contains a key that follows a '>'
    followed by its gene sequence.
    :param codon_freq_table_file_path: a file path to the codon frequency table which starts
    with a line 'Codon,Amino Acid,Frequency' followed by
    codon, amino acid, frequency data in that order.
    :return: the optimized gene sequence or -1 if the gene id wasn't found.
    """
    codon_list = get_gene_sequence(gene_id, fasta_file_path)
    if codon_list != -1:
        dictionary_list = create_dictionaries(codon_freq_table_file_path, "optimize")
        codon_to_amino_acid = dictionary_list[0]
        amino_acid_to_codon = dictionary_list[1]
        optimized_gene_sequence = convert_gene_sequence(
            codon_list, codon_to_amino_acid, amino_acid_to_codon)
        create_gene_sequence_file(gene_id, "optimized", optimized_gene_sequence)
        return optimized_gene_sequence
    else:
        return -1


def deoptimize_gene(gene_id, fasta_file_path, codon_freq_table_file_path):
    """
    deoptimize_gene takes in a file path to a fasta file with gene sequences, a file path
    to a codon frequency table, and a specific gene id. It will find that gene id and deoptimize
    it by replacing its codons with the least common codon resulting in the same amino acid based
    on the codon frequency table. This function will return the new gene sequence,
    if the gene id wasn't found, it will return -1.

    :param gene_id: the gene id to search for.
    :param fasta_file_path: a file path to a fasta file which contains a key that follows a '>'
    followed by its gene sequence.
    :param codon_freq_table_file_path: a file path to the codon frequency table which starts
    with a line 'Codon,Amino Acid,Frequency' followed by
    codon, amino acid, frequency data in that order.
    :return: the deoptimized gene sequence or -1 if the gene id wasn't found.
    """
    codon_list = get_gene_sequence(gene_id, fasta_file_path)
    if codon_list != -1:
        dictionary_list = create_dictionaries(codon_freq_table_file_path, "deoptimize")
        codon_to_amino_acid = dictionary_list[0]
        amino_acid_to_codon = dictionary_list[1]
        deoptimized_gene_sequence = convert_gene_sequence(
            codon_list, codon_to_amino_acid, amino_acid_to_codon)
        create_gene_sequence_file(gene_id, "deoptimized", deoptimized_gene_sequence)
        return deoptimized_gene_sequence
    else:
        return -1


def get_gene_sequence(gene_id, fasta_file_path):
    """
    get_gene_sequence gets the gene sequence from the fasta file given a gene id
    and returns it as a list with each element being a separate codon in order of the
    original gene sequence.

    :param gene_id: the gene id to search for.
    :param fasta_file_path: a file path to a fasta file which contains a key that follows a '>'
    followed by its gene sequence.
    :return: a list where each element is a separate codon in the order of the original gene sequence
    or -1 if the gene id couldn't be found.
    """
    fasta_file = list(open(fasta_file_path))
    for i in range(len(fasta_file)):
        fasta_file[i] = fasta_file[i].rstrip()
    if ">" + gene_id in fasta_file:
        gene_sequence_index = fasta_file.index(">" + gene_id) + 1
        gene_sequence = fasta_file[gene_sequence_index]
        codon_list = []
        for i in range(0, len(gene_sequence), 3):
            codon_list.append(gene_sequence[i:i + 3])
        return codon_list
    else:
        return -1


def create_dictionaries(codon_freq_table_file_path, swap_type):
    """
    create_dictionaries creates the dictionaries used to convert the original gene sequence
    to either an optimized or deoptimized version. Returns the dictionaries as a list where the first
    item in the list is a dictionary that maps codon to amino acid and the second item in the list
    is a dictionary that maps amino acid to a list that contains 2 elements, the frequency of the
    optimized/deoptimized codon as well as that specific codon.

    :param codon_freq_table_file_path: a file path to the codon frequency table which starts
    with a line 'Codon,Amino Acid,Frequency' followed by
    codon, amino acid, frequency data in that order.
    :param swap_type: 'optimize' or 'deoptimize' depending on which operation is being performed.
    :return: a list where the first item in the list is a dictionary that maps
    codon to amino acid and the second item in the list is a dictionary that maps
    amino acid to a list that contains 2 elements,
    the frequency of the optimized/deoptimized codon as well as that specific codon.
    """
    ecol_codon_freq_file = list(open(codon_freq_table_file_path))
    codon_to_amino_acid = {}
    amino_acid_to_codon = {}
    for i in range(1, len(ecol_codon_freq_file)):
        line = ecol_codon_freq_file[i].rstrip().split(",")
        codon_to_amino_acid[line[0]] = line[1]
        if ((line[1] not in amino_acid_to_codon.keys())
                or (amino_acid_to_codon[line[1]][0] < float(line[2]) and swap_type == "optimize")
                or (amino_acid_to_codon[line[1]][0] > float(line[2]) and swap_type == "deoptimize")):
            amino_acid_to_codon[line[1]] = (float(line[2]), line[0])
    return [codon_to_amino_acid, amino_acid_to_codon]


def convert_gene_sequence(codon_list, codon_to_amino_acid, amino_acid_to_codon):
    """
    convert_gene_sequence converts the original gene sequence to the newly optimized/deoptimized
    gene sequence. This function uses the codon to amino acid and amino acid to codon dictionaries
    to convert the original list of codons to a string representation of the newly
    optimized/deoptimized gene sequence.

    :param codon_list: a list of codons in order of the original gene sequence.
    :param codon_to_amino_acid: a dictionary that maps codon to amino acid.
    :param amino_acid_to_codon: a dictionary that maps amino acid to a list that contains 2 elements, the frequency of the
    optimized/deoptimized codon as well as that specific codon.
    :return: a string representation of the newly optimized/deoptimized gene sequence.
    """
    amino_acid_list = []
    for i in range(len(codon_list)):
        amino_acid_list.append(codon_to_amino_acid[codon_list[i]])
    new_codon_list = []
    for i in range(len(amino_acid_list)):
        new_codon_list.append(amino_acid_to_codon[amino_acid_list[i]][1])
    new_gene_sequence = ""
    for i in range(len(new_codon_list)):
        new_gene_sequence += new_codon_list[i]
    return new_gene_sequence


def create_gene_sequence_file(gene_id, swap_type, gene_sequence):
    """
    create_gene_sequence_file writes a file named <gene_id>_<swap_type>.fasta which
    contains the gene_id preceded by a '>' on the first line
    and the gene_sequence on the second line.

    :param gene_id: the id of the gene that was optimized/deoptimized
    :param swap_type: 'optimize' or 'deoptimize' depending on which operation is being performed.
    :param gene_sequence: the newly optimized/deoptimized gene sequence.
    :return: None
    """
    file = open(gene_id + "_" + swap_type + ".fasta", "w")
    file.write(">" + gene_id + "\n")
    file.write(gene_sequence)
    file.close()


def main():
    running = True
    while running:
        fasta_file_path = str(input("Please enter a .fasta file path: "))
        codon_freq_table_file_path = str(input("Please enter a codon frequency file path: "))
        gene_id = str(input("Please enter a gene ID: "))
        swap_type = str(input("Would you like to 'optimize' or 'deoptimize' this gene sequence: "))
        if swap_type == "optimize":
            if optimize_gene(gene_id, fasta_file_path, codon_freq_table_file_path) != -1:
                print("Gene sequence is optimized. Gene sequence has been saved to",
                      gene_id + "_" + swap_type + "d.fasta\n")
            else:
                print("Gene ID does not exist. Try again.\n")
        elif swap_type == "deoptimize":
            if deoptimize_gene(gene_id, fasta_file_path, codon_freq_table_file_path) != -1:
                print("Gene sequence is deoptimized. Gene sequence has been saved to",
                      gene_id + "_" + swap_type + "d.fasta\n")
            else:
                print("Gene ID does not exist. Try again.\n")
        else:
            print("Please enter either 'optimize' or 'deoptimize'. Try again.\n")
        running = str(input("Would you like to continue? 'y' to continue, any other key to cancel: ")) == 'y'


if __name__ == "__main__":
    main()
