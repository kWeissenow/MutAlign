from Bio import SeqIO
import numpy as np
import string
import os
import sys


def parse_vespa_output(filename):
    aa_list = list("ACDEFGHIKLMNPQRSTVWY")

    # scan length of file
    with open(filename, "r") as f:
        for line in f:
            pass
        last_index = line.strip().split(";")[0][1:-1]
        length = int(last_index) + 1

    # fill mutation matrix
    mut_matrix = np.full((length, 20), 0.5)
    with open(filename, "r") as f:
        for line in f:
            # skip header line
            if line.startswith("Mutant"):
                continue
            id,score = line.strip().split(";")
            from_aa = id[0]
            to_aa = id[-1]
            idx = int(id[1:-1])
            mut_matrix[idx, aa_list.index(to_aa)] = float(score)

    return mut_matrix


def gradient_color(minval, maxval, val, color_palette=((0,0,255), (255,255,255), (255,0,0))):
    """ Computes intermediate RGB color of a value in the range of minval
        to maxval (inclusive) based on a color_palette representing the range.
    """
    max_index = len(color_palette)-1
    delta = maxval - minval
    if delta == 0:
        delta = 1
    v = float(val-minval) / delta * max_index
    i1, i2 = int(v), min(int(v)+1, max_index)
    (r1, g1, b1), (r2, g2, b2) = color_palette[i1], color_palette[i2]
    f = v - i1
    return "#{:02x}{:02x}{:02x}".format(int(r1 + f*(r2-r1)), int(g1 + f*(g2-g1)), int(b1 + f*(b2-b1)))


def main():
    out_dir = "./html/"
    if len(sys.argv) > 1:
        out_dir = sys.argv[1]

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    # parse identifiers
    identifiers = []
    with open("identifiers.txt", "r") as f:
        for line in f:
            identifiers.append(line.strip())
    identifiers = sorted(identifiers)

    for identifier in identifiers:
        print("Processing: " + identifier)

        # read MSA
        msa_records = list(SeqIO.parse("alignment_files/" + identifier + ".a3m", "fasta"))
        num_sequences = len(msa_records)

        # parse VESPA output and create mutation matrices
        mut_matrices = []
        for i,record in enumerate(msa_records):
            mut_matrix = parse_vespa_output("vespa_predictions/" + record.id + ".csv")
            if i != 0:
                # for all sequences (except query/human), remove non-aligned parts
                parts = record.description.split()
                tStart = int(parts[7])
                tEnd = int(parts[8])
                mut_matrix = mut_matrix[tStart:tEnd+1]
            mut_matrices.append(mut_matrix)

        # we are building a combined matrix for all sequences in the MSA (removing insertions)
        # determine total size of combined matrix
        max_len = 0
        table = str.maketrans('', '', string.ascii_lowercase)
        for record in msa_records:
            length = len(str(record.seq).translate(table))
            if length > max_len:
                max_len = length
        combined_matrix = np.full((max_len, 20 * num_sequences), 0.5)

        # iterate over rows of combined matrix and fill in values from individual mutation matrices
        combined_index = 0
        running_indices = [0 for i in range(num_sequences)]
        matrix_indices = [0 for i in range(num_sequences)]
        sequence_characters = ["" for i in range(num_sequences)]
        for combined_index in range(max_len):
            for i in range(num_sequences):
                while running_indices[i] < len(msa_records[i].seq) and msa_records[i].seq[running_indices[i]].islower():
                    running_indices[i] += 1
                    matrix_indices[i] += 1
                if running_indices[i] < len(msa_records[i].seq):
                    sequence_characters[i] += msa_records[i].seq[running_indices[i]]
                    if msa_records[i].seq[running_indices[i]] != "-":
                        combined_matrix[combined_index, i*20:i*20+20] = mut_matrices[i][matrix_indices[i], :]
                        matrix_indices[i] += 1
                    running_indices[i] += 1

        # render combined matrix as HTML files
        with open(os.path.join(out_dir, identifier + ".html"), "w") as f:
            f.write("<!doctype html><html><head><title>" + identifier + "</title></head>\n<body><table>")

            f.write("<tr><td>#</td>")
            for i in range(num_sequences):
                f.write("<td colspan=22 align=center>" + msa_records[i].id + "</td>")
            f.write("</tr>")

            for combined_index in range(max_len):
                f.write("<tr><td>" + str(combined_index+1) + "</td>")
                for i in range(num_sequences):
                    # write amino-acid character (or gap)
                    character = "-"
                    if combined_index < len(sequence_characters[i]):
                        character = sequence_characters[i][combined_index]
                    f.write("<td>" + character + "</td>")

                    # compute delta to human for non-human sequences (if position not gapped)
                    if i != 0 and character != "-":
                        delta = np.sum(np.abs(combined_matrix[combined_index, 0:20] - combined_matrix[combined_index, i*20:i*20+20]))
                        f.write("<td width=\"20\" style=\"background-color: " + gradient_color(0.0, 20.0, delta, color_palette=((255,255,255), (30,30,30), (25,25,25), (0,0,0))) + "\">&nbsp;</td>")
                    else:
                        f.write("<td>&nbsp;</td>")

                    # write cells with individual SAV prediction values
                    for idx in range(20):
                        f.write("<td style=\"background-color: " + gradient_color(0.0, 1.0, combined_matrix[combined_index, i*20+idx]) + "\">&nbsp;</td>")
                f.write("</tr>")

            f.write("</table></html>")

    # write index HTML file
    with open(os.path.join(out_dir, "index.html"), "w") as f:
        f.write("<!doctype html><html><head><title>VESPAl predictions for Jonathan proteins</title></head><body><div align=center><p><b>VESPAl predictions for Jonathan proteins</b></p><p>")
        for id in identifiers:
            f.write("<a href=\"{}.html\">{}</a><br/>".format(id, id))
        f.write("</p></div></body></html>")


if __name__ == "__main__":
    main()
