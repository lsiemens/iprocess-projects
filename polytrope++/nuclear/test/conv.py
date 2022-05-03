"""Convert data copied from tables NACRE II to files of the form
(T9 rate rate_min rate_max)

Xu, Y. al, et al. "NACRE II: an update of the NACRE compilation of
    charged-particle-induced thermonuclear reaction rates for nuclei
    with mass number A< 16." Nuclear Physics A 918 (2013): 61-169.
"""

fnames = ["he3-he3pp-he4_NACREII", "n15-pg-o16_NACREII", "b11-an-n14_NACREII"]

for fname in fnames:
    with open(fname, "r") as fin:
        text = fin.read()

    text = text.replace("−", "-")
    text = text.split("\n")
    data = ["#NACRE II: " + fname.split("_")[0]]

    if text[0].startswith("#NACRE"):
        print(f"{fname}: Already converted")
        continue

    for line in text:
        if line.strip() != "":
            line = line.split()
            data.append(f"{line[0]} {line[1]} {line[2]} {line[3]}")

    for line in text:
        if line.strip() != "":
            line = line.split()
            try:
                data.append(f"{line[4]} {line[5]} {line[6]} {line[7]}")
            except IndexError:
                pass

    with open(fname, "w") as fout:
        fout.write("\n".join(data))
