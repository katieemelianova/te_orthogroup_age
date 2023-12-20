



##################################
#           READ IN GFF          #
##################################


yahouensis_dist = []
impolita_dist = []
sandwicensis_dist = []


yahouensis_gff = open("yahouensis.fasta.mod.EDTA.intact.gff3")
yahouensis_gff = yahouensis_gff.readlines()
for i in yahouensis_gff:
    if "LTR" in i.split("\t")[2]:
        coords = i.split("\t")[3:5]
        info = i.split("\t")[8].split(";")
        name = [i.lstrip("Name=") for i in info if i.startswith("Name")]
        ltr_identity = [i.lstrip("ltr_identity=") for i in info if i.startswith("ltr_identity")]
        yahouensis_dist.append(name + coords + ltr_identity)

impolita_gff = open("impolita.fasta.mod.EDTA.intact.gff3")
impolita_gff = impolita_gff.readlines()
for i in impolita_gff:
    if "LTR" in i.split("\t")[2]:
        coords = i.split("\t")[3:5]
        info = i.split("\t")[8].split(";")
        name = [i.lstrip("Name=") for i in info if i.startswith("Name")]
        ltr_identity = [i.lstrip("ltr_identity=") for i in info if i.startswith("ltr_identity")]
        impolita_dist.append(name + coords + ltr_identity)


sandwicensis_gff = open("sandwicensis.fasta.mod.EDTA.intact.gff3")
sandwicensis_gff = sandwicensis_gff.readlines()
for i in sandwicensis_gff:
    if "LTR" in i.split("\t")[2]:
        coords = i.split("\t")[3:5]
        info = i.split("\t")[8].split(";")
        name = [i.lstrip("Name=") for i in info if i.startswith("Name")]
        ltr_identity = [i.lstrip("ltr_identity=") for i in info if i.startswith("ltr_identity")]
        sandwicensis_dist.append(name + coords + ltr_identity)


orthogroups = open("Orthogroups_noheader.tsv")
orthogroups = orthogroups.readlines()

orthogroup_distances = open("orthogroup_distances.tsv", "a")

##################################
#           YAHOUENSIS           #
##################################


for orthogroup in orthogroups:
 line = orthogroup.split("\t")
 line = [i for i in line if i]
 orthogroup = line[0]
 yahouensis_TEs = [i for i in line[1].split(", ") if i.startswith("yahouensis")]
 yahouensis_TEs = [i.lstrip("yahouensis_") for i in yahouensis_TEs]
 identified_tes = [i for i in yahouensis_TEs if i.startswith("TE_")]
 unidentified_tes = [i for i in yahouensis_TEs if not i.startswith("TE_")]
 for i in unidentified_tes:
     name = i.split("#")[0].replace("_", ":")
     start = name.split(":")[1].split("..")[0]
     end = name.split(":")[1].split("..")[1]
     unidentified_dist = ["yahouensis", orthogroup] + [i for i in yahouensis_dist if i[0:3] == [name, start, end]][0]
     print(unidentified_dist)
     orthogroup_distances.write("\t".join(unidentified_dist) + "\n")

 for i in identified_tes:
     name = i.split("|")[0]
     start = i.split("|")[1].split("#")[0].split("_")[1].split("..")[0]
     end = i.split("|")[1].split("#")[0].split("_")[1].split("..")[1]
     identified_dist = [["yahouensis", orthogroup] + [i for i in yahouensis_dist if i[0:3] == [name, start, end]][0]]

 
##################################
#           IMPOLITA           #
##################################

for orthogroup in orthogroups:
 line = orthogroup.split("\t")
 line = [i for i in line if i]
 orthogroup = line[0]
 impolita_TEs = [i for i in line[1].split(", ") if i.startswith("impolita")]
 impolita_TEs = [i.lstrip("impolita_") for i in impolita_TEs]
 identified_tes = [i for i in impolita_TEs if i.startswith("TE_")]
 unidentified_tes = [i for i in impolita_TEs if not i.startswith("TE_")]
 for i in unidentified_tes:
     name = i.split("#")[0].replace("_", ":")
     start = name.split(":")[1].split("..")[0]
     end = name.split(":")[1].split("..")[1]
     unidentified_dist = [["impolita", orthogroup] + [i for i in impolita_dist if i[0:3] == [name, start, end]][0]]
     print(unidentified_dist)
 for i in identified_tes:
     name = i.split("|")[0]
     start = i.split("|")[1].split("#")[0].split("_")[1].split("..")[0]
     end = i.split("|")[1].split("#")[0].split("_")[1].split("..")[1]
     identified_dist = [["impolita", orthogroup] + [i for i in impolita_dist if i[0:3] == [name, start, end]][0]]


##################################
#           SANDWICENSIS           #
##################################


for orthogroup in orthogroups:
 line = orthogroup.split("\t")
 line = [i for i in line if i]
 orthogroup = line[0]
 sandwicensis_TEs = [i for i in line[1].split(", ") if i.startswith("sandwicensis")]
 sandwicensis_TEs = [i.lstrip("sandwicensis_") for i in sandwicensis_TEs]
 identified_tes = [i for i in sandwicensis_TEs if i.startswith("TE_")]
 unidentified_tes = [i for i in sandwicensis_TEs if not i.startswith("TE_")]
 for i in unidentified_tes:
     name = i.split("#")[0].replace("_", ":")
     start = name.split(":")[1].split("..")[0]
     end = name.split(":")[1].split("..")[1]
     unidentified_dist = [["sandwicensis", orthogroup] + [i for i in sandwicensis_dist if i[0:3] == [name, start, end]][0]]
     print(unidentified_dist)
 for i in identified_tes:
     name = i.split("|")[0]
     start = i.split("|")[1].split("#")[0].split("_")[1].split("..")[0]
     end = i.split("|")[1].split("#")[0].split("_")[1].split("..")[1]
     identified_dist = [["sandwicensis", orthogroup] + [i for i in sandwicensis_dist if i[0:3] == [name, start, end]][0]]






