



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


orthogroups = open("Orthogroups.txt")
orthogroups = orthogroups.readlines()

orthogroup_distances = open("orthogroup_distances.tsv", "a")

##################################
#           YAHOUENSIS           #
##################################


for orthogroup in orthogroups:
 orthogroup = orthogroup.rstrip("\n")
 line = orthogroup.split(" ")
 line = [i for i in line if i]
 orthogroup = line[0].rstrip(":")
 yahouensis_TEs = [i for i in line[1:] if i.startswith("yahouensis")]
 yahouensis_TEs = [i.lstrip("yahouensis_") for i in yahouensis_TEs]
 identified_tes = [i for i in yahouensis_TEs if i.startswith("TE_")]
 unidentified_tes = [i for i in yahouensis_TEs if not i.startswith("TE_")]
 for i in unidentified_tes:
     name = i.split("#")[0].replace("_", ":")
     start = name.split(":")[1].split("..")[0]
     end = name.split(":")[1].split("..")[1]
     unidentified_dist = ["yahouensis", orthogroup] + [i for i in yahouensis_dist if i[0:3] == [name, start, end]][0]
     orthogroup_distances.write("\t".join(unidentified_dist) + "\n")

 for i in identified_tes:
     name = i.split("|")[0]
     start = i.split("|")[1].split("#")[0].split("_")[1].split("..")[0]
     end = i.split("|")[1].split("#")[0].split("_")[1].split("..")[1]
     identified_dist = ["yahouensis", orthogroup] + [i for i in yahouensis_dist if i[0:3] == [name, start, end]][0]
     orthogroup_distances.write("\t".join(identified_dist) + "\n")


##################################
#           IMPOLITA           #
##################################


for orthogroup in orthogroups:
 orthogroup = orthogroup.rstrip("\n")
 line = orthogroup.split(" ")
 line = [i for i in line if i]
 orthogroup = line[0].rstrip(":")
 impolita_TEs = [i for i in line[1:] if i.startswith("impolita")]
 impolita_TEs = [i.lstrip("impolita_") for i in impolita_TEs]
 identified_tes = [i for i in impolita_TEs if i.startswith("TE_")]
 unidentified_tes = [i for i in impolita_TEs if not i.startswith("TE_")]
 for i in unidentified_tes:
     name = i.split("#")[0].replace("_", ":")
     start = name.split(":")[1].split("..")[0]
     end = name.split(":")[1].split("..")[1]
     unidentified_dist = ["impolita", orthogroup] + [i for i in impolita_dist if i[0:3] == [name, start, end]][0]
     orthogroup_distances.write("\t".join(unidentified_dist) + "\n")
 for i in identified_tes:
     name = i.split("|")[0]
     start = i.split("|")[1].split("#")[0].split("_")[1].split("..")[0]
     end = i.split("|")[1].split("#")[0].split("_")[1].split("..")[1]
     identified_dist = ["impolita", orthogroup] + [i for i in impolita_dist if i[0:3] == [name, start, end]][0]
     orthogroup_distances.write("\t".join(identified_dist) + "\n")


##################################
#           SANDWICENSIS           #
##################################


for orthogroup in orthogroups:
 orthogroup = orthogroup.rstrip("\n")
 line = orthogroup.split(" ")
 line = [i for i in line if i]
 orthogroup = line[0].rstrip(":")
 sandwicensis_TEs = [i for i in line[1:] if i.startswith("sandwicensis")]
 sandwicensis_TEs = [i.lstrip("sandwicensis_") for i in sandwicensis_TEs]
 identified_tes = [i for i in sandwicensis_TEs if i.startswith("TE_")]
 unidentified_tes = [i for i in sandwicensis_TEs if not i.startswith("TE_")]
 for i in unidentified_tes:
     # need to do this to allow proper reformatting of the colon, will undo in following lines
     i = i.replace("Scaffolds_", "Scaffolds")
     name = i.split("#")[0].replace("_", ":")
     name = name.replace("Scaffolds", "Scaffolds_")
     start = name.split(":")[1].split("..")[0]
     end = name.split(":")[1].split("..")[1]
     unidentified_dist_attempt = [i for i in sandwicensis_dist if i[0:3] == [name, start, end]]
     if unidentified_dist_attempt:
         unidentified_dist = ["sandwicensis", orthogroup] + unidentified_dist_attempt[0]
         orthogroup_distances.write("\t".join(unidentified_dist) + "\n")
 for i in identified_tes:
     i = i.replace("Scaffolds_", "Scaffolds")
     name = i.split("|")[0]
     start = i.split("|")[1].split("#")[0].split("_")[1].split("..")[0]
     end = i.split("|")[1].split("#")[0].split("_")[1].split("..")[1]
     print(i)
     print(name)
     print(start)
     print(end)
     identified_dist_attempt = [i for i in sandwicensis_dist if i[0:3] == [name, start, end]]
     if identified_dist_attempt:
         identified_dist = ["sandwicensis", orthogroup] + identified_dist_attempt[0]
         orthogroup_distances.write("\t".join(identified_dist) + "\n")

orthogroup_distances.close()


