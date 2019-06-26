###################################################################################
#
#							seq processiong functions
#
###################################################################################

#############################################################################
#extract secondary linked node of nodes from whole network
#parental_network should be a list
#nodes of the parental_network are in first two columns

#parental_network = read_lines_to_list(handle,sep_string="\t")

#take out the row contains target string
def search_linked_nodes(seed_nodes,parental_network):
	linked_nodes = []
	for node in seed_nodes:
		for lines in parental_network:
			if node in lines[0:2]:
				linked_nodes.append(tuple(lines[0:2]))
	return linked_nodes
#search for the nodes linked to nodes that linked to seed_nodes
#find nodes in linked_nodes that are not in seed_nodes
#used them as seed_nodes search for nodes linked to them
#return the results
def search_linked_secondary_nodes(seed_nodes,parental_network):
	linked_nodes = []
	first_layer_nodes = []
	raw_first_layer_nodes = search_linked_nodes(seed_nodes,parental_network)
	for lines in raw_first_layer_nodes:
#		if lines[0] in seed_nodes && lines[1] in seed_nodes:
#			continue
		if lines[0] in seed_nodes:
			first_layer_nodes.append(lines[1])
		elif lines[1] in seed_nodes:
			first_layer_nodes.append(lines[0])
		else:
			print("error, node don't linked to seed_nodes appeared")
	linked_nodes = search_linked_nodes(first_layer_nodes,parental_network)
	return linked_nodes

#########################################################################
#network sif file rename (first column is node name)
# replace first column by value in a dictionary
def sif_rename(sif_handle,dict_handle,sep_string="\t",sepstring="\t"):
	sif_file = read_table(sif_handle,sepstring)
	dict_list = read_lines_to_dict(dict_handle,sep_string)
	out_file = []
	for i in sif_file:
		if i[1] in dict_list.keys():
			out_file.append(dict_list[i[1]])
		else:
			out_file.append(i)
	return out_file



####################################################################################
#judge if a variable exists
def isset(v):
	try :
		type (eval(v))
	except :
		return  0
	else :
		return  1
######################################################
#read fasta sequence into a list
def read_fasta2dict(handle):
	i = 0
	matr = {}
	while True:
		line = handle.readline()
		if line == "":
			i +=1
			if i>100:
				print ("100 empty line")
				break
		if not line:
			break
		if line[0] =='>':
			break
	while True:
		if not line:
			return matr
		if line[0] != '>':
			print ("Records in Fasta files should start with '>' character")
			break
		if line[0] == '>':
			title = line[1:].rstrip()
			lines = []
			line = handle.readline()
			while True:
				if not line:
					break
				if line[0] =='>':
					break
				lines.append(line.rstrip())
				line = handle.readline()
			matr[title]=("".join(lines).replace(" ", "").replace("\r", ""))
######################################################
#read fasta sequence into an list
def read_fasta(handle):
	i = 0
	matr = []
	while True:
		line = handle.readline()
		if line == "":
			i +=1
			if i>100:
				print ("100 empty line")
				break
		if not line:
			break
		if line[0] =='>':
			break
	while True:
		if not line:
			return matr
		if line[0] != '>':
			print ("Records in Fasta files should start with '>' character")
			break
		if line[0] == '>':
			title = line[1:].rstrip()
			lines = []
			line = handle.readline()
			while True:
				if not line:
					break
				if line[0] =='>':
					break
				lines.append(line.rstrip())
				line = handle.readline()
			matr.append((title, "".join(lines).replace(" ", "").replace("\r", "")))
#############################################################################################
#delete duplicated sequences
#input should be [[title,seq],[title,seq],...,[title,seq]]
#outname is a sting file handle
def del_dup(seq_list,out_name):
	i=0
	a=len(seq_list)
	while True:
		if i >= a:
			break
		else:
			if isset(out_name):
				out_name.write(str(seq_list[i][0])+'_')
			n = i + 1

		while True:

			if n >= a:
				break
			if seq_list[i][1] == seq_list[n][1]:

				del seq_list[n]
				a -= 1
				print ("OK!")

			else:
				n +=1
		a = len(seq_list)
		i += 1
	return seq_list
#################################################################################
#delete similar sequences base on Bio.pairwise.alignment scores
def del_sim_pa(seq,iden_ran,out_n):
	from Bio import pairwise2
	i=0
	a=len(seq)
	scor_por = 0.0
	while True:
		if i >= a:
			break
		n = i + 1
		print (i)
	#	print a
		while True:
	#		print a
			print (n)
			print (seq[i][1])
			if n > a-1:
				break
			seq_l=[len(seq[i][1]),len(seq[n][1])]
			align=pairwise2.align.globalxx(seq[i][1],seq[n][1],score_only=1)
			print (align)
			scor_por = float(align/min(seq_l))
			if  scor_por > iden_ran:
				if seq_l.index(min(seq_l)) == 0:
					out_n.write('%s\n' % str(seq[i][0]))
					del seq[i]

					n = i + 1
					a -=1
					print ("del f")
				if seq_l.index(min(seq_l)) == 1:
					out_n.write('%s\n' % str(seq[n][0]))
					del seq[n]
					a -= 1
					print ("del l")

	#			print n
	#			print a
	#			print i
			else:
				n += 1
		a = len(seq)
		i += 1

###############read lines of strings into a list ######################
#read lines in a file into a list
#sepatate every line into several column by sep_strings(like ","","\t"... )
def read_lines_to_list(handle,sep_string="No"):
	i = 0
	out_list = []
	if sep_string=="No":
		while True:
			line = handle.readline()
			if line == "":
				i +=1
				if i>100:
					print ("100 empty line")
					break
			if not line:
				break
			if line !="":
				break
		out_list.append(line.rstrip())
		while True:
			if not line:
				return out_list
				break
			else:
				if line == "":
					return out_list
					break
				else:
					line = []
					line = handle.readline()
					out_list.append(line.rstrip())
		return out_list
	else:
		while True:
			line = handle.readline()
			if line == "":
				i +=1
				if i>100:
					print ("100 empty line")
					break
			if not line:
				break
			if line !="":
				break
		out_list.append(line.rstrip().split(sep_string))
		while True:
			if not line:
				return out_list
				break
			else:
				if line == "":
					return out_list
					break
				else:
					line = []
					line = handle.readline()
					out_list.append(line.rstrip().split(sep_string))
		return out_list

###############read two column file into a dictionary##################
#input should be text file with two separated columns
#or only the first two columns will be readed

def read_lines_to_dict(handle,sep_string,head="T"):
	i=0
	out_dict = {}
	while True:
		line = handle.readline()
		if line == "":
			i +=1
			if i>100:
				print ("100 empty line")
				break
		if not line:
			break
		if line !="":
			break
	if head == "T":
		out_dict["head"]=line.rstrip().split(sep_string)
	else:
		out_dict[line.rstrip().split(sep_string)[0]]=line.rstrip().split(sep_string)[1:]
	while True:
		line = []
		line = handle.readline()
		if not line:
			return out_dict
			break
		else:
			if line == "":

				continue
			else:
				out_dict[line.rstrip().split(sep_string)[0]]=line.rstrip().split(sep_string)[1:]
	return out_dict
################################## read text file  list#######################
#iput is files contains multi rows and columns
#seperated by "\t" "," etc.
#return is [[col1,col2,...,coln],[col1,col2,...,coln],...,[col1,col2,...,coln]]
def read_table(handle,sep_string):
	data_sheet = []
	i = 0
	while True:
		line = handle.readline()
		if line == "":
			i +=1
			if i>100:
				print ("100 empty line")
				break
		if not line:
			break
		if line !="":
			break
	data_sheet.append(line.rstrip().split(sep_string))
	while True:
		line = []
		line = handle.readline()
		if not line:
			return data_sheet
			break
		else:
			if line == "":
				continue
			else:
				data_sheet.append(line.rstrip().split(sep_string))
	return data_sheet

############## BIC gene cluser list to DAVID multi list file ################
#BIC_skmeans kgg file as input
#DAVID multi list format file as output
#def Bsk_kgg2DAVID_ml(handle):
#	import sys
#	sys.path.append("D:\\Program Files\\Anaconda3\\Lib")
#	import pandas
#	string_list = read_table(handle,"\t")
#


###########################################################################

#####################excel csv######################################





################Evolutionary trees##################################
#change nwk file into names
#use regular expression
#replace all the round brackets and commas ....
#
def read_nwk(handle):
	line = handle.readline().rstrip().rstrip("\n")
	return line


def nwk2name(lin_input):
	import re
#	out_m = lin_input.replace("\(","").replace("\)","").replace(",","\n").replace(":-0",":0").replace(":0.\d*.\d*:0.\d*.\d*","").replace(":0.\d*.\d*","").replace("'","")
#	out_m = sub("\(","",)
	out_m = re.sub("\(","",lin_input)
	out_m = re.sub("\)","",out_m)
	out_m = re.sub(",","\n",out_m)
	out_m = re.sub(":-0",":0",out_m)
	out_m = re.sub(":0.\d*.\d*:0.\d*.\d*","",out_m)
	out_m = re.sub(":0.\d*.\d*","",out_m)
	out_m = re.sub("'","",out_m)

	return out_m





################deal with files contain mutiple columns
#file name could contain path
def readFileIntoList(filename,sep=" "):
    contentList = map(lambda x:x.rstrip().split(sep), open(filename).readlines())
    return list(contentList)

#only return the first column
def readFile1thColIntoList(filename,sep=" "):
    contentList = map(lambda x:x.rstrip().split(sep)[0], open(filename).readlines())
    return list(contentList)

#This only return the first column of a filename
#Can be used to read GeneID
def readFileIntoUniqueList(filename):
    contentList = map(lambda x:x.rstrip().split(), open(filename).readlines())
    return list(set(map(tuple,contentList)))

def readLinesIntoList(filename):
    contentList = open(filename).readlines()
    return contentList



class ZipExhausted(Exception):
    pass
def izip_longest(*args, **kwds):
    # izip_longest('ABCD', 'xy', fillvalue='-') --> Ax By C- D-
    fillvalue = kwds.get('fillvalue')
    counter = [len(args) - 1]
    def sentinel():
        if not counter[0]:
            raise ZipExhausted
        counter[0] -= 1
        yield fillvalue
    fillers = repeat(fillvalue)
    iterators = [chain(it, sentinel(), fillers) for it in args]
    try:
        while iterators:
            yield tuple(map(next, iterators))
    except ZipExhausted:
        pass
#find the differences between two list
#input is two list of strings
#output is the interscection of two lists
def intersectionBetweenTwoLists(list1,list2,interscectionFile='F',unique1='F',unique2='F'):
    if type(list1) == str and  type(list2) == str:
        contentList1 = readLinesIntoList(list1)
        contentList2 = readLinesIntoList(list2)
        intersection = set(contentList1).intersection(list2)
        list1Unique = [x for x in contentList1 if x not in intersection]
        list2Unique = [x for x in contentList2 if x not in intersection]
        if intersectionFile != 'F':
            outHandle = open(intersectionFile,'w')
            for iterm in intersection:
                outHandle.write("%s\n"% (str(iterm)))
            outHandle.close()
        if unique1 != 'F' and unique2 != 'F':
            outUnique1 = open(unique1,'w')
            outUnique2 = open(unique2,'w')
            for iterm1,iterm2 in izip_longest(list1Unique,list2Unique,fillvalue=''):
                outUnique1.write(("%s\n"% (str(iterm1))))
                outUnique2.write(("%s\n"% (str(iterm2))))
            outUnique1.close()
            outUnique2.close()
    elif type(list1) == list and type(list2) == list:
        intersection=set(list1).intersection(list2)
        print("Overlap elements number: %f"% (float(len(intersection))))
        print("First list unique elements number: %f" %(float(len(list1) - len(intersection))))
        print("Second list unique elements number: %f " % (float(len(list2) - len(intersection))))
        list1Unique = [x for x in list1 if x not in intersection]
        list2Unique = [x for x in list2 if x not in intersection]
        return {"intersection":intersection,"lsit1Unique":list1Unique,"list2Unique":list2Unique}



##################################################################################
#                       file system
#make directory
def mkdir(folderName):
    import os
    if not os.path.exists(folderName):
        os.makedirs(folderName)
#delete directory
def rmdir(folderName):
    import os
    if os.path.exists(folderName):
        os.system("rm -r " + folderName)

#copy specific files under the same parent folder
def copyThrough(path,pattern,dest="./"):
    import fnmatch, os
    matches = []
    for root, dirnames, filenames in os.walk(path):
        for filename in fnmatch.filter(filenames,pattern):
            matches.append(os.path.join(root,filename))
    for item in matches:
        os.system("cp %s %s" % (item,dest))





