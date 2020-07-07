import os
import shutil
import re

def change_fasta_name(data_path):
	'''用于更改fasta文件名称为ncbi编号'''
	filelist=os.listdir(data_path)
	new_filelist=[]
	for i in filelist:
		old_name = data_path + i
		with open(old_name, "r") as fasta:
			fasta_name = fasta.readline()
		number = re.search("\w{2}\d+", fasta_name)
		if number is None:
			print(fasta_name, "\nfilenameis: ", i)
		else:
			new_name=data_path+number.group()
			os.rename(old_name, new_name)
			new_filelist.append(new_name)

	return new_filelist

change_fasta_name(data_path="pydata/NCBI/")

data_path="pydata/lab/"

def fasta_to_file(data_path, dir_path, pattern):
	"""用于把多序列fasta文件分裂成单序列fasta文件"""
	print(data_path)
	with open(data_path, "r") as fastafile:
		total="".join(fastafile.readlines()).split(">")
	#print(total)
	for i in total[1:]:
		i=">"+i
		#用于去除序列中的\n
		i = i.split("\n")
		i = "\n".join([i[0], "".join(i[1:])])
		filename=re.search(pattern, i)
		if filename is None:
			print(i, filename)
		else:
			singlefastaname=dir_path+filename.group(1)
			print(singlefastaname)
			with open(singlefastaname, "w") as  singlefasta:
				singlefasta.write(i)

def total_fasta_to_file(data_dir_path, dir_path):
	if os.path.isdir(dir_path):
		#注意当删除非空文件夹时使用shutil.retree
		shutil.rmtree(dir_path)
		os.makedirs(dir_path)
	else:
		os.makedirs(dir_path)
	fasta_list=os.listdir(data_dir_path)
	#print(fasta_list)
	for i in fasta_list:
		data_path=data_dir_path+i
#		"""根据文件类型选择正则表达式"""
		#if not i=='BOLDSYSTEM.fasta':
		#	fasta_to_file(data_path=data_path, dir_path=dir_path, pattern=">(.+)\s[A-Z][a-z]")
		#else:
			#fasta_to_file(data_path=data_path, dir_path=dir_path, pattern=">(.+\d\d)")
		fasta_to_file(data_path=data_path, dir_path=dir_path, pattern=">(.+?)\s")
#total_fasta_to_file("pydata/lab/", "pydata/lab-single/")

#total_fasta_to_file("trans/03-06/", "trans/03-06-single/")

total_fasta_to_file("trans/01/", "trans/01-single/")
total_fasta_to_file("trans/02/", "trans/02-single/")