#_*_ coding:utf-8 _*_

# 该脚本用于整理物种的fasta,gff文件,同时生成jbrowse数据脚本和配置文件

import os
import shutil

# rootdir为实际存放物种数据的文件夹,每次运行需重置该文件夹,切记备份!
rootdir="/home/yhy/Documents/my-project/生物学数据库/白蚁数据库/20190123客户线粒体资料/线粒体基因组可视化"
files=os.listdir(rootdir)

# dockerdir_rootpath为jbrowse容器内部存放整理好的数据的文件夹
dockerdir_rootpath="/jbrowse/jbrowse-data"
names=[]
jbrowse_sh=open("jbrowse-baiyi.sh", "w") #批量制作参考基因组和注释序列的脚本,需在jbrowse容器内运行
jbrowse_config=open("jbrowse.conf", "w") #jbrowse.conf是配置文件, 后期把内容cat到真正配置文件的后面就可以


for	i in files :
	filepath=rootdir + "/" + i
	dirname=os.path.splitext(i)[0]
	dirname=dirname.replace(" ", "-")
	renamefile=filepath.replace(" ", "-")
	os.rename(filepath, renamefile)
	if os.path.isfile(renamefile):
		dirpath=os.path.splitext(renamefile)[0]
		extname=os.path.splitext(renamefile)[1]
		if not dirpath in names :
			names.append(dirpath)
			os.makedirs(dirpath)
		shutil.move(renamefile, dirpath)
		#编写jbrowse脚本
		dockerfilepath=dockerdir_rootpath + "/" + dirname + "/" + i.replace(" ", "-")
		dockerdirpath=os.path.split(dockerfilepath)[0]
		if extname == ".fasta" :
			jbrowse_sh.write("prepare-refseqs.pl --fasta %s --out %s \n" %(dockerfilepath, dockerdirpath))
		else :
			jbrowse_sh.write("flatfile-to-json.pl --gff %s --trackLabel rna --out %s \n" %(dockerfilepath, dockerdirpath) )
			urlpath="/".join(dockerdirpath.split("/")[2:])
			jbrowse_config.write("[ %s ]\nurl = ?data=%s\nname = %s\n" %(dirname, urlpath, dirname))

jbrowse_sh.close()
jbrowse_config.close()


