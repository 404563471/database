#!/usr/bin/python
import os
def loadGeneInfo(geneinfo):
	geneinfo_data={}
	geneinfo_file=open(geneinfo, 'r')
	geneinfo_file.readline()
	while True:
		line_geneinfo = geneinfo_file.readline()
		if line_geneinfo:
			temp_geneinfo=[]
			temp_geneinfo= line_geneinfo.strip("\n").replace('"','').split(",")
			gene_symbol=temp_geneinfo[2]
			gene_id=temp_geneinfo[1]
			gene_descrip=temp_geneinfo[8]
			if not gene_symbol in geneinfo_data:
				geneinfo_data[gene_symbol] = gene_id+","+gene_descrip
		else:
			break
	print len(geneinfo_data)
	#print edges
	geneinfo_file.close()
	return geneinfo_data

def __main__():
	geneinfo="Homo_sapiens.gene_info.txt"
	os.mkdir("genes_edges_json")
	geneinfo_data={}
	geneinfo_data=loadGeneInfo(geneinfo)
	
	for filename in os.listdir("genes_edges"):
		genes_edges_file=open("genes_edges/"+filename, 'r')
		query_gene=filename.strip("_genes_edges.txt")
		genes_edges_json_file=open("genes_edges_json/"+query_gene+"_ppi.json",'w')
		temp_info=geneinfo_data[query_gene].split(",")
		genes_edges_json_file.write('{"genes": [{"query": true, "sys_name": "'+query_gene+'","std_name": "'+query_gene+'","id":'+temp_info[0]+',"descrip":"'+temp_info[1]+'"}')
		genes_edges_json_file.flush()
		genes={}
		edges={}
		genes[query_gene]=0
		genes_edges_file.readline()
		i=1
		while True:
			line = genes_edges_file.readline()
			if line:
				temp=[]
				temp= line.strip("\n").replace('"','').split(",")
				if not temp[1] in genes:
					genes[temp[1]] = i
					i=i+1
				if not temp[2] in genes:
					genes[temp[2]] = i
					i=i+1
				source_temp=str(genes[temp[1]])
				target_temp=str(genes[temp[2]])
				edge_temp=source_temp+","+target_temp
				if not edge_temp in edges:
					if temp[3] == "NA":
						edges[edge_temp]=0.0
					else:
						edges[edge_temp]=float(temp[3])
			else:
				break
		print len(edges)
		genes_edges_file.close()
		for gene in genes.keys():
			if gene != query_gene:
				temp_info=geneinfo_data[gene].split(",")
				genes_edges_json_file.write(',{"query": false, "sys_name": "'+gene+'","std_name": "'+gene+'","id":'+temp_info[0]+',"descrip":"'+temp_info[1]+'"}')
				genes_edges_json_file.flush()
		genes_edges_json_file.write('],"edges": [')
		i=0
		#print edges[0]
		#sorted_edges=sorted(edges.iteritems(),key=lambda d:d[1],reverse=True)
		for edge in sorted(edges.items(), lambda x, y: cmp(x[1], y[1]), reverse=True):
			source_temp=edge[0].split(",")[0]
			target_temp=edge[0].split(",")[1]
			if i<1274:
				genes_edges_json_file.write('{"source":'+source_temp+',"target":'+target_temp+',"weight":'+str(edge[1])+',"id":'+str(i)+'},')
				genes_edges_json_file.flush()
			else:
				genes_edges_json_file.write('{"source":'+source_temp+',"target":'+target_temp+',"weight":'+str(edge[1])+',"id":'+str(i)+'}]}')
				genes_edges_json_file.flush()
			i=i+1
		genes_edges_json_file.close()

if __name__=="__main__": __main__()

