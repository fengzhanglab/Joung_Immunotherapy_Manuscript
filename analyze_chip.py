import csv
import numpy as np

refseq_file = "hg38_refseq.csv"
genes = ["JUNB_352", "JUNB_353"]
chip_files = [("ChIP_"+g+"_peaks_fixed") for g in genes]
chrom_list = ['chr' + str(x+1) for x in range(22)] + ['chrX', 'chrY']
upstream = 2000
downstream = 500

#read refseq annotations and return the promoter regions
def read_annotation(filename):
	gene_tss = {}
	with open(filename, 'rb') as gf:
		f = [row for row in csv.reader(gf.read().splitlines())]
		for i,l in enumerate(f): # i is index, l is entry
			if i == 0: 
				columns = l
				continue

			#fetch the current transcript
			transcript = dict([(columns[i],e) for i,e in enumerate(l)])
			chrom = transcript["chrom"]
			gene = transcript["name2"]
			if chrom not in chrom_list:
				continue
			if transcript["strand"] == "+":
				tss = [chrom, int(transcript["txStart"])-upstream, int(transcript["txStart"])+downstream]
			else:
				tss = [chrom, int(transcript["txEnd"])-downstream, int(transcript["txEnd"])+upstream]
			if gene in gene_tss.keys():
				gene_tss[gene].append(tss)
			else:
				gene_tss[gene] = [tss]
	return gene_tss

#determine whether genomic regions defined by window1 and window2 overlap
def window_overlap(window1, window2):
	if window1[0] == window2[0]:
		if window2[1] < window1[1] and window2[2] > window1[1]:
			return True
		elif window2[1] > window1[1] and window2[2] < window1[2]:
			return True
		elif window2[1] < window1[2] and window2[2] > window1[2]:
			return True
		else:
			return False
	else:
		return False

#determine whether the peak summit is contained within the genomic region defined by window
def summit_overlap(summit, window):
	if summit[0] == window[0]:
		if summit[1] > window[1] and summit[1] < window[2]:
			return True
		else:
			return False
	else:
		return False

#read N- and C-term flag tag ChIP peaks and output peaks that overlapped
with open(chip_files[0]+'.csv', 'rb') as gf:
	data1 = [row for row in csv.reader(gf.read().splitlines())]
with open(chip_files[1]+'.csv', 'rb') as gf:
	data2 = [row for row in csv.reader(gf.read().splitlines())]
with open("ChIP_JUNB_overlap_peaks.csv", 'w') as csvfile:
	csvwriter = csv.writer(csvfile)
	csvwriter.writerow(data1[0])
	for row1 in data1[1:]:
		if row1[0] not in chrom_list:
			continue
		peak1 = row1[:3]
		for row2 in data2[1:]:
			peak2 = row2[:3]
			if window_overlap(peak1, peak2):
				csvwriter.writerow(row1)
				break

#map overlapping ChIP peaks to respective genes
gene_tss = read_annotation(refseq_file)
with open("ChIP_JUNB_overlap_peaks.csv", 'rb') as gf:
	f = [row for row in csv.reader(gf.read().splitlines())]
	annotations = []
	for i,l in enumerate(f): # i is index, l is entry
		if i == 0: 
			columns = l
			continue
		peak = dict([(columns[i],e) for i,e in enumerate(l)])
		if peak["chr"] not in chrom_list:
			continue
		summit = [peak["chr"], int(peak["start"]) + int(peak["summit"])]
		for g in gene_tss.keys():
			for w in gene_tss[g]:
				if summit_overlap(summit, w):
					if g not in annotations:
						annotations.append(g)
	with open("ChIP_JUNB_overlap_genes_2000_500.csv", 'w') as csvfile:
		csvwriter = csv.writer(csvfile)
		for a in annotations:
			csvwriter.writerow([a])

#output the genes that overlapped between ChIP- and RNA-seq
with open("ChIP_JUNB_overlap_genes_2000_500.csv", 'rb') as gf:
	chip_genes = [row[0] for row in csv.reader(gf.read().splitlines())]
with open("../RNAseq/190113/X342_0.01_exp10.csv", 'rb') as gf:
	rna_data = [row for row in csv.reader(gf.read().splitlines())]
	rna_genes = [row[0] for row in rna_data[1:]]
overlap = list(set(chip_genes).intersection(set(rna_genes)))
with open("ChIP_JUNB_rna_overlap.csv", 'w') as csvfile:
	csvwriter = csv.writer(csvfile)
	csvwriter.writerow(rna_data[0])
	for row in rna_data[1:]:
		if row[0] in overlap:
			csvwriter.writerow(row)
