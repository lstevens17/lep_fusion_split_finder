#!/usr/bin/env python3
import sys
import argparse

def parse_table(table_file):
	with open(table_file, 'r') as table:
		table_dict, chr2busco_dict = {}, {}
		all_chr_list = []
		for line in table:
			if not line.startswith("#"):
				cols = line.rstrip("\n").split()
				buscoID, status = cols[0], cols[1]
				if status == 'Complete':
					chr, start, stop = cols[2], int(cols[3]), int(cols[4])
					table_dict[buscoID] = [chr, start, stop]
					try:
						chr2busco_dict[chr].append(buscoID)
					except KeyError:
						chr2busco_dict[chr] = [buscoID]
					if not chr in all_chr_list:
						all_chr_list.append(chr)
				elif status == 'Duplicated':
					chr = cols[2]
					if not chr in all_chr_list:
						all_chr_list.append(chr)
		for chr in all_chr_list:
			if not chr in chr2busco_dict:
				print("[+] WARNING: " + chr + " has no complete BUSCOs and cannot be assigned to a chromosome. It will be omitted from the 'chrosome_assignments.tsv' file.")
	return table_dict, chr2busco_dict

def find_fusions_and_splits(sp1_table_dict, sp2_table_dict, sp1_chr2busco_dict, sp2_chr2busco_dict, min_proportion):
	non_ancestral_dict = {}
	for sp1_chr, buscoID_list in sp1_chr2busco_dict.items():
		sp2_chr_list = []
		for buscoID in buscoID_list:
			try:
				sp2_chr_list.append(sp2_table_dict[buscoID][0])
			except KeyError:
				pass
		top_sp2_chr = max(set(sp2_chr_list), key=sp2_chr_list.count)
		proportion = sp2_chr_list.count(top_sp2_chr)/len(sp2_chr_list)
		if not proportion > min_proportion:
			non_ancestral_dict[sp1_chr] = []
			for sp2_chr in set(sorted(sp2_chr_list)):
				proportion = sp2_chr_list.count(sp2_chr)/len(sp2_chr_list)
				if proportion > 0.05:
					non_ancestral_dict[sp1_chr].append([sp2_chr, proportion])
	return(non_ancestral_dict)

def write_fusions_and_splits_files(fusion_split_dict, prefix):
	with open(prefix + "_chromosomes.tsv", 'w') as fusion_split_file:
		if len(fusion_split_dict) != 0:
			for query_chr, list in fusion_split_dict.items(): 
				outlist = []
				for ancestral_chr_prop in list:
					ancestral_chr, prop = ancestral_chr_prop[0], round(ancestral_chr_prop[1], 2)
					outlist.append(ancestral_chr + ":" + str(prop))
				fusion_split_file.write(("%s\t%s\n") % (query_chr, ",".join(outlist)))
			print("[+]\tSuccessfully written " + prefix + "_chromosomes.tsv")
		else:
			print("[+]\tWARNING: " + prefix + "_chromosomes.tsv is empty!")


def assign_chromosomes(reference_table_dict, query_table_dict, reference_chr2busco_dict, query_chr2busco_dict, split_dict, fusion_dict):
	with open("chromosome_assignments.tsv", "w") as chromosome_assignment_file:
		print("[+] Assigning unfused/unsplit chromosomes to a reference chromosome...")
		chromosome_assignment_file.write(("%s\t%s\t%s\t%s\t%s\t%s\n") % ("query_chr", "status", "assigned_ref_chr", "assigned_ref_BUSCOs", "total_BUSCOs", "prop_BUSCOs"))
		for query_chr, buscoID_list in query_chr2busco_dict.items():
			if not query_chr in split_dict and not query_chr in fusion_dict:
				reference_chr_list = []
				for buscoID in buscoID_list:
					try:
						reference_chr_list.append(reference_table_dict[buscoID][0])
					except KeyError:
						pass
				top_reference_chr = max(set(reference_chr_list), key=reference_chr_list.count)
				assigned_ref_BUSCOs, total_BUSCOs = reference_chr_list.count(top_reference_chr), len(reference_chr_list)
				proportion = assigned_ref_BUSCOs/total_BUSCOs
				chromosome_assignment_file.write(("%s\t%s\t%s\t%s\t%s\t%s\n") % (query_chr, "ancestral", top_reference_chr, assigned_ref_BUSCOs, total_BUSCOs, proportion))
			elif query_chr in split_dict:
				chromosome_assignment_file.write(("%s\t%s\t%s\t%s\t%s\t%s\n") % (query_chr, "split", "-", "-", "-", "-"))
			elif query_chr in fusion_dict:
				chromosome_assignment_file.write(("%s\t%s\t%s\t%s\t%s\t%s\n") % (query_chr, "fusion", "-", "-", "-", "-"))
		print("[+]\tSuccessfully written " + str(len(query_chr2busco_dict)) + " chromosomes to chromosome_assignments.tsv")

if __name__ == "__main__":
	SCRIPT = "buscopainter.py"
	# argument set up
	parser = argparse.ArgumentParser()
	parser.add_argument("-r", "--reference_table", type=str, help = "full_table.tsv file for reference species", required=True)
	parser.add_argument("-q", "--query_table", type=str, help = "full_table.tsv for query species", required=True)
	parser.add_argument("-m", "--min_proportion", type=float, help = "Minimum proportion of BUSCO genes used to infer fusions/splits", default=0.9)
	args = parser.parse_args()
	reference_table_file = args.reference_table
	query_table_file = args.query_table
	min_proportion = args.min_proportion
	# run the functions
	reference_table_dict, reference_chr2busco_dict = parse_table(reference_table_file)
	query_table_dict, query_chr2busco_dict = parse_table(query_table_file)
	print("[+] Identifying fused and split chromosomes...")
	split_dict = find_fusions_and_splits(reference_table_dict, query_table_dict, reference_chr2busco_dict, query_chr2busco_dict, min_proportion)
	print("[+]\tIdentified " + str(len(split_dict)) + " split chromosomes")
	write_fusions_and_splits_files(split_dict, "split")
	fusion_dict = find_fusions_and_splits(query_table_dict, reference_table_dict, query_chr2busco_dict, reference_chr2busco_dict, min_proportion)
	print("[+]\tIdentified " + str(len(fusion_dict)) + " fused chromosomes")
	write_fusions_and_splits_files(fusion_dict, "fused")
	assign_chromosomes(reference_table_dict, query_table_dict, reference_chr2busco_dict, query_chr2busco_dict, split_dict, fusion_dict)