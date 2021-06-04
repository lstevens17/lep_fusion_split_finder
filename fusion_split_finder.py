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
					chr, start, stop = cols[2].split(":")[0], int(cols[3]), int(cols[4])
					table_dict[buscoID] = [chr, start, stop]
					try:
						chr2busco_dict[chr].append(buscoID)
					except KeyError:
						chr2busco_dict[chr] = [buscoID]
					if not chr in all_chr_list:
						all_chr_list.append(chr)
				elif status == 'Duplicated':
					chr = cols[2].split(":")[0]
					if not chr in all_chr_list:
						all_chr_list.append(chr)
		for chr in all_chr_list:
			if not chr in chr2busco_dict:
				print("[+]\tWARNING: " + chr + " has no complete BUSCOs and cannot be assigned to a chromosome. It will be omitted from all output files!")
	return table_dict, chr2busco_dict

def find_fusions_and_splits(sp1_table_dict, sp2_table_dict, sp1_chr2busco_dict, sp2_chr2busco_dict, min_proportion, report_proportion):
	non_ancestral_dict = {}
	for sp1_chr, buscoID_list in sp1_chr2busco_dict.items():
		sp2_chr_list = []
		for buscoID in buscoID_list:
			try:
				sp2_chr_list.append(sp2_table_dict[buscoID][0])
			except KeyError:
				pass
		try:
			top_sp2_chr = max(set(sp2_chr_list), key=sp2_chr_list.count)
			proportion = sp2_chr_list.count(top_sp2_chr)/len(sp2_chr_list)
			if not proportion > min_proportion:
				non_ancestral_dict[sp1_chr] = []
				for sp2_chr in set(sorted(sp2_chr_list)):
					proportion = sp2_chr_list.count(sp2_chr)/len(sp2_chr_list)
					if proportion >= report_proportion:
						non_ancestral_dict[sp1_chr].append([sp2_chr, proportion, sp2_chr_list.count(sp2_chr)])
		except ValueError: # for the rare case that none of the genes on an ancestral chromosome are found ANYWHERE
			# WTF why does this happen???
			print("[+]\tWARNING: All BUSCOs from " + sp1_chr + " are missing in the other species!")
	return(non_ancestral_dict)

def write_fusions_and_splits_files(fusion_split_dict, fusion_split_prefix, prefix):
	with open(prefix + "_" + fusion_split_prefix + "_chromosomes.tsv", 'w') as fusion_split_file:
		if len(fusion_split_dict) != 0:
			for query_chr, list in fusion_split_dict.items(): 
				outlist = []
				for ancestral_chr_prop in list:
					ancestral_chr, prop, count = ancestral_chr_prop[0], round(ancestral_chr_prop[1], 2), ancestral_chr_prop[2]
					outlist.append(ancestral_chr + ":" + str(count) + ":" + str(prop))
				fusion_split_file.write(("%s\t%s\n") % (query_chr, ",".join(outlist)))
			print("[+]\tSuccessfully written '" + prefix + "_" + fusion_split_prefix + "_chromosomes.tsv'")
		else:
			print("[+]\tWARNING: '" + prefix + "_" + fusion_split_prefix + "_chromosomes.tsv' is empty!")


def assign_chromosomes(reference_table_dict, query_table_dict, reference_chr2busco_dict, query_chr2busco_dict, split_dict, fusion_dict, prefix):
	with open(prefix + "_chromosome_assignments.tsv", "w") as chromosome_assignment_file:
		print("[+] Assigning unfused/unsplit chromosomes to a reference chromosome...")
		chromosome_assignment_file.write(("%s\t%s\t%s\t%s\t%s\t%s\n") % ("query_chr", "status", "assigned_ref_chr", "assigned_ref_BUSCOs", "total_BUSCOs", "prop_BUSCOs"))
		### ALLOW SEARCHING OF SPLITS
		split_list = []
		for ref_chr, component_list in split_dict.items():
			for item in component_list:
				split_list.append(item[0])
		## ALLOW SEARCHING OF FUSIONS
		fusion_list = []
		for query_chr, component_list in fusion_dict.items():
			for item in component_list:
				fusion_list.append(item[0])
		# ASSIGN REF TO QUERY
		reference2query_dict = {}
		for reference_chr, buscoID_list in reference_chr2busco_dict.items():
			total_BUSCOs = len(buscoID_list)
			if not reference_chr in split_dict and not reference_chr in fusion_list:
				query_chr_list = []
				for buscoID in buscoID_list:
					try:
						query_chr_list.append(query_table_dict[buscoID][0])
					except KeyError:
						pass
				try:
					top_query_chr = max(set(query_chr_list), key=query_chr_list.count)
					reference2query_dict[reference_chr] = top_query_chr
				except ValueError:
					print("[+]\tWARNING: All BUSCOs from " + reference_chr + " are missing in the other species!")
		# ASSIGN QUERY TO REF
		for query_chr, buscoID_list in query_chr2busco_dict.items():
			total_BUSCOs = len(buscoID_list) # note to self: this means the proprotion of assigned BUSCOs is relative to ALL BUSCOs, not just those that are also found in ref species
			if not query_chr in split_list and not query_chr in fusion_dict:
				reference_chr_list = []
				for buscoID in buscoID_list:
					try:
						reference_chr_list.append(reference_table_dict[buscoID][0])
					except KeyError:
						pass
				try:
					top_reference_chr = max(set(reference_chr_list), key=reference_chr_list.count) # NEED TO ADD RECIPRICOL!
					assigned_ref_BUSCOs = reference_chr_list.count(top_reference_chr)
					proportion = assigned_ref_BUSCOs/total_BUSCOs
					try:
						if query_chr == reference2query_dict[top_reference_chr]:
							chromosome_assignment_file.write(("%s\t%s\t%s\t%s\t%s\t%s\n") % (query_chr, "ancestral", top_reference_chr, assigned_ref_BUSCOs, total_BUSCOs, proportion))
						else:
							print("[+]\tWARNING: Cannot assign " + query_chr + " to a reference chromosome")
							chromosome_assignment_file.write(("%s\t%s\t%s\t%s\t%s\t%s\n") % (query_chr, "unassigned", "-", "-", total_BUSCOs, "-"))
					except KeyError:
						print("[+]\tWARNING: Cannot assign " + query_chr + " to a reference chromosome")
						chromosome_assignment_file.write(("%s\t%s\t%s\t%s\t%s\t%s\n") % (query_chr, "unassigned", "-", "-", total_BUSCOs, "-"))
				except ValueError: 
					print("[+]\tWARNING: Cannot assign " + query_chr + " to reference chromosome - all BUSCOs are missing in other species!")
			elif query_chr in split_list:
				chromosome_assignment_file.write(("%s\t%s\t%s\t%s\t%s\t%s\n") % (query_chr, "split", "-", "-", total_BUSCOs, "-"))
			elif query_chr in fusion_dict:
				chromosome_assignment_file.write(("%s\t%s\t%s\t%s\t%s\t%s\n") % (query_chr, "fusion", "-", "-", total_BUSCOs, "-"))
		print("[+]\tSuccessfully written " + str(len(query_chr2busco_dict)) + " chromosomes to '" + prefix + "_chromosome_assignments.tsv'")

if __name__ == "__main__":
	SCRIPT = "buscopainter.py"
	# argument set up
	parser = argparse.ArgumentParser()
	parser.add_argument("-r", "--reference_table", type=str, help = "full_table.tsv file for reference species", required=True)
	parser.add_argument("-q", "--query_table", type=str, help = "full_table.tsv for query species", required=True)
	parser.add_argument("-p", "--min_proportion", type=float, help = "Minimum proportion of BUSCO genes used to identify ancestral chromosomes", default=0.9)
	parser.add_argument("-m", "--report_proportion", type=float, help = "Minimum proportion of BUSCOs required report for each fused/split chromosome", default=0.05)
	parser.add_argument("-f", "--prefix", type=str, help = "Prefix for all output files", default="fsf")
	args = parser.parse_args()
	reference_table_file = args.reference_table
	query_table_file = args.query_table
	min_proportion = args.min_proportion
	report_proportion = args.report_proportion
	prefix = args.prefix
	# run the functions
	print("[+] Parsing input files...")
	reference_table_dict, reference_chr2busco_dict = parse_table(reference_table_file)
	query_table_dict, query_chr2busco_dict = parse_table(query_table_file)
	print("[+] Identifying fused and split chromosomes...")
	split_dict = find_fusions_and_splits(reference_table_dict, query_table_dict, reference_chr2busco_dict, query_chr2busco_dict, min_proportion, report_proportion)
	print("[+]\tIdentified " + str(len(split_dict)) + " split chromosomes")
	write_fusions_and_splits_files(split_dict, "split", prefix)
	fusion_dict = find_fusions_and_splits(query_table_dict, reference_table_dict, query_chr2busco_dict, reference_chr2busco_dict, min_proportion, report_proportion)
	print("[+]\tIdentified " + str(len(fusion_dict)) + " fused chromosomes")
	write_fusions_and_splits_files(fusion_dict, "fused", prefix)
	assign_chromosomes(reference_table_dict, query_table_dict, reference_chr2busco_dict, query_chr2busco_dict, split_dict, fusion_dict, prefix)