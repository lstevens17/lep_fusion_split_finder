#!/usr/bin/env python3
import sys
import argparse

def parse_table(table_file):
	with open(table_file, 'r') as table:
		table_dict, chr2busco_dict = {}, {} # table_dict will record coordinates of each BUSCO; chr2busco_dict will record which BUSCOs are found on which chromosome
		all_chr_list = [] # to keep track of which chromosomes have no BUSCOs
		for line in table:
			if not line.startswith("#"):
				cols = line.rstrip("\n").split()
				buscoID, status = cols[0], cols[1]
				if status == 'Complete': # ignore fragmented or duplicated BUSCOs
					chr, start, stop = cols[2].split(":")[0], int(cols[3]), int(cols[4])
					table_dict[buscoID] = [chr, start, stop] # save coordinates in dict
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
			if not proportion >= min_proportion: # if the proportion of BUSCOs is lower than specified threshold
				non_ancestral_dict[sp1_chr] = [] 
				for sp2_chr in set(sorted(sp2_chr_list)): # for every chromosome in sp2 with these BUSCOs
					proportion = sp2_chr_list.count(sp2_chr)/len(sp2_chr_list) # get the proportion 
					if proportion >= report_proportion: # and if this is higher than the threshold for reporting
						non_ancestral_dict[sp1_chr].append([sp2_chr, proportion, sp2_chr_list.count(sp2_chr)]) # store it in a dict
		except ValueError: # to catch rare cases where none of the BUSCOs on sp1_chr are complete in sp2 (these are usually small chunks of chromosomes)
			print("[+]\tWARNING: None of the BUSCOs from " + sp1_chr + " are complete in other species (BUSCO count: " + str(len(buscoID_list)) + ")")
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
			print("[+]\tWARNING: There are no " + fusion_split_prefix + " chromosomes. '" + prefix + "_" + fusion_split_prefix + "_chromosomes.tsv' is empty.")

def assign_chromosomes(reference_table_dict, query_table_dict, reference_chr2busco_dict, query_chr2busco_dict, prefix):
	print("[+] Assigning unfused/unsplit chromosomes to a reference chromosome...")
	# assign query chr to reference
	query2reference_dict = {}
	for query_chr, buscoID_list in query_chr2busco_dict.items():
		reference_chr_list = []
		for buscoID in buscoID_list:
			try:
				reference_chr_list.append(reference_table_dict[buscoID][0])
			except KeyError:
				pass
		try:
			top_reference_chr = max(set(reference_chr_list), key=reference_chr_list.count)
			proportion = reference_chr_list.count(top_reference_chr)/len(reference_chr_list)
			if proportion >= min_proportion:
				query2reference_dict[query_chr] = top_reference_chr
		except ValueError:
			print("[+]\tWARNING: None of the BUSCOs from " + query_chr + " are complete in other species (BUSCO count: " + str(len(buscoID_list)) + ")")
	# assign reference chr to query
	reference2query_dict = {}
	for reference_chr, buscoID_list in reference_chr2busco_dict.items():
		query_chr_list = []
		for buscoID in buscoID_list:
			try:
				query_chr_list.append(query_table_dict[buscoID][0])
			except KeyError:
				pass
		try:
			top_query_chr = max(set(query_chr_list), key=query_chr_list.count)
			proportion = query_chr_list.count(top_query_chr)/len(query_chr_list)
			if proportion >= min_proportion:
				reference2query_dict[reference_chr] = top_query_chr
		except ValueError:
			print("[+]\tWARNING: None of the BUSCOs from " + reference_chr + " are complete in other species (BUSCO count: " + str(len(buscoID_list)) + ")")
	# check reciprocol 
	reciprocol_list = []
	for query_chr, reference_chr in query2reference_dict.items():
		try:
			if reference2query_dict[reference_chr] == query_chr:
					reciprocol_list.append(query_chr)
			else:
					print("[+]\tWARNING: Cannot assign " + query_chr + " to a reference chromosome (non-reciprocol relationship with top reference chromosome)")
		except KeyError:
			pass # these are the splits and fusions
	return query2reference_dict, reference2query_dict, reciprocol_list

def write_chromosome_assignment_file(reference_table_dict, query_table_dict, reference_chr2busco_dict, query_chr2busco_dict, split_dict, fusion_dict, reciprocol_list):
	# collect info about splits and fusions
	split_list = []
	for ref_chr, component_list in split_dict.items():
		for item in component_list:
			split_list.append(item[0])
	fusion_list = []
	for query_chr, component_list in fusion_dict.items():
		for item in component_list:
			fusion_list.append(item[0])
	# write output file
	with open(prefix + "_chromosome_assignments.tsv", "w") as chromosome_assignment_file:
		chromosome_assignment_file.write(("%s\t%s\t%s\t%s\t%s\t%s\n") % ("query_chr", "status", "assigned_ref_chr", "assigned_ref_BUSCOs", "total_BUSCOs", "prop_BUSCOs"))
		for query_chr, buscoID_list in query_chr2busco_dict.items():
			total_BUSCOs = len(buscoID_list)
			if not query_chr in split_list and not query_chr in fusion_dict and query_chr in reciprocol_list:
				assigned_reference_chr = query2reference_dict[query_chr]
				reference_chr_list = []
				for buscoID in buscoID_list:
					try:
						reference_chr_list.append(reference_table_dict[buscoID][0])
					except KeyError:
						pass
				assigned_reference_BUSCOs = reference_chr_list.count(assigned_reference_chr)
				proportion = assigned_reference_BUSCOs/total_BUSCOs
				chromosome_assignment_file.write(("%s\t%s\t%s\t%s\t%s\t%s\n") % (query_chr, "ancestral", assigned_reference_chr, assigned_reference_BUSCOs, total_BUSCOs, proportion))
			elif query_chr in split_list:
				chromosome_assignment_file.write(("%s\t%s\t%s\t%s\t%s\t%s\n") % (query_chr, "split", "-", "-", total_BUSCOs, "-"))
			elif query_chr in fusion_dict:
				chromosome_assignment_file.write(("%s\t%s\t%s\t%s\t%s\t%s\n") % (query_chr, "fusion", "-", "-", total_BUSCOs, "-"))
		print("[+]\tSuccessfully written chromosome assignments to '" + prefix + "_chromosome_assignments.tsv'")


if __name__ == "__main__":
	SCRIPT = "fusion_split_finder.py"
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
	query2reference_dict, reference2query_dict, reciprocol_list = assign_chromosomes(reference_table_dict, query_table_dict, reference_chr2busco_dict, query_chr2busco_dict, prefix)
	write_chromosome_assignment_file(reference_table_dict, query_table_dict, reference_chr2busco_dict, query_chr2busco_dict, split_dict, fusion_dict, reciprocol_list)
