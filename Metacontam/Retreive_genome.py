import os
import random
from collections import defaultdict
import copy
from Bio.Seq import Seq
from Bio import SeqIO

class TaxBag:
    def __init__(self):
        self.bag_dict = defaultdict(list)
    def add_tax(self,taxid,ptaxid):
        self.bag_dict[ptaxid]+=[taxid]
    def find_all_tax(self,input_tax):
        self.new_bag=copy.deepcopy(self.bag_dict)
        removed_tax=[]
        determinant = False
        while not determinant:
            for tax_element in self.new_bag[input_tax]:
                if tax_element in self.new_bag.keys():
                    #print(tax_element)
                    self.new_bag[input_tax].remove(tax_element)
                    removed_tax.append(tax_element)
                    self.new_bag[input_tax]+=self.bag_dict[tax_element]
            #print([item not in self.bag_dict.keys() for item in new_bag[input_tax]],new_bag[input_tax])
            determinant=all(item not in self.bag_dict.keys() for item in self.new_bag[input_tax])
        removed_tax=list(set(removed_tax))
        self.new_bag[input_tax]+=removed_tax
        
        return self.new_bag[input_tax]


def find_contaminant_genome(KRAKEN_DB, contaminant_community_taxid, output, seed=5):
    random.seed(seed)
    
    contaminant_community_taxid = [str(i) for i in contaminant_community_taxid]
    TaxBag_Instance = TaxBag()
    
    with open(os.path.join(KRAKEN_DB, "taxonomy", "nodes.dmp"), "r") as file:
        for line in file:
            taxid = line.split("|")[0].strip()
            ptaxid = line.split("|")[1].strip()
            rank = line.split("|")[2].strip()
            TaxBag_Instance.add_tax(taxid, ptaxid)

    # Find species taxid that do not exist in Library fasta file
    taxid_exist = []
    with open(os.path.join(KRAKEN_DB, "seqid2taxid.map"), "r") as mapping:
        for line in mapping:
            try:
                info = line.split("\t")[0]
                taxid = info.split("|")[1].strip()
                genome = info.split("|")[2]
                taxid_exist.append(taxid)
            except:
                pass

    not_exist_list = set(contaminant_community_taxid) - set(taxid_exist)

    not_exist_dict = {}
    for taxid in sorted(not_exist_list):  # Sort to ensure consistent order
        valid_taxids = sorted(set(taxid_exist).intersection(TaxBag_Instance.find_all_tax(taxid)))
        if valid_taxids:
            not_exist_dict[taxid] = random.choice(valid_taxids)
        else:
            not_exist_dict[taxid] = None
    

    with open(os.path.join(output,"Network_Output","Contaminant_community_taxid"),"w") as file:
        for taxid in contaminant_community_taxid:
            file.writelines([taxid,"\n"])


    if len(not_exist_list) > 0:
        print(f"{len(not_exist_list)} species taxid do not exist in kraken2 library files")


    for speciestax, straintax in not_exist_dict.items():
        print(f"Replace species taxid {speciestax} to strain taxid {straintax}")
    
    contaminant_community_taxid = [not_exist_dict.get(taxid, taxid) for taxid in contaminant_community_taxid]
    
    return not_exist_dict, contaminant_community_taxid


    




#####################################################################



def retreive_genome_from_library(KRAKEN_DB, output_dir, not_exist_dict, contaminant_community_taxid):
    bac_sequences = os.path.join(KRAKEN_DB, "library", "bacteria", "library.fna")
    vir_sequences = os.path.join(KRAKEN_DB, "library", "viral",    "library.fna")
    arc_sequences = os.path.join(KRAKEN_DB, "library", "archaea",  "library.fna")
    sequence_dbs = [bac_sequences, vir_sequences, arc_sequences]
    copyed_contam_candidate = copy.deepcopy(contaminant_community_taxid)
    genome_path = os.path.join(output_dir,"Genome_dir", "Candidate.fasta")
    not_exist_dict_reverse={l:i for i,l in not_exist_dict.items()}


    with open(genome_path, "w") as filtered:
        for sequence_db in sequence_dbs:
            origin = sequence_db
            records = SeqIO.parse(origin, "fasta")
            for record in  records:
                tax_id = record.id.split("|")[1]
                if (tax_id in copyed_contam_candidate) and ("plasmid" not in record.description) :
                    if tax_id in not_exist_dict_reverse:
                        strain_taxid=tax_id
                        species_taxid=not_exist_dict_reverse[tax_id]
                        ID=str(record.id)
                        record.id=ID.replace(f'taxid|{strain_taxid}',f'taxid|{species_taxid}')
                    copyed_contam_candidate.remove(tax_id)
                    # Convert Seq to string, replace 'x' with 'N'
                    new_sequence = str(record.seq).replace("x", "N")
                    # Convert back to Seq object
                    record.seq = Seq(new_sequence)
                    # Write record to output file
                    SeqIO.write(record, filtered, "fasta")
        return copyed_contam_candidate
    
    
#####################################################################
