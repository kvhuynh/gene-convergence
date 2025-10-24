from Bio import SeqIO;
import json
import pprint;

# Genes must be directly adjacent for a convergent pair to be counted
# Fixed sliding window of length 2
# 1 == positive strand, -1 == negative strand

def read_record(genome_record):
    genes = [];
    count = 0;
    for feature in genome_record.features:
        if feature.type == "gene":
            genes.append(feature);
            count += 1
    print(count);
    return genes;
                
def find_convergence(genes):
    res = {}
    left = 0;
    right = 0;
    while right < len(genes):
        strand_left = genes[left].location.strand;
        location_end_left = genes[left].location.end;
        strand_right = genes[right].location.strand;
        location_start_right = genes[right].location.start;
        if strand_left == 1 and strand_right == -1:
            # replace locus_tag with "gene" if key error (APMV, ASFV)
            gene_name_left = genes[left].qualifiers["locus_tag"][0];
            gene_name_right = genes[right].qualifiers["locus_tag"][0];
            res[gene_name_left] = {
                "right": gene_name_right,
                "left_end_position": int(location_end_left),
                "right_start_position": int(location_start_right),           "converging_distance": location_start_right - location_end_left           
            }
        right += 1;
        left = right - 1
    pprint.pprint(res)
    return res;

path = "./input/ASFV.gb";
genome_record = SeqIO.read(path, "genbank");
genes = read_record(genome_record);
res = find_convergence(genes);
with open("./output/ASFV_convergence.json", 'w') as f:
    json.dump(res, f, indent=4)