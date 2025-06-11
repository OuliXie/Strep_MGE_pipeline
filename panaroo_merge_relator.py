import csv
import re
import itertools

def parse_gene_pres_abs(pres_abs_file):
    with open(pres_abs_file, 'r', newline='', ) as gene_presence_absence:
            # Read column header line
            gff_file_names = gene_presence_absence.readline()

            # Strip for whitespace
            gff_file_names = gff_file_names.strip()
            # split column names
            gff_file_names = gff_file_names.split(',')

            if 'Genome Fragment' in gff_file_names:
                n_skip_lines = 14
            else:
                n_skip_lines = 5
    
            # Index gff filenames and column position in dict for better search
            gff_file_dict = {}
            for i, gff_name in enumerate(gff_file_names[n_skip_lines:]):
                if 'Phage' == gff_name:
                    break
                
                gff_file_dict[gff_name] = i
    
            # Initialise reader object to read remaining lines
            reader = csv.reader(gene_presence_absence, delimiter=',')

            gene_dict = {item: {} for item in gff_file_dict}
            annotation_dict = {}
            no_iso_dict = {}
    
    		# Read gene presence absence file
            for line in reader:
                if n_skip_lines == 5:
                    # Add classification
                    annotation_dict[line[0]] = line[-1]
                    # Add the number of times a gene is observed
                    no_iso_dict[line[0]] = line[3]
                else:
                    annotation_dict[line[0]] = line[2]

                for genome in gene_dict.keys():
                    # Check if there is an annotation for the given genome
                    if len(line[n_skip_lines + gff_file_dict[genome]]) > 0:
                        gene_dict[genome][line[n_skip_lines+gff_file_dict[genome]]] = line[0]
            
            return gene_dict, annotation_dict, no_iso_dict


def relator(original_dict, merged_dict, return_dict):
    # Go through all genes in a genome and match it to its merged counterpart
    for entry in original_dict:
        # See if the pan gneome cluster is already added
        try:
            return_dict[merged_dict[entry]][original_dict[entry]] += 1
        except KeyError:

            # Try adding it
            try:
                return_dict[merged_dict[entry]][original_dict[entry]] = 1
                
                # if locus_tag is not found try regex searching for it.
            except KeyError:
                keys = [key for key in merged_dict.keys() if re.search(entry+';', key) or re.search(entry+'$', key)]

                # If locus_tag is still not found see if it is a fragment and rearrange into combinations
                if len(keys) == 0 and ';' in entry:
                    perms = list(itertools.permutations(entry.split(';')))
                    perms = [';'.join(iteration) for iteration in perms]

                    keys_dict = {}
                    for name in perms:
                        possible_keys = [key for key in merged_dict.keys() if bool(re.search(name+';', key)) or bool(re.search(name+'$', key)) or name == key]
                        if possible_keys:
                            keys_dict[name] = possible_keys
                    # Extarct and flatten keys found
                    keys = [key for key_list in list(keys_dict.values()) for key in key_list]

                    # If there still is no keys search for each of the locus_tags themselves (likely seperated by a third tag)
                    if len(keys) == 0:
                        keys_dict = {}
                        individ_tags = entry.split(';')
                        for tag in individ_tags:
                            possible_keys = [key for key in merged_dict.keys() if bool(re.search(tag+';', key)) or bool(re.search(tag+'$', key))]
                            keys = list(set(possible_keys))

                # Check that only one key is matched
                if len(keys) != 1:
                    print(keys)
                    print("More or less than one key was matched")
                    print(entry)
                    exit()
                else:
                    try:
                        return_dict[merged_dict[keys[0]]][original_dict[entry]] += 1
                    except KeyError:
                        return_dict[merged_dict[keys[0]]][original_dict[entry]] = 1
    return return_dict


def collector(SDSE_related, SPY_related, merged_pan_clusters, SDSE_classification, SPY_classification, merged_annotation, sdse_no_iso_dict, pyo_no_iso_dict):
    output_related = []

    for cluster in merged_pan_clusters:
        cluster_info = {}

        cluster_info['Merged_cluster'] = cluster

        cluster_info['S_pyo-cluster'] = ';'.join(SPY_related[cluster].keys())

        cluster_info['SDSE-cluster'] = ';'.join(SDSE_related[cluster].keys())

        cluster_info['S_pyo-classification'] = ';'.join([SPY_classification[key] for key in SPY_related[cluster].keys()])

        cluster_info['SDSE-classification'] = ';'.join([SDSE_classification[key] for key in SDSE_related[cluster].keys()])

        cluster_info['S_pyo-No.isolates'] = ';'.join([pyo_no_iso_dict[key] for key in SPY_related[cluster].keys()])

        cluster_info['SDSE-No.isolates'] = ';'.join([sdse_no_iso_dict[key] for key in SDSE_related[cluster].keys()])

        cluster_info['Merged_annotation'] = merged_annotation[cluster]

        output_related.append(cluster_info)

    return output_related


def write_output(rel_dict, out_file_path):
    output_columns = ['Merged_cluster', 'Merged_annotation','SDSE-cluster','S_pyo-cluster', 'SDSE-classification', 'S_pyo-classification', 'SDSE-No.isolates', 'S_pyo-No.isolates']
    with open(out_file_path, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=output_columns)
        writer.writeheader()
        for line in rel_dict:
            writer.writerow(line)
        


if __name__ == '__main__':
    print("Reading SDSE pan")
    sdse_ori_pan, sdse_classification, sdse_no_iso_dict = parse_gene_pres_abs('classified_genes_sdse.csv')
    print("Reading S.pyo pan")
    pyo_ori_pan, pyo_classification, pyo_no_iso_dict = parse_gene_pres_abs('classified_genes_pyo.csv')
    print("Reading merged pan")
    merged_pan, merged_annotations, _ = parse_gene_pres_abs('gene_presence_absence_roary_90.csv')

    # Find all pan-genome clusters in the merged pan-genome
    pan_genome_clusters = set([pan_cluster for genome in merged_pan for pan_cluster in list(merged_pan[genome].values())])
    # Construct dicts to hold the relation to the original genomes
    SDSE_reldict = {key: {} for key in pan_genome_clusters}
    SPY_reldict = {key: {} for key in pan_genome_clusters}

    print("Relating SDSE to merged pan")
    for genome_name in sdse_ori_pan:
       SDSE_reldict = relator(sdse_ori_pan[genome_name], merged_pan[genome_name], SDSE_reldict)
    
    print("Relating Pyo to merged pan")
    for genome_name in pyo_ori_pan:
       SPY_reldict = relator(pyo_ori_pan[genome_name], merged_pan[genome_name], SPY_reldict)

    print("Collecting relations")
    out_rel_dict = collector(SDSE_reldict, SPY_reldict, pan_genome_clusters, sdse_classification, pyo_classification, merged_annotations, sdse_no_iso_dict, pyo_no_iso_dict)

    print("Writing output")
    write_output(out_rel_dict, 'test_out_file.csv')

    print('DONE')


    # pan_genome_clusters = set(['cas', 'pur','uvrA'])
    # SDSE_reldict = {
    #     'cas': {'cas_SDSE': 1},
    #     'pur': {},
    #     'uvrA': {'uvrA_SDSE_1': 1, 'uvrA_SDSE_2': 1}
    # }
    # SPY_reldict = {
    #     'cas': {},
    #     'pur': {'pur_PYO': 1},
    #     'uvrA': {},
    # }

    # sdse_classification = {
    #     'cas': 'non-mge',
    #     'pur': 'phage',
    #     'uvrA': 'ICE'
    # }

    # SPY_classification = {
    #     'cas': 'non-mge',
    #     'pur': 'ICE',
    #     'uvrA': 'phage'
    # }

    # merged_annotation = {
    #     'cas': 'CRISPRcas loci',
    #     'pur': 'purine kinar',
    #     'uvrA': 'uuuuhhhh',
    # }

    # out_rel_dict = collector(SDSE_reldict, SPY_reldict, pan_genome_clusters, sdse_classification, SPY_classification, merged_annotation)
    # for i in out_rel_dict:
    #     print(i)

    # write_output(out_rel_dict, 'test_out_file.csv')
    

    #for key in SDSE_reldict:
    #    if len(SDSE_reldict[key]) > 1:
    #        print(key)
    #        print(SDSE_reldict[key])
    #print(SDSE_reldict)




