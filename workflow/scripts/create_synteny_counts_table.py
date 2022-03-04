import pickle
import statistics

orthogroups_2_species_2_genes = pickle.load(open(snakemake.input.orthogroups_2_species_2_genes, 'rb'))
ordered_gff_genes = snakemake.input.ordered_gff_genes
output_synteny_dict = open(snakemake.output.synteny_dict, 'wb')
output_synteny_counts = open(snakemake.output.synteny_counts, 'w')

ordered_gff_species_2_orthogroups = {}
ordered_gff_orthogroups_2_species = {}

with open(ordered_gff_genes) as f:
    next(f)
    for line in f:
        splitline = line.strip().split('\t')
        spec = splitline[0]
        orthogroup = splitline[1]

        if spec not in ordered_gff_species_2_orthogroups.keys():
            ordered_gff_species_2_orthogroups[spec] = [orthogroup]
        # Avoiding duplicates
        elif orthogroup not in ordered_gff_species_2_orthogroups[spec]:
            ordered_gff_species_2_orthogroups[spec].append(orthogroup)

        if orthogroup not in ordered_gff_orthogroups_2_species.keys():
            ordered_gff_orthogroups_2_species[orthogroup] = [spec]
        # Avoiding duplicates
        elif spec not in ordered_gff_orthogroups_2_species[orthogroup]:
            ordered_gff_orthogroups_2_species[orthogroup].append(spec)


def get_left_orthogroup(spec, index):
    if index <= len(ordered_gff_species_2_orthogroups[spec]) and index > 0:
        return ordered_gff_species_2_orthogroups[spec][index-1]
    else:
        return None

def get_right_orthogroup(spec, index):
    if index <= len(ordered_gff_species_2_orthogroups[spec])-2 and index > -1:
        return ordered_gff_species_2_orthogroups[spec][index+1]
    else:
        return None


def get_closest_different_left_orthogroup(spec, index, left_orthogroup):
    left_orthogroup_index = index-1
    closest_different_left_orthogroup = get_left_orthogroup(spec, left_orthogroup_index)
    while closest_different_left_orthogroup == left_orthogroup:
        left_orthogroup_index -= 1
        closest_different_left_orthogroup = get_left_orthogroup(spec, left_orthogroup_index)
    return closest_different_left_orthogroup


def get_closest_different_right_orthogroup(spec, index, right_orthogroup):
    right_orthogroup_index = index+1
    closest_different_right_orthogroup = get_right_orthogroup(spec, right_orthogroup_index)
    while closest_different_right_orthogroup == right_orthogroup:
        right_orthogroup_index += 1
        closest_different_right_orthogroup = get_right_orthogroup(spec, right_orthogroup_index)
    return closest_different_right_orthogroup


synteny_counts = {}

i = 0
n_ref_spec = len(ordered_gff_species_2_orthogroups.keys())
for reference_species in ordered_gff_species_2_orthogroups.keys():
    i += 1
    print('\nComputing SYN with ' + reference_species + ' as reference species. ' + str(i) + '/' + str(n_ref_spec))

    n_reference_orthogroups = len(ordered_gff_species_2_orthogroups[reference_species])
    is_left_duplicated = False
    is_right_duplicated = False

    j = 0
    n_genes = len(ordered_gff_species_2_orthogroups[reference_species])
    reference_orthogroup_index = 0
    for reference_orthogroup in ordered_gff_species_2_orthogroups[reference_species]:
        j += 1
        if j%100 == 0:
            progress = round(100*j/n_genes, 2)
            print('... ' + str(progress) + ' %')

        left_reference_orthogroup = get_left_orthogroup(reference_species, reference_orthogroup_index)
        right_reference_orthogroup = get_right_orthogroup(reference_species, reference_orthogroup_index)

        if left_reference_orthogroup == reference_orthogroup and left_reference_orthogroup != None:
            closest_different_left_reference_orthogroup = get_closest_different_left_orthogroup(reference_species, reference_orthogroup_index, left_reference_orthogroup)
        else:
            closest_different_left_reference_orthogroup = left_reference_orthogroup

        if right_reference_orthogroup == reference_orthogroup and right_reference_orthogroup != None:
            closest_different_right_reference_orthogroup = get_closest_different_right_orthogroup(reference_species, reference_orthogroup_index, right_reference_orthogroup)
        else:
            closest_different_right_reference_orthogroup = right_reference_orthogroup


        SYN = 0
        for target_species in ordered_gff_orthogroups_2_species[reference_orthogroup]:

            if target_species != reference_species:
                target_orthogroup_index = ordered_gff_species_2_orthogroups[target_species].index(reference_orthogroup)
                left_target_orthogroup = get_left_orthogroup(target_species, target_orthogroup_index)
                right_target_orthogroup = get_right_orthogroup(target_species, target_orthogroup_index)


                if left_target_orthogroup == reference_orthogroup and left_target_orthogroup != None:
                    closest_different_left_target_orthogroup = get_closest_different_left_orthogroup(target_species, target_orthogroup_index, left_target_orthogroup)
                else:
                    closest_different_left_target_orthogroup = left_target_orthogroup

                if right_target_orthogroup == reference_orthogroup and right_target_orthogroup != None:
                    closest_different_right_target_orthogroup = get_closest_different_right_orthogroup(target_species, target_orthogroup_index, right_target_orthogroup)
                else:
                    closest_different_right_target_orthogroup = right_target_orthogroup


                if left_target_orthogroup == left_reference_orthogroup:
                    SYN += 0.5
                elif left_target_orthogroup == closest_different_left_reference_orthogroup:
                    SYN += 0.5
                elif closest_different_left_target_orthogroup == closest_different_left_reference_orthogroup:
                    SYN += 0.5
                elif closest_different_left_target_orthogroup == left_reference_orthogroup:
                    SYN += 0.5

                if right_target_orthogroup == right_reference_orthogroup:
                    SYN += 0.5
                elif right_target_orthogroup == closest_different_right_reference_orthogroup:
                    SYN += 0.5
                elif closest_different_right_target_orthogroup == closest_different_right_reference_orthogroup:
                    SYN += 0.5
                elif closest_different_right_target_orthogroup == right_reference_orthogroup:
                    SYN += 0.5
            else:
                continue



        if reference_orthogroup not in synteny_counts.keys():
            synteny_counts[reference_orthogroup] = [SYN]
        else:
            synteny_counts[reference_orthogroup].append(SYN)

        reference_orthogroup_index += 1


# Process output files
pickle.dump(synteny_counts, output_synteny_dict, protocol=pickle.HIGHEST_PROTOCOL)

output_synteny_counts.write('orthogroup\tgff_n_species\tmax_syn\tavg_syn\tmed_syn\tsd_synn')
for orthogroup in sorted(synteny_counts.keys()):
    n = len(ordered_gff_orthogroups_2_species[orthogroup])
    max_counts = max(synteny_counts[orthogroup])
    avg = statistics.mean(synteny_counts[orthogroup])
    med = statistics.median(synteny_counts[orthogroup])
    if len(synteny_counts[orthogroup]) > 1:
        sd = statistics.stdev(synteny_counts[orthogroup])
    else:
        sd = 'NA'
    output_synteny_counts.write(f"{orthogroup}\t{n}\t{max_counts}\t{avg}\t{med}\t{sd}\n")

# Close files
output_synteny_dict.close()
output_synteny_counts.close()
