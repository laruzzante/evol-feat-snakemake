def get_MRCA_branchlength_from_species_list(species_list, MRCA_branchlengths_dict):

    AGE_max = 0

    for spec1 in sorted(species_list):
        for spec2 in sorted(species_list):
            if(spec1 != spec2):
                if MRCA_branchlengths_dict[spec1][spec2]:
                    AGE = MRCA_branchlengths_dict[spec1][spec2]
                elif MRCA_branchlengths_dict[spec2][spec1]:
                    AGE = MRCA_branchlengths_dict[spec2][spec1]
                else:
                    print('ERROR: species combination not present in MRCA branchlengths dictionary.')
                    sys.close()

                if AGE > AGE_max:
                    AGE_max = AGE

    if AGE_max == 0:
        print('ERROR: species list ' + ' '.join(species_list) + ' returns a MRCA branchlength of 0.')

    return AGE_max



def get_MRCA_ntips_from_species_list(species_list, mrca_ntips_dict):

    ntips_max = 0

    for spec1 in sorted(species_list):
        for spec2 in sorted(species_list):
            if(spec1 != spec2):
                if mrca_ntips_dict[spec1][spec2]:
                    ntips = mrca_ntips_dict[spec1][spec2]
                elif mrca_ntips_dict[spec2][spec1]:
                    ntips = mrca_ntips_dict[spec2][spec1]
                else:
                    print('ERROR: species combination not present in MRCA ntips dictionary.')
                    sys.close()

                if ntips > ntips_max:
                    ntips_max = ntips

    if ntips_max == 0:
        print('ERROR: species list ' + ' '.join(species_list) + ' returns a MRCA ntips of 0.')

    return ntips_max
