#!/usr/bin/env python3

import pandas as pd
import json
import requests
import sys
from datetime import datetime
import time

host = 'https://www.orthodb.org/'
search = 'group?id='
outpath = '/home/lruzzant/data/evolrates/6656_evolrates.tab'
part_step = 500

odb = pd.read_csv('/home/lruzzant/data/odb10v1_OG2genes_6656.tab', sep='\t', header=None)
orthogroups = sorted(set(odb[0]))
n_orthogroups = len(orthogroups)
evol_rates = {}

i_0 = 0 # Last orthgroup row index at which the part download was interrupted
i = i_0

print(f'Host: {host} ; search: {search}$orthogroup_id')
print( f'{n_orthogroups} orthogroups to process')
print(f'START: {datetime.now()}') ## Get starting time of script

for orthogroup in orthogroups[i_0:]:
    query = host+search+orthogroup
    #print(f'Your query: {query}')
    response = requests.get(query)
    while response.status_code == 504: ## 504 Gateway Timeout error is an HTTP status code that means that one server did not receive a timely response from another server. So we try again until the query works.
        response = requests.get(query)
        print(f'Gateway Timeout error. Response code: {response.status_code}\n')
        print('Trying again in 5 seconds ...')
        time.sleep(5)
    if response.status_code != 200: ## Any other response code than 200 is an error
        print(f'Bad query: {query}\n')
        print(f'Response code: {response.status_code}\n')
        print(f'Response content: {response.text}\n')
        sys.exit() ## Quits the script, because we have a bad query

    jsonData = json.loads(response.content.decode('utf-8')) ## Convert the response we got from the api to python json format
    EVR = jsonData['data']['evolutionary_rate']

    #print(f'orthogroup: {orthogroup}; evolutionary_rate: {EVR}')

    evol_rates[orthogroup] = EVR

    i += 1
    if i%10 == 0:
        print('.', end='', flush=True)
    if i%500 == 0:
        print( f'{i}/{n_orthogroups} --> {round(i/n_orthogroups*100, 2)}%')

    # Writing to file parts of download, by part_step size
    if i in range(part_step,n_orthogroups+part_step,part_step):
        outfile_part = outpath+str(i_0+1)+'_'+str(i)
        with open(outfile_part, 'w') as f:
            for orthogroup in sorted(evol_rates.keys()):
                f.write(orthogroup + '\t' + str(evol_rates[orthogroup]) + '\n')
        print(f'Part download saved in {outfile_part}')

# And to get all of the last remaining orthogroups:
outfile_part = outpath+str(i_0+1)+'_'+str(i)
with open(outpath, 'w') as f:
    outfile_part = outfile_part+str(i_0+1)+'_'+str(i)
    for orthogroup in sorted(evol_rates.keys()):
        f.write(orthogroup + '\t' + str(evol_rates[orthogroup]) + '\n')

print(f'File written to: {outpath}')
