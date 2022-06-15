from urllib.request import urlopen
import json
import pandas as pd

df = pd.DataFrame(columns=['gene_curie','phenotypes','total_phenotypes','homologs','total_homologs','ortholog_phenotypes','total_ortholog_phenotypes'])

nf = 'not found'

with open('gene curies.tsv', 'r') as gc:
	next(gc)
	for line in gc:
		line = line.strip()
		row = line.split('\t')
		curie = row[1]
		print(curie)
		if curie == 'not found':
			df = df.append({'gene_curie':curie,'phenotypes':nf,'total_phenotypes':0,'homologs':nf,'total_homologs':0,'ortholog_phenotypes':nf,'total_ortholog_phenotypes':0}, ignore_index=True)
		else:
			phen_url = 'https://api.monarchinitiative.org/api/bioentity/gene/' + curie + '/phenotypes?rows=100&facet=false&unselect_evidence=false&exclude_automatic_assertions=false&fetch_objects=false&use_compact_associations=false&direct=false&direct_taxon=false'
			phen_raw = urlopen(phen_url)
			phen_json = json.loads(phen_raw.read())
			phens = []
			tot_phen = 0
			if phen_json['numFound'] != 0:
				tot_phen = len(phen_json['associations'])
				for p in phen_json['associations']:
					phen_label = p['object']['label']
					phen_id = p['object']['id']
					source = p['provided_by']
					phen = [phen_label,phen_id,source]
					phens.append(phen)
			else:
				phens.append(nf)
			hom_url = 'https://api.monarchinitiative.org/api/bioentity/gene/' + curie + '/homologs?rows=100&facet=false&unselect_evidence=false&exclude_automatic_assertions=false&fetch_objects=false&use_compact_associations=false&direct=false&direct_taxon=false'
			hom_raw = urlopen(hom_url)
			hom_json = json.loads(hom_raw.read())
			homs = []
			tot_hom = 0
			if hom_json['numFound'] != 0:
				tot_hom = len(hom_json['associations'])
				for h in hom_json['associations']:
					hom_label = h['object']['label']
					hom_id = h['object']['id']
					hom_taxon = h['object']['taxon']['label']
					hom_source = h['provided_by']
					hom = [hom_label,hom_id,hom_taxon,hom_source]
					homs.append(hom)
			else:
				homs.append(nf)
			ophen_url = 'https://api.monarchinitiative.org/api/bioentity/gene/' + curie + '/ortholog/phenotypes?rows=100&facet=false&unselect_evidence=false&exclude_automatic_assertions=false&fetch_objects=false&use_compact_associations=false&direct=false&direct_taxon=false'
			ophen_raw = urlopen(ophen_url)
			ophen_json = json.loads(ophen_raw.read())
			ophens = []
			tot_ophen = 0
			if ophen_json['numFound']!= 0:
				tot_ophen = len(ophen_json['associations'])
				for o in ophen_json['associations']:
					ophen_label = o['object']['label']
					ophen_id = o['object']['id']
					ophen_taxon = o['subject']['taxon']['label']
					ophen_source = o['provided_by']
					ophen = [ophen_label,ophen_id,ophen_taxon,ophen_source]
					ophens.append(ophen)
			else:
				ophens.append(nf)
			df = df.append({'gene_curie':curie,'phenotypes':phens,'total_phenotypes':tot_phen,'homologs':homs,'total_homologs':tot_hom,'ortholog_phenotypes':ophens,'total_ortholog_phenotypes':tot_ophen}, ignore_index=True)
df.to_csv('gene_data.tsv', sep="\t")