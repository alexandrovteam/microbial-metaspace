import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from metaspace import sm_annotation_utils as smau
from pathlib import Path


def clip(img, bottom=0, top=99):
    """Remove bottom and top percentiles of image pixels"""
    lo = np.min(img)
    lo_thresh = np.max(img) / 256
    hi = np.percentile(img[img >= lo_thresh], top)
    img_clipped = np.clip(img, lo, hi)
    return img_clipped


def rescale(data, bottom=0, top=255):
    """Rescale data to lie between [min, max] values"""
    compressed = (data - np.min(data)) / np.max(data - np.min(data))
    rescaled = compressed * (top - bottom) + bottom
    return rescaled


# User inputs
p = Path(r"C:\Users\Veronica\Documents\LAB\projects\microbes\microbial_metaspace")
p_npa = p / r"np_atlas_2019_08_corrected.tsv"
p_proj = p / r"microbial_metaspace"
NPA_version = 'v.2019_08'
ds_id = '2020-01-07_11h15m21s'
fdr = 0.5

# Load NPA database
npa = pd.read_csv(str(p_npa), sep='\t', header=0)

# Get results table from Metaspace
ms = smau.SMInstance()
d = ms.dataset(id=ds_id)
databases = d.databases
unfiltered_ms_res = pd.DataFrame()
filtered_ms_res = pd.DataFrame()

for db in databases:
    if db == 'NPA-2019-08':
        unfiltered_ms_res = pd.concat(
            [unfiltered_ms_res, ms.msm_scores([d], d.annotations(fdr, database=db), db_name=db).T])
        filtered_ms_res = unfiltered_ms_res.drop_duplicates()

# Filter in batch mode
batch_f = pd.read_csv(str(p_proj / 'batch_filtration.csv'))
for i in batch_f.index:
    folder = batch_f.loc[i, 'folder']
    subfolder = batch_f.loc[i, 'subfolder']
    origin = [batch_f.loc[i, 'type']]
    genus = [batch_f.loc[i, 'genus']]
    species = [batch_f.loc[i, 'species']]
    if species == ['no_filter']: species = ['']
    species_match = batch_f.loc[i, 'filter']

    p_out = p_proj / str(folder) / str(subfolder)
    p_out.mkdir(parents=True, exist_ok=True)

# Select results that correspond only to organism of interest in NPA
    by_metadata_ms_res = pd.DataFrame()
    for formula, adduct in filtered_ms_res.index:
        match = npa.loc[npa['Molecular Formula'] == formula]
        if origin != ['']:
            filter1 = match.loc[match['Origin Type'].isin(origin)]
        else:
            print('Error message: organism type not specified')
            break
        if genus != ['']:
            filter2 = filter1.loc[filter1['Genus'].isin(genus)]
        else:
            filter2 = filter1
        if (species != ['']) & (species_match == 'exact_search'):
            filter3 = filter2.loc[filter2['Origin Species'].isin(species)]
        elif (species != ['']) & (species_match == 'substring_search'):
            filter3 = filter2.loc[filter2['Origin Species'].str.contains('|'.join(species))]
        else:
            filter3 = filter2
        by_metadata_ms_res = pd.concat([by_metadata_ms_res, filter3])

    by_metadata_ms_res[['NPAID', 'Names', 'Molecular Formula', 'Accurate Mass', 'Origin Type', 'Genus', 'Origin Species', 'NPA URL']].drop_duplicates().to_csv(str(p_out / 'filtered_molecules.csv'), index=False)
    stat1 = len(set(by_metadata_ms_res['Molecular Formula']))

    # Visualise selected molecules
    ion_hits = pd.DataFrame(columns=['formula', 'adduct', 'msm'])
    all_res = filtered_ms_res.reset_index()

    for formula in set(by_metadata_ms_res['Molecular Formula']):
        instances = all_res.loc[all_res['formula'] == formula]
        for j in instances.index:
            adduct = instances.loc[j, 'adduct']
            msm = instances.loc[j, instances.columns[2]]
            ion_hits.loc[j] = [formula, adduct, msm]

    ion_hits.to_csv(str(p_out / 'ion_images.csv'), index=False)
    stat2 = len(ion_hits)

    # Save images
    for k in ion_hits.index:
        formula = ion_hits.loc[k, 'formula']
        adduct = ion_hits.loc[k, 'adduct']
        msm = ion_hits.loc[k, 'msm']
        image = d.isotope_images(formula, adduct)[0]  # main isotope image
        plt.imsave(str(p_out / (formula + adduct + '.png')), rescale(clip(image)))

    # Write a text file with parameters and results
    txt = open(p_out / "analysis_info.txt", "a")
    txt.write('Dataset ID: {} \n\
    FDR: {} \n\
    NPA_downloaded: {} \n\
    Origin Type: {} \n\
    Genus: {} \n\
    Origin Species: {} \n\
    Species were filtered using {} \n\
    {} unique molecular formulas match filter criteria \n\
    These correspond to {} ion images (e.g. different adducts)'.format(ds_id, fdr, NPA_version, origin, genus, species, species_match, stat1, stat2))
    txt.close()