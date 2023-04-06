#!/usr/bin/env python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import cosine_similarity


# Read disease/phenotype data

df = pd.read_table('disease_name_hpo_list.tsv')
diseases = df['hpo']

# Create TFIDF disease vectors
vectorizer = TfidfVectorizer()
disease_vec = vectorizer.fit_transform(diseases)

# Create TFIDF query vector
query = (
    'HP_0012531,HP_0001263,HP_0000152,HP_0000234,HP_0000707,HP_0001250,HP_0001298,HP_0002011,HP_0002384,HP_0012443,HP_0001649,HP_0002349,HP_0007359,HP_0004372,HP_0001254,HP_0002329,HP_0004897,HP_0011024,HP_0025031,HP_0002060,HP_0003011,HP_0003198,HP_0003458,HP_0004303,HP_0012759,HP_0040290,HP_0200134,HP_0000708,HP_0000818,HP_0001438,HP_0001742,HP_0001943,HP_0003202,HP_0011013,HP_0000765,HP_0001626,HP_0001627,HP_0001697,HP_0002910,HP_0012393,HP_0045026,HP_0000153,HP_0000163,HP_0000486,HP_0000737,HP_0001270,HP_0001633,HP_0001653,HP_0001654,HP_0001698,HP_0001941,HP_0002242,HP_0002576,HP_0002814,HP_0003128,HP_0004796,HP_0005214,HP_0006705,HP_0010676,HP_0010877,HP_0040066,HP_0003111,HP_0012378,HP_0000577,HP_0001549,HP_0001637,HP_0001942,HP_0002086,HP_0002244,HP_0011703,HP_0012735,HP_0000565,HP_0000716,HP_0000739,HP_0000750,HP_0000958,HP_0001251,HP_0001252,HP_0001289,HP_0001290,HP_0001392,HP_0001635,HP_0001638,HP_0001946,HP_0001974,HP_0002087,HP_0002104,HP_0002133,HP_0002474,HP_0002480,HP_0002684,HP_0002715,HP_0002813,HP_0003201,HP_0003558,HP_0008119,HP_0009900,HP_0011105,HP_0011343,HP_0011439,HP_0011663,HP_0012103,HP_0012639,HP_0012847,HP_0000019,HP_0000194,HP_0000271,HP_0000315,HP_0000366,HP_0000478,HP_0000490,HP_0000504,HP_0000520,HP_0000593,HP_0000646,HP_0000819,HP_0000995,HP_0001276,HP_0001367,HP_0001623,HP_0001701,HP_0002013,HP_0002015,HP_0002067,HP_0002079,HP_0002098,HP_0002353,HP_0002375,HP_0002487,HP_0002509,HP_0002789,HP_0002883,HP_0003040,HP_0003270,HP_0003698,HP_0003764,HP_0004328,HP_0005335,HP_0005957,HP_0005978,HP_0006872,HP_0006895,HP_0007239,HP_0010808,HP_0011885,HP_0011886,HP_0011947,HP_0012819,HP_0025085,HP_0025406,HP_0030048,HP_0031417,HP_0100649,HP_0100962,HP_0200055,HP_0200136'
)
re.sub(query, ':', '_')
disease_name = 'Metabolic encephalomyopathic crises, recurrent, with rhabdomyolysis, cardiac arrhythmias, and neurodegeneration'
query_vec = vectorizer.transform([query])

# Calculate cosine similarity between query and disease vectors
cos_sim = cosine_similarity(disease_vec, query_vec)

# Create DataFrame of cosine similarity matches
cs_df = pd.DataFrame(cos_sim, index=df['name'], columns=['similarity'])
cs_df.index.rename('disease', inplace=True)

# Sort DataFrame of matches by descending similarity
cs_df.sort_values(by='similarity', ascending=False, inplace=True)


cs_dfx = cs_df.reset_index()
cs_dfx.index[cs_dfx['disease'] == dis_name].values.tolist()


# Bar plot of top 25 disease candidates
# plt.figure(figsize=(15,8))
# sns.barplot(data=cs_df.iloc[0:25,:].reset_index(), x='similarity', y='disease')

# KDE plot of cosing similarity scores for all diseases
# sns.kdeplot(data=cs_df)
