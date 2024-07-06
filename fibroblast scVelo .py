#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd


# In[3]:


base_path = "/Users/eaiscui/Desktop/fibroblastforscVelo/"

# Load the sparse matrix
X = io.mmread(base_path + "counts.mtx")

# Create an AnnData object
adata = anndata.AnnData(
    X=X.transpose().tocsr()
)

# Load cell metadata
cell_meta = pd.read_csv(base_path + "metadata.csv")

# Load gene names
with open(base_path + "gene_names.csv", 'r') as f:
    gene_names = f.read().splitlines()


# In[4]:


# set anndata observations and index obs by barcodes, var by gene names
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode']
adata.var.index = gene_names

# load dimensional reduction:
# Load dimensional reduction
pca = pd.read_csv(base_path + "pca.csv")
pca.index = adata.obs.index

# set pca and umap
adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_UMAP'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T


# In[5]:


adata.obsm['X_umap'] = adata.obsm['X_UMAP']


# In[6]:


# plot a UMAP colored by sampleID to test:
sc.pl.umap(adata, color=['cell_type'], frameon=False, save=True)


# In[7]:


# save dataset as anndata format
adata.write('my_data2.h5ad')


# In[8]:


import os
print(os.getcwd())


# In[9]:


# reload dataset
adata = sc.read_h5ad('my_data2.h5ad')


# In[2]:


import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad


# In[3]:


scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
cr.settings.verbosity = 2


# In[4]:


path_to_loom_files = "/Volumes/1TB/MOUSEAV/GSM5457081/velocyto/"

ldata1 = scv.read(path_to_loom_files + 'GSM5457081.loom', cache=True)


# In[12]:


adata = sc.read_h5ad('my_data2.h5ad')
# load loom files for spliced/unspliced matrices for each sample:
path_to_loom_files = "/Users/eaiscui/Desktop/loomfiles/"

ldata1 = scv.read(path_to_loom_files + 'Dminus_H2B.loom', cache=True)
ldata2 = scv.read(path_to_loom_files + 'Dplus_H2B.loom', cache=True)


# In[13]:


# rename barcodes in order to merge:
barcodes = [bc.split(':')[1] for bc in ldata1.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_10' for bc in barcodes]
ldata1.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata2.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_11' for bc in barcodes]
ldata2.obs.index = barcodes


# In[14]:


# make variable names unique
ldata1.var_names_make_unique()
ldata2.var_names_make_unique()


# In[15]:


# concatenate the three loom
ldata = ldata1.concatenate([ldata2])


# In[16]:


scv.utils.clean_obs_names(adata)
scv.utils.clean_obs_names(ldata)

adata = scv.utils.merge(adata, ldata)


# In[17]:


# plot umap to check
sc.pl.umap(adata, color='cell_type', frameon=False, legend_loc='on data', title='', save='_cell_types.pdf')


# In[5]:


scv.pl.proportions(ldata1)


# In[19]:


# pre-process
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)


# In[20]:


# compute velocity
scv.tl.velocity(adata, mode='stochastic')


# In[21]:


scv.tl.velocity_graph(adata)


# In[22]:


scv.pl.velocity_embedding(adata, basis='umap', frameon=False, save='embedding.pdf', color=['cell_type'])


# In[23]:


scv.pl.velocity_embedding_grid(adata, basis='umap', color='cell_type', save='embedding_grid.pdf', title='', scale=0.25)


# In[24]:


scv.pl.velocity_embedding_stream(adata, basis='umap', color=['cell_type'], save = '/Users/eaiscui/Desktop/fibroblastforscVelo/basicflow.pdf')


# In[25]:


adata.obs['cell_type']


# In[ ]:


# plot velocity of a selected gene
scv.pl.velocity(adata, var_names=['Col1a1'], color='cell_type', save = '/Users/eaiscui/Desktop/fibroblastforscVelo/col1a1.pdf')


# In[ ]:


# plot velocity of a selected gene
scv.pl.velocity(adata, var_names=['Col3a1'], color='cell_type', save = '/Users/eaiscui/Desktop/fibroblastforscVelo/col3a1.pdf')


# In[ ]:


# plot velocity of a selected gene
scv.pl.velocity(adata, var_names=['Eln'], color='cell_type', save = '/Users/eaiscui/Desktop/fibroblastforscVelo/eln.pdf')


# In[ ]:


# plot velocity of a selected gene
scv.pl.velocity(adata, var_names=['Ccl2'], color='cell_type', save = '/Users/eaiscui/Desktop/fibroblastforscVelo/ccl2.pdf')


# In[ ]:


# plot velocity of a selected gene
scv.pl.velocity(adata, var_names=['Ugdh'], color='cell_type', save = '/Users/eaiscui/Desktop/fibroblastforscVelo/ugdh.pdf')


# In[ ]:


# plot velocity of a selected gene
scv.pl.velocity(adata, var_names=['Col1a2'], color='cell_type', save = '/Users/eaiscui/Desktop/fibroblastforscVelo/col1a2.pdf')


# In[ ]:


# plot velocity of a selected gene
scv.pl.velocity(adata, var_names=['Has2'], color='cell_type', save = '/Users/eaiscui/Desktop/fibroblastforscVelo/Has2.pdf')


# In[ ]:


# plot velocity of a selected gene
scv.pl.velocity(adata, var_names=['Ugp2'], color='cell_type', save = '/Users/eaiscui/Desktop/fibroblastforscVelo/Ugp2.pdf')


# In[26]:


scv.pl.velocity_embedding(adata, arrow_length=3, arrow_size=2, dpi=120, color='cell_type')


# In[27]:


scv.pl.velocity(adata, ['Col1a1', 'Col1a2', 'Col3a1', 'Col5a1', 'Col8a1','Eln','Ugdh','Ugp2','Has2'], ncols=1, color='cell_type', save = '/Users/eaiscui/Desktop/fibroblastforscVelo/all.pdf')


# In[28]:


scv.tl.rank_velocity_genes(adata, groupby='cell_type', min_corr=.3)

df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head()


# In[29]:


scv.tl.recover_dynamics(adata)


# In[30]:


scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)


# In[31]:


scv.pl.velocity_embedding_stream(adata, basis='umap',color='cell_type', save = '/Users/eaiscui/Desktop/fibroblastforscVelo/dynamicflow.pdf')


# In[32]:


df = adata.var
df = df[(df['fit_likelihood'] > .1) & df['velocity_genes'] == True]

kwargs = dict(xscale='log', fontsize=16)
with scv.GridSpec(ncols=3) as pl:
    pl.hist(df['fit_alpha'], xlabel='transcription rate', **kwargs)
    pl.hist(df['fit_beta'] * df['fit_scaling'], xlabel='splicing rate', xticks=[.1, .4, 1], **kwargs)
    pl.hist(df['fit_gamma'], xlabel='degradation rate', xticks=[.1, .4, 1], **kwargs)

scv.get_df(adata, 'fit*', dropna=True).head()


# In[33]:


scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80, save = '/Users/eaiscui/Desktop/fibroblastforscVelo/latenttime.pdf')


# In[34]:


top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata, var_names=top_genes, sortby='latent_time', col_color='cell_type', n_convolve=100, yticklabels=True,figsize=(4,50),save = '/Users/eaiscui/Desktop/fibroblastforscVelo/heatmapflow.pdf')


# In[35]:


top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata, var_names=top_genes, sortby='latent_time', col_color='cell_type', n_convolve=100,save = '/Users/eaiscui/Desktop/fibroblastforscVelo/heatmapflow.pdf')


# In[36]:


top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index
scv.pl.scatter(adata, basis=top_genes[:15], ncols=5,color='cell_type',frameon=False, save = '/Users/eaiscui/Desktop/fibroblastforscVelo/likelihood.pdf')


# In[37]:


var_names = ['Col1a1', 'Col1a2', 'Col3a1', 'Col5a1', 'Col8a1','Eln','Ugdh','Ugp2','Has2']
scv.pl.scatter(adata, var_names,color='cell_type', frameon=False, save = '/Users/eaiscui/Desktop/fibroblastforscVelo/gene1.pdf')
scv.pl.scatter(adata, x='latent_time', y=var_names,color='cell_type', frameon=False, save = '/Users/eaiscui/Desktop/fibroblastforscVelo/latenttime2.pdf')


# In[38]:


scv.tl.rank_dynamical_genes(adata, groupby='cell_type')
df = scv.get_df(adata, 'rank_dynamical_genes/names')
df.head(5)


# In[ ]:





# In[ ]:





# In[ ]:





# In[39]:


vk = cr.kernels.VelocityKernel(adata)


# In[40]:


vk.compute_transition_matrix()


# In[41]:


ck = cr.kernels.ConnectivityKernel(adata)
ck.compute_transition_matrix()

combined_kernel = 0.8 * vk + 0.2 * ck


# In[42]:


print(combined_kernel)


# In[48]:


vk.plot_projection(basis='umap',color='cell_type', save = '/Users/eaiscui/Desktop/fibroblastforscVelo/cellrankflow.pdf')


# In[54]:


vk.plot_random_walks(start_ixs={"cell_type": "1-Fibroblast"}, max_iter=200, seed=0,basis='umap',save = '/Users/eaiscui/Desktop/fibroblastforscVelo/randomwalkfrom1.pdf')


# In[55]:


vk.plot_random_walks(start_ixs={"cell_type": "2-Fibroblast"}, max_iter=200, seed=0,basis='umap',save = '/Users/eaiscui/Desktop/fibroblastforscVelo/randomwalkfrom2.pdf')


# #### vk.plot_random_walks(start_ixs={"cell_type": "1-fib"}, max_iter=200, seed=0,basis='umap')

# In[56]:


vk.plot_random_walks(start_ixs={"cell_type": "3-Fibroblast"}, max_iter=200, seed=0,basis='umap',save = '/Users/eaiscui/Desktop/fibroblastforscVelo/randomwalkfrom3.pdf')


# In[57]:


vk.plot_random_walks(start_ixs={"cell_type": "4-Fibroblast"}, max_iter=200, seed=0,basis='umap',save = '/Users/eaiscui/Desktop/fibroblastforscVelo/randomwalkfrom4.pdf')


# In[58]:


vk.plot_random_walks(start_ixs={"cell_type": "5-Fibroblast"}, max_iter=200, seed=0,basis='umap',save = '/Users/eaiscui/Desktop/fibroblastforscVelo/randomwalkfrom5.pdf')


# In[59]:


vk.plot_random_walks(start_ixs={"cell_type": "6-Fibroblast"}, max_iter=200, seed=0,basis='umap',save = '/Users/eaiscui/Desktop/fibroblastforscVelo/randomwalkfrom6.pdf')


# In[60]:


vk.plot_random_walks(start_ixs={"cell_type": "7-Fibroblast"}, max_iter=200, seed=0,basis='umap',save = '/Users/eaiscui/Desktop/fibroblastforscVelo/randomwalkfrom7.pdf')


# In[61]:


vk.plot_random_walks(start_ixs={"cell_type": "8-Fibroblast"}, max_iter=200, seed=0,basis='umap',save = '/Users/eaiscui/Desktop/fibroblastforscVelo/randomwalkfrom8.pdf')


# In[62]:


g = cr.estimators.GPCCA(vk)
print(g)


# In[65]:


g.fit(cluster_key="cell_type", n_states=[4, 12])
g.plot_macrostates(which="all", discrete=True, legend_loc="right", s=100,save = '/Users/eaiscui/Desktop/fibroblastforscVelo/celltypestate.pdf')


# In[66]:


g.predict_terminal_states()
g.plot_macrostates(which="terminal", legend_loc="right", s=100,save = '/Users/eaiscui/Desktop/fibroblastforscVelo/terminal1.pdf')


# In[68]:


g.plot_macrostates(which="terminal", discrete=False,save = '/Users/eaiscui/Desktop/fibroblastforscVelo/terminalstate.pdf')


# In[69]:


g.predict_initial_states()
g.plot_macrostates(which="initial", legend_loc="right", s=100,save = '/Users/eaiscui/Desktop/fibroblastforscVelo/initialstate.pdf')


# In[70]:


g


# In[72]:


sc.pl.embedding(
    adata,
    basis="umap",
    color=["Col1a1", "Fap", "Nt5e", "Thy1", "Eng", "Cd44"],
    size=50,
)


# In[73]:


g.plot_coarse_T()


# In[74]:


g2 = cr.estimators.GPCCA(vk)
print(g2)


# In[75]:


g2.compute_schur()
g2.plot_spectrum(real_only=True)


# In[76]:


g2.compute_macrostates(n_states=11, cluster_key="cell_type")
g2.plot_macrostates(which="all", legend_loc="right", s=100)


# In[77]:


g2.plot_macrostate_composition(key="cell_type", figsize=(7, 4))


# In[78]:


g2.plot_coarse_T(annotate=False)


# In[79]:


g2.predict_terminal_states()
g2.plot_macrostates(which="terminal", legend_loc="right", s=100)


# In[80]:


g2.predict_initial_states()
g2.plot_macrostates(which="initial", s=100)


# In[84]:


g.compute_fate_probabilities()
g.plot_fate_probabilities(same_plot=False,save = '/Users/eaiscui/Desktop/fibroblastforscVelo/termianlfate1.pdf')


# In[85]:


g.plot_fate_probabilities(same_plot=True,save = '/Users/eaiscui/Desktop/fibroblastforscVelo/termianlfate2.pdf')


# In[86]:


g.fit(cluster_key="cell_type", n_states=8)

g.set_terminal_states(states=["4-Fibroblast"])
g.plot_macrostates(which="terminal", legend_loc="right", size=100)


# In[87]:


g.compute_fate_probabilities()
g.plot_fate_probabilities(same_plot=False)


# In[88]:


vk.plot_projection(basis="umap",color="cell_type")


# In[90]:


top_genes = adata.var["fit_likelihood"].sort_values(ascending=False).index
scv.pl.scatter(adata, basis=top_genes[:10], ncols=5, frameon=False,color="cell_type",save = '/Users/eaiscui/Desktop/fibroblastforscVelo/likehoodgenes2.pdf')


# In[91]:


sc.tl.diffmap(adata)


# In[92]:


adata.obsm['X_diffmap'][:, 1].argmax()


# In[107]:


root_ixs = 4000  # has been found using `adata.obsm['X_diffmap'][:, 3].argmax()`
scv.pl.scatter(
    adata,
    basis="diffmap",
    c=["cell_type", root_ixs],
    legend_loc="right",
    components=["2, 3"],
)

adata.uns["iroot"] = root_ixs


# In[111]:


sc.tl.dpt(adata)
sc.pl.embedding(
    adata,
    basis="umap",
    color=["dpt_pseudotime"],
    color_map="gnuplot2"
)


# In[112]:


pk = cr.kernels.PseudotimeKernel(adata, time_key="dpt_pseudotime")
pk.compute_transition_matrix()

print(pk)


# In[120]:


pk.plot_projection(basis="umap", recompute=True, color="cell_type", save='/Users/eaiscui/Desktop/fibroblastforscVelo/11psedotime111.pdf')


# In[ ]:




