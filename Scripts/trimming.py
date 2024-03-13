import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

tab = pd.read_csv("openrdp_results")

filtered_tab = tab[tab['Pvalue'] < 0.05]
print(filtered_tab)

density = gaussian_kde(filtered_tab['Pvalue'])
x = np.linspace(0, max(filtered_tab['Pvalue']), 1000)
plt.title('PDF of p-value on the curated results table')
plt.xlabel('P-value')
plt.ylabel('Probability Density')
plt.plot(x, density(x))
plt.show()

max_gencv = filtered_tab[filtered_tab['Method'] == "Geneconv"]['Pvalue'].max()
max_boot = filtered_tab[filtered_tab['Method'] == "Bootscan"]['Pvalue'].max()
max_maxchi = filtered_tab[filtered_tab['Method'] == "Maxchi"]['Pvalue'].max()
max_sisc = filtered_tab[filtered_tab['Method'] == "Siscan"]['Pvalue'].max()
max_chim = filtered_tab[filtered_tab['Method'] == "Chimaera"]['Pvalue'].max()
max_3seq = filtered_tab[filtered_tab['Method'] == "3Seq"]['Pvalue'].max()
max_rdp = filtered_tab[filtered_tab['Method'] == "Rdp"]['Pvalue'].max()

print('\n----------------------------------------')
print(f'\nGeneconv highest p-value: {max_gencv}')
print(f'Bootscan highest p-value: {max_boot}')
print(f'MaxChi highest p-value: {max_maxchi}')
print(f'SiScan highest p-value: {max_sisc}')
print(f'Chimaera highest p-value: {max_chim}')
print(f'3Seq highest p-value: {max_3seq}')
print(f'RDP highest p-value: {max_rdp}\n')

filtered_tab.to_csv("curated_openrdp_results", index=False)

print('The OpenRDP results were curated to p-value < 0.05 and saved on: "curated_openrdp_results"\n')
