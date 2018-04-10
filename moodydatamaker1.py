import pandas as pd

df2 = pd.read_csv('fixedlistall.csv')

#percaa2 = df2.pivot_table('Percentage of Amino Acid', ['Strain', 'Domain','Actual OGT'], 'Amino Acid')
percaa2 = df2.pivot(index=['Strain','Alignment file'], columns='Amino Acid', values='Percentage of Amino Acid')


percaa2.to_csv('whatever.csv')
