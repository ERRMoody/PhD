import pandas as pd
df = (pd.read_csv('fixed5.csv'))
df2 = (pd.read_csv('UBERLIST.csv'))

df3 = pd.merge(df,df2, on='Strain', how='left')

df3.to_csv('finalfixedencompassing.csv')
