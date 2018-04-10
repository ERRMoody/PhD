import pandas as pd

F123 = pd.read_csv ('nameid.csv')

df1 = F123.ix[:,['0', '2']]

df1.columns = ['Species', 'LocalID']

df1.to_csv('nameid.csv', index = False)

result = df1.sort_values (by='Species', ascending=1)



result1 = result.drop_duplicates(['Species'])

result2 = result1.drop('LocalID', 1)

result2.to_csv('specieslist.csv', index = False)