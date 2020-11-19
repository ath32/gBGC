''' Perform logistic regression analysis. Requires CSV containing PxAbundance and Recombination data. '''

### IMPORT ###
import pandas as pd
import numpy as np
from sklearn import preprocessing
import matplotlib.pyplot as plt
plt.rc("font", size=14)
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import RFE
import seaborn as sns
sns.set(style="white")
sns.set(style="whitegrid", color_codes=True)
import statsmodels.api as sm

### CHOOSE SOURCE FOLDER ###

source = 'Expression/gene_information.csv'

### MAIN ###

def main():

    non_switches = ['TAA_stay', 'TGA_stay', 'TAG_stay']

    #Import data
    data = pd.read_csv(source, header = 0)
    data = data.set_index('Gene_id')
    data = data[data['Switch'] != 'TAA_stay']
    data = data[data['Switch'] != 'TGA_stay']
    data = data[data['Switch'] != 'TAG_stay']
    print (data.head())

    #Convert switches to dummies
    df = pd.get_dummies(data, columns=['Switch'])
    print (df.head())
    print(df.groupby('Switch_TAA>TGA').mean())

    #Choose columns
    dependent = 'Switch_TAA>TGA'
    to_include = ['PxAbundance']
    X = df[to_include]
    y = df[dependent]

    #Model
    logreg = LogisticRegression()
    logit_model=sm.Logit(y,X)
    result=logit_model.fit()
    print(result.summary2())

    #Get coefficient and intercept
    clf = LogisticRegression(random_state=0).fit(X, y)
    print(clf.coef_, clf.intercept_)

    #Get top feature
    predictors = X
    selector = RFE(clf, n_features_to_select=1)
    selector = selector.fit(predictors, y)
    order = selector.ranking_

    feature_ranks = []
    for i in order:
        feature_ranks.append(f"{i}. {to_include[i-1]}")
    print (feature_ranks)

### RUN ###

if __name__ == '__main__':
    main()
