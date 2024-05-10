from typing import List, Union
import numpy.typing as npt
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
import csv
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
from sklearn.model_selection import train_test_split
from sklearn.model_selection import KFold
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import mean_squared_error

def readCSV(csvfile):
    with open("4_compound_machine_learning/"+csvfile) as f:
        reader = csv.DictReader(f)
        data = {}
        for r in reader:
            key = r['No.']
            value = {'SMILES':r['SMILES'].strip('\t\n'),'LogP app':float(r['LogP app'])}
            data[key]=value
        return data
def draw_molecule(csvfile: str) -> None:
    # 課題 4-1
    data = readCSV(csvfile)
    m = Chem.MolFromSmiles(data['6']["SMILES"])
    Draw.MolToFile(m,'4_compound_machine_learning/images/CHEMBL540227.png')
    pass

def create_2d_descriptors(smiles: str) -> Union[npt.NDArray[np.float_], List[float]]:
    # 課題 4-2
    m = Chem.MolFromSmiles(smiles)
    res_dir:dict = Descriptors.CalcMolDescriptors(m)
    res = []
    for key in res_dir.keys():
        res.append(res_dir[key])
    return res

def predict_logpapp(csvfile: str) -> Union[npt.NDArray[np.float_], pd.Series, List[float]]:
    # 課題 4-3
    np.random.seed(0) # 出力を固定するためにseedを指定
    rfr = RandomForestRegressor(random_state=0) # 出力を固定するためにrandom_stateを指定
    data = readCSV(csvfile)
    X = []
    y = []
    for key in data.keys():
        X.append(create_2d_descriptors(data[key]["SMILES"]))
        y.append(data[key]["LogP app"])
    X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=700, random_state=0)
    rfr.fit(X_train, y_train)
    return rfr.predict(X_test)

def grid_search(csvfile: str) -> float:
    # 課題 4-4
    # こちらも出力を固定するためにseedやrandom_stateを指定すること
    np.random.seed(0)
    def param():
        ret = {
            'n_estimators':[100, 200, 400],
            'max_depth':[5, 10, 15],
        }
        return ret

    # # Xは説明変数、yは目的変数
    data = readCSV(csvfile)
    X = []
    y = []
    for key in data.keys():
        X.append(create_2d_descriptors(data[key]["SMILES"]))
        y.append(data[key]["LogP app"])
    X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=700, random_state=0)
    rfr = RandomForestRegressor(random_state=0)
    gscv = GridSearchCV(rfr, param(), cv=4, verbose=2)
    gscv.fit(X_train, y_train)
    pv = gscv.predict(X_test)
    return float(np.sqrt(mean_squared_error(pv,y_test)))

if __name__ == "__main__":
    smiles = "C(=O)(c1ccc(OCCCCCC)cc1)CCNc1cc(Cl)ccc1"
    filepath = "data/fukunishi_data.csv"
    # 課題 4-1
    draw_molecule(filepath)
    # 課題 4-2
    print(create_2d_descriptors(smiles))
    # 課題 4-3
    print(predict_logpapp(filepath))
    # 課題 4-4
    print(grid_search(filepath))
