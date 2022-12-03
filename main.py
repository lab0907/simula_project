import numpy as np
import networkx as nx
from scipy import sparse as sp
from matplotlib import pyplot as plt
from tqdm import tqdm
import pandas as pd
import random
import ssl

kinetics = [[0, 0, 0], [0, 0, 0.01], [0, 0.01, 0]]                      # 速率因子
mol = np.array([100000])                                                 # 分子数
group_type = np.array([[1, 1, 1]])                          # 每个分子的官能数
interval = 10                                                                  # 循环间隔
cycles = 6*(10**8)                                                     # 循环次数
column = ['cycle', 'P', 'Mn', 'Mw', 'P_theory', 'PMw']
k = 0.3


N_mol = sum(mol)
l_gro = len(group_type[0])
N_branch = range(l_gro)
network = sp.dok_matrix((N_mol, N_mol))
group = np.repeat(group_type, mol, axis=0)
sum_gro = np.array([np.sum(group == i) for i in range(1, np.max(group)+1)])
N_gro = sum_gro.copy()
data = pd.DataFrame(columns = column)
data_size = pd.DataFrame()
mol_list = range(N_mol)

def react():
    if k > random.random():
        mol1, mol2 = random.choices(mol_list, k = 2)                                    # 抽两个分子
        rct1, rct2 = random.choices(N_branch, k = 2)
        gro1, gro2 = group[mol1][rct1], group[mol2][rct2]                               # 该官能度
        if mol1 != mol2 and gro1 > 0 and gro2 > 0  and  network[mol1, mol2] != 1:  # 判断反应是否进行  and gro1 != gro2
            network[mol1, mol2], network[mol2, mol1] = 1, 1
            group[mol1][rct1], group[mol2][rct2] = -group[mol1][rct1], -group[mol2][rct2]


def cal(cycle):
    graph = nx.from_scipy_sparse_array(network)                                              # 构造nx图
    graph_component = nx.connected_components(graph)                                        # 获得graph_node所有子图
    tree_size = np.sort(np.array([graph.subgraph(i).number_of_nodes() for i in graph_component], dtype='int64'))   # n_node储存了所有子图的节点数目
    P = min((sum_gro - np.array([np.sum(group == i) for i in range(1, np.max(group) + 1)])) / sum_gro)
    P_theory = 1- (1 / ( 1 + k * cycle ))
    Mw = sum( tree_size ** 2 ) / N_mol                                                  # 重均分子量
    PMw = sum( tree_size[:-1] ** 2 ) / N_mol
    Mn = N_mol / len(tree_size)                                                        # 数均分子量
    # Mn_theory =
    # Mw_theory =
    index = int(cycle)                                                        # 数据的列数
    data.loc[index, column] = [cycle, P, Mn,  Mw, P_theory, PMw]
    if P < 0.5:
        data_size[str(cycle)] = np.append(tree_size, (N_mol - len(tree_size)) * [None])

# def draw():
#     plt.subplot(2, 2, 1)
#     y1 = np.array(data.loc[:, 'Mw']).T
#     y2 = np.array(data.loc[:, 'P']).T
#     y3 = np.array(data.loc[:, 'PMw']).T #/ np.array(data.loc[:, 'Mn'])
#     x = np.array(data.loc[:, 'cycle'])
#     plt.plot(y2, y1,   label='P')
#     plt.legend()
#     plt.subplot(2, 2, 2)
#     plt.plot(x, y2,  label='P_theory')
#     plt.legend()
#     plt.subplot(2, 2, 3)
#     plt.plot(y2, y3, label='Mn')
#     plt.legend()
#     plt.subplot(2, 2, 4)
#     y4 = np.array(data.loc[:, 'P'])/    np.array(data.loc[:, 'P_theory']  )
#     print(y4)
#     plt.plot(x, y4, label='')
#     plt.legend()
#     plt.show()


for i in tqdm(range(cycles)):
    react()
    if (i/N_mol) % interval == 0:
        cal(i/N_mol)
print(data_size)
print(data)
# draw()
data_size.to_csv("A:\data\交联k={0}.csv".format(k), index=0)
data.to_csv("A:\data\交联k={0}_size.csv".format(k), index=0)









