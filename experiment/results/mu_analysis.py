import numpy as np
import pandas as pd
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")


# load R results
def load_data():
    net_1 = np.loadtxt(Path(__file__).parent / "coeff_1.txt", delimiter=",", dtype=float)
    net_2 = np.loadtxt(Path(__file__).parent / "coeff_2.txt", delimiter=",", dtype=float)
    return net_1, net_2


def get_dropout(net: np.ndarray):
    return np.sum(net == 0)/len(net.flatten())


# load network
def load_ref_data():
    genes = np.loadtxt(Path(__file__).parent / "genes.txt", delimiter="\t", dtype=str)
    genes_dic = dict(zip(genes, np.arange(len(genes))))
    net_ref = pd.read_csv(Path(__file__).parent.parent / "data/Curated/GSD/GSD-2000-1/refNetwork.csv", 
                        )
    net_ref["ref"] = -1
    net_ref["ref"][net_ref["Type"] == "+"] = 1
    net_ref["pair_idx"] = pd.concat([net_ref["Gene1"].map(genes_dic), 
                                    net_ref["Gene2"].map(genes_dic)], 
                                    axis=1).apply(
                                    lambda x: (x[0], x[1]), axis=1)
    # adjacent
    adj_ref = np.zeros((len(genes), len(genes)), dtype=int)
    for r, i in zip(net_ref["ref"], net_ref["pair_idx"]):
        adj_ref[i] = r
    # print("dropout of adj:", np.sum(adj_ref==0)/len(adj_ref.flatten()))		
    return net_ref, adj_ref

    
#muSuSiE
# map adj to ref network
def net_mapping(net:np.ndarray, # weighted adj
                net_ref:pd.DataFrame, # ref network
                mask=None):
    idx = (net_ref["pair_idx"]).tolist()
    net_res = np.array([1 if net[i]>=0 else -1 for i in idx])
    if mask is not None:
        acc = np.sum(net_ref["ref"][mask] == net_res[mask]) / len(net_ref[mask])
    else:
        acc = np.sum(net_ref["ref"] == net_res) / len(net_ref)
    
    return acc # acc on network only


def adj_stacking(net: np.ndarray,
                adj_ref: np.ndarray,
                ):
    net_t = net.copy()
    net_t[net>0] = 1
    net_t[net<0] = -1
    net_t = net_t.astype(int)
    return np.sum(net_t == adj_ref)/len(net_t.flatten()) # acc on adj


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-L', '--L_max', type = int, default = 10)
    parser.add_argument('--file', type = str, default = 'sum_filter')
    args = parser.parse_args()

    net_1, net_2 = load_data()
    drop_1, drop_2 = get_dropout(net_1), get_dropout(net_2)
    net_ref, adj_ref = load_ref_data()
    acc_1, acc_2 = net_mapping(net_1, net_ref), net_mapping(net_2, net_ref)
    acc_1_adj, acc_2_adj = adj_stacking(net_1, adj_ref), adj_stacking(net_2, adj_ref)
    try:
        f = open(Path(__file__).parent / f'{args.file}.txt', 'r')
    except IOError:
        f = open(Path(__file__).parent / f'{args.file}.txt', 'w')
        f.writelines("L_max\tdropout\tacc_network_only\tacc_adj\n")
    finally:  
        f = open(Path(__file__).parent / f'{args.file}.txt', 'a')
        f.writelines(f"{args.L_max}, ({drop_1:.2f}, {drop_2:.2f}),"
                     f"({acc_1:.2f}, {acc_2:.2f}),"
                     f"({acc_1_adj:.2f}, {acc_2_adj:.2f})\n")						
        f.close()
        