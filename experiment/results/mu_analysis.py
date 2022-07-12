import numpy as np
import pandas as pd
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")


# load R results
def load_data(paths: list):
    if not isinstance(paths, list):
        raise TypeError("input a list of .txt file paths")
    return [np.loadtxt(Path(__file__).parent / path, delimiter=",", dtype=float) for path in paths]


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
    if mask is not None: # pos/neg
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
    import glob
    parser = argparse.ArgumentParser()
    parser.add_argument('-L', '--L_max', type = int, default = 10) # only for write
    parser.add_argument('--file_out', type = str, default = 'sum_filter')
    parser.add_argument('--file_in_dir', type = str, default = 'multi_coeff')
    args = parser.parse_args()

    files = glob.glob(f"{args.file_in_dir}/*.txt")
    nets = load_data(files)
    drops = [get_dropout(net) for net in nets]
    net_ref, adj_ref = load_ref_data()
    # get accuracies on network
    accs = [net_mapping(net, net_ref) for net in nets]
    mask = net_ref["ref"] == 1
    accs_pos = [net_mapping(net, net_ref, mask) for net in nets]
    mask = net_ref["ref"] == -1
    accs_neg = [net_mapping(net, net_ref, mask) for net in nets]
    # accuracy on adj
    accs_adj = [adj_stacking(net, adj_ref) for net in nets]

    drop, acc, acc_pos, acc_neg, acc_adj = np.mean([drops, accs, accs_pos, accs_neg, accs_adj], axis=1)
    try:
        f = open(Path(__file__).parent / f'{args.file_out}.txt', 'r')
    except IOError:
        f = open(Path(__file__).parent / f'{args.file_out}.txt', 'w')
        f.writelines("L_max,dropout,acc_network,acc_network_pos,acc_network_neg,acc_adj\n")
    finally:  
        f = open(Path(__file__).parent / f'{args.file_out}.txt', 'a')
        f.writelines(f"{args.L_max},"
                     f"{drop:.2f},{acc:.2f},{acc_pos:.2f},{acc_neg:.2f},{acc_adj:.2f}\n")
                    #  f"{drop_1:.2f},{drop_2:.2f},"
                    #  f"{acc_1:.2f},{acc_2:.2f},"
                    # dd f"{acc_pos_1:.2f},{acc_pos_2:.2f},"
                    #  f"{acc_neg_1:.2f},{acc_neg_2:.2f},"
                    #  f"{acc_adj_1:.2f},{acc_adj_2:.2f}\n")						
        f.close()

    # python mu_analysis.py -L 10 --file_out sum_filter
        