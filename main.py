import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
import scipy as scp

def comb(n, k):
    return math.factorial(n) / (math.factorial(k) * math.factorial(n - k))
def case1_exact(N, P, column_names, index_names):
    # exact data's
    data = []
    for p in P:
        temp = []
        for n in N:
            q = 1 - p
            prob = 0
            exception = False
            for k in range(math.ceil(n / 2 - np.sqrt(n * p * q)), int(n / 2 + np.sqrt(n * p * q)) + 1):
                try:
                    prob += comb(n, k) * p ** k * q ** (n - k)
                except:
                    exception = True
                    temp.append("Can't calculate")
                    break
            if exception:
                continue
            temp.append(prob)
        data.append(temp)
    df = pd.DataFrame(data, columns=column_names, index=index_names)
    df.style.set_properties(**{'vertical-align': 'center',
                               'horizontal-align': 'center',
                               'text-align': 'center'})
    fig = plt.figure(figsize=(6.2, 1.2))
    ax = fig.add_subplot(111, frame_on=False)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    pd.plotting.table(ax, df, loc="center")
    plt.savefig("report/images/case1_exact.png")
def case1_Pois(N, P, column_names, index_names):
    data = []
    for p in P:
        temp = []
        for n in N:
            q = 1 - p
            prob = 0
            lam = n * p
            exception = False
            for k in range(math.ceil(n / 2 - np.sqrt(n * p * q)), int(n / 2 + np.sqrt(n * p * q)) + 1):
                try:
                    prob += np.exp(-lam) * lam ** k / math.factorial(k)
                except:
                    exception = True
                    temp.append("Can't calculate")
                    break
            if exception:
                continue
            temp.append(prob)
        data.append(temp)
    df = pd.DataFrame(data, columns=column_names, index=index_names)
    df.style.set_properties(**{'vertical-align': 'center',
                               'horizontal-align': 'center',
                               'text-align': 'center'})
    fig = plt.figure(figsize=(6.2, 1.2))
    ax = fig.add_subplot(111, frame_on=False)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    pd.plotting.table(ax, df, loc="center")
    plt.savefig("report/images/case1_pois.png")
def case1_LocalTh(N, P, column_names, index_names):
    data = []
    for p in P:
        temp = []
        for n in N:
            q = 1 - p
            prob = 0
            exception = False
            for k in range(math.ceil(n / 2 - np.sqrt(n * p * q)), int(n / 2 + np.sqrt(n * p * q)) + 1):
                try:
                    x = (k - n * p) / np.sqrt(n * p * q)
                    prob += np.exp(-x ** 2 / 2) / np.sqrt(2 * np.pi) / np.sqrt(n * p * q)
                except:
                    exception = True
                    temp.append("Can't calculate")
                    break
            if exception:
                continue
            temp.append(prob)
        data.append(temp)
    df = pd.DataFrame(data, columns=column_names, index=index_names)
    df.style.set_properties(**{'vertical-align': 'center',
                               'horizontal-align': 'center',
                               'text-align': 'center'})
    fig = plt.figure(figsize=(6.2, 1.2))
    ax = fig.add_subplot(111, frame_on=False)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    pd.plotting.table(ax, df, loc="center")
    plt.savefig("report/images/case1_local.png")
def case1_IntTh(N, P, column_names, index_names):
    data = []
    for p in P:
        temp = []
        for n in N:
            q = 1 - p
            try:
                k1 = n / 2 - np.sqrt(n * p * q)
                k2 = n / 2 + np.sqrt(n * p * q)
                a = (k1 - n * p) / np.sqrt(n * p * q)
                b = (k2 - n * p) / np.sqrt(n * p * q)
                F = lambda t: scp.integrate.quad(lambda x: np.exp(-x ** 2 / 2) / np.sqrt(2 * np.pi), -np.inf, t)[0]
                prob = F(b) - F(a)
                temp.append(prob)
            except:
                temp.append("Can't calculate")
        data.append(temp)
    df = pd.DataFrame(data, columns=column_names, index=index_names)
    df.style.set_properties(**{'vertical-align': 'center',
                               'horizontal-align': 'center',
                               'text-align': 'center'})
    fig = plt.figure(figsize=(6.2, 1.2))
    ax = fig.add_subplot(111, frame_on=False)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    pd.plotting.table(ax, df, loc="center")
    plt.savefig("report/images/case1_norm.png")
def case2_exact(N, P, column_names, index_names):
    data = []
    for p in P:
        temp = []
        for n in N:
            q = 1 - p
            prob = 0
            exception = False
            for k in range(0, 6):
                try:
                    prob += comb(n, k) * p ** k * q ** (n - k)
                except:
                    exception = True
                    temp.append("Can't calculate")
                    break
            if exception:
                continue
            temp.append(prob)
        data.append(temp)
    df = pd.DataFrame(data, columns=column_names, index=index_names)
    df.style.set_properties(**{'vertical-align': 'center',
                               'horizontal-align': 'center',
                               'text-align': 'center'})
    fig = plt.figure(figsize=(6.2, 1.2))
    ax = fig.add_subplot(111, frame_on=False)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    pd.plotting.table(ax, df, loc="center")
    plt.savefig("report/images/case2_exact.png")
def case2_Pois(N, P, column_names, index_names):
    data = []
    for p in P:
        temp = []
        for n in N:
            q = 1 - p
            prob = 0
            lam = n * p
            exception = False
            for k in range(0, 6):
                try:
                    prob += np.exp(-lam) * lam ** k / math.factorial(k)
                except:
                    exception = True
                    temp.append("Can't calculate")
                    break
            if exception:
                continue
            temp.append(prob)
        data.append(temp)
    df = pd.DataFrame(data, columns=column_names, index=index_names)
    df.style.set_properties(**{'vertical-align': 'center',
                               'horizontal-align': 'center',
                               'text-align': 'center'})
    fig = plt.figure(figsize=(6.2, 1.2))
    ax = fig.add_subplot(111, frame_on=False)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    pd.plotting.table(ax, df, loc="center")
    plt.savefig("report/images/case2_pois.png")
def case2_LocalTh(N, P, column_names, index_names):
    data = []
    for p in P:
        temp = []
        for n in N:
            q = 1 - p
            prob = 0
            exception = False
            for k in range(0, 6):
                try:
                    x = (k - n * p) / np.sqrt(n * p * q)
                    prob += np.exp(-x ** 2 / 2) / np.sqrt(2 * np.pi) / np.sqrt(n * p * q)
                except:
                    exception = True
                    temp.append("Can't calculate")
                    break
            if exception:
                continue
            temp.append(prob)
        data.append(temp)
    df = pd.DataFrame(data, columns=column_names, index=index_names)
    df.style.set_properties(**{'vertical-align': 'center',
                               'horizontal-align': 'center',
                               'text-align': 'center'})
    fig = plt.figure(figsize=(6.2, 1.2))
    ax = fig.add_subplot(111, frame_on=False)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    pd.plotting.table(ax, df, loc="center")
    plt.savefig("report/images/case2_local.png")
def case2_IntTh(N, P, column_names, index_names):
    data = []
    for p in P:
        temp = []
        for n in N:
            q = 1 - p
            try:
                k1 = 0
                k2 = 5
                a = (k1 - n * p) / np.sqrt(n * p * q)
                b = (k2 - n * p) / np.sqrt(n * p * q)
                F = lambda t: scp.integrate.quad(lambda x: np.exp(-x ** 2 / 2) / np.sqrt(2 * np.pi), -np.inf, t)[0]
                prob = F(b) - F(a)
                temp.append(prob)
            except:
                temp.append("Can't calculate")

        data.append(temp)
    df = pd.DataFrame(data, columns=column_names, index=index_names)
    df.style.set_properties(**{'vertical-align': 'center',
                               'horizontal-align': 'center',
                               'text-align': 'center'})
    fig = plt.figure(figsize=(6.2, 1.2))
    ax = fig.add_subplot(111, frame_on=False)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    pd.plotting.table(ax, df, loc="center")
    plt.savefig("report/images/case2_norm.png")
def case3_exact(N, P, column_names, index_names):
    data = []
    for p in P:
        temp = []
        for n in N:
            q = 1 - p
            k = int(np.floor((n+1) * p))
            try:
                prob = comb(n, k) * p ** k * q ** (n - k)
            except:
                temp.append("Can't calculate")
                continue
            temp.append(prob)
        data.append(temp)
    df = pd.DataFrame(data, columns=column_names, index=index_names)
    df.style.set_properties(**{'vertical-align': 'center',
                               'horizontal-align': 'center',
                               'text-align': 'center'})
    fig = plt.figure(figsize=(6.2, 1.2))
    ax = fig.add_subplot(111, frame_on=False)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    pd.plotting.table(ax, df, loc="center")
    plt.savefig("report/images/case3_exact.png")
def case3_Pois(N, P, column_names, index_names):
    data = []
    for p in P:
        temp = []
        for n in N:
            lam = n * p
            k = int(np.floor((n+1) * p))
            try:
                prob = np.exp(-lam) * lam ** k / math.factorial(k)
            except:
                temp.append("Can't calculate")
                continue
            temp.append(prob)
        data.append(temp)
    df = pd.DataFrame(data, columns=column_names, index=index_names)
    df.style.set_properties(**{'vertical-align': 'center',
                               'horizontal-align': 'center',
                               'text-align': 'center'})
    fig = plt.figure(figsize=(6.2, 1.2))
    ax = fig.add_subplot(111, frame_on=False)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    pd.plotting.table(ax, df, loc="center")
    plt.savefig("report/images/case3_pois.png")
def case3_LocalTh(N, P, column_names, index_names):
    data = []
    for p in P:
        temp = []
        for n in N:
            q = 1 - p
            k = np.floor((n + 1) * p)
            try:
                x = (k - n * p) / np.sqrt(n * p * q)
                prob = (np.exp(-x ** 2 / 2) /
                        np.sqrt(2 * np.pi) / np.sqrt(n * p * q))
            except:
                temp.append("Can't calculate")
                continue

            temp.append(prob)
        data.append(temp)
    df = pd.DataFrame(data, columns=column_names, index=index_names)
    df.style.set_properties(**{'vertical-align': 'center',
                               'horizontal-align': 'center',
                               'text-align': 'center'})
    fig = plt.figure(figsize=(6.2, 1.2))
    ax = fig.add_subplot(111, frame_on=False)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    pd.plotting.table(ax, df, loc="center")
    plt.savefig("report/images/case3_local.png")

def task5():
    N = [100, 1000, 10000]
    P = [0.001, 0.01, 0.1, 0.25, 0.5]
    index_names = ["p=0.001", "p=0.01", "p=0.1", "p=0.25", "p=0.5"]
    column_names = ["n=100", "n=1000", "n=10000"]
    # First point
    case1_exact(N, P, column_names, index_names)
    case1_Pois(N, P, column_names, index_names)
    case1_LocalTh(N, P, column_names, index_names)
    case1_IntTh(N, P, column_names, index_names)
    # Second point
    case2_exact(N, P, column_names, index_names)
    case2_Pois(N, P, column_names, index_names)
    case2_LocalTh(N, P, column_names, index_names)
    case2_IntTh(N, P, column_names, index_names)
    # Third point
    case3_exact(N, P, column_names, index_names)
    case3_Pois(N, P, column_names, index_names)
    case3_LocalTh(N, P, column_names, index_names)

def main():
    task5()
if __name__ == "__main__":
    main()
