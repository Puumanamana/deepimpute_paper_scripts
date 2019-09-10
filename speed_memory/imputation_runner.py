import pandas as pd
from time import time

def timer(func):
    def wrapper(*args,**kwargs):
        t0 = time()
        func(*args,**kwargs)
        duration = time()-t0
        
        print("{}: {} s".format(func.__name__, duration))
        return duration
        
    return wrapper

class ImputationRunner:

    def __init__(self,path,method,ncores=20):
        self.data_path = path
        self.method = method
        self.ncores = ncores        
        self.writeback = True

    def load(self):
        data = pd.read_csv(self.data_path,index_col=0)

        return data

    def magic_impute(self,data):
        import magic
        
        model = magic.MAGIC(n_jobs=self.ncores)
        imputed = model.fit_transform(data.values)
        return pd.DataFrame(imputed)

    def di_impute(self,data):
        from deepimpute.multinet import MultiNet

        model = MultiNet(ncores=self.ncores)
        imputed = model.predict(data)
        return imputed

    def dca_impute(self,data):
        from dca.api import dca
        import scanpy as sc

        data = self.load()
        adata = sc.AnnData(data.values, obs=data.index, var=data.columns)
        imputed = dca(adata, threads=self.ncores).X
        return pd.DataFrame(imputed)

    @timer
    def run(self):
        data = self.load()

        if self.method.lower()=='dca':
            imputed = self.dca_impute(data)
        elif self.method.lower()=='deepimpute':
            imputed = self.di_impute(data)
        elif self.method.lower()=='magic':
            imputed = self.magic_impute(data)
        else:
            print('Unknown imputation method. Aborting.')
            exit(1)

        if self.writeback:
            imputed.to_csv('imputed.csv')
        

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--method', type=str, default='magic')
    parser.add_argument('--path', type=str)
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--trial', type=int, default=1)    
    args = parser.parse_args()

    runner = ImputationRunner(args.path,args.method,args.threads)
    dt = runner.run()

    with open("time_{}_{}.txt".format(args.method,args.trial),'w') as handle:
        handle.write('{} {} {}'.format(args.method,args.trial,dt))
