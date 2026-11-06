"""
exact_posterior_4cs.py -- exact posterior for every Condition 4 examinee.

For each examinee, on the items administered and the responses given: the posterior
under the N(0, 1) prior on a grid, (a) with the bank's item parameters and (b) with
the likelihood averaged over each item's 500 parameter draws. Prints RMSE, mean
posterior SD and coverage for the engine and for both exact posteriors.

Use:  python3 exact_posterior_4cs.py      (from the folder that holds results/)
"""
import csv,json,math,statistics as st
import numpy as np
R='results/2026-10-30_162500_brm4cs'
D=json.load(open(R+'/testdef_brm4cs.json'))
bank={it['itm_id']:it for it in D['item_bank'] if it['status']=='OP'}
grid=np.arange(-6,6.0001,0.01); prior=np.exp(-grid**2/2)
def probs(a,b,d,model,x):

    if model=='GPCM':
        steps=np.concatenate([np.zeros((a.shape[0],1)),d],axis=1) if d is not None else None
        m=steps.shape[1]
        z=np.cumsum(a[:,None,None]*(grid[None,None,:]-b[:,None,None]+steps[:,:,None]),axis=1)
        z=np.concatenate([np.zeros((a.shape[0],1,grid.size)),z[:,1:,:]],axis=1) if False else z

        z=z - (a[:,None,None]*(grid[None,None,:]-b[:,None,None]))
        z=z - z.max(axis=1,keepdims=True); e=np.exp(z); p=e/e.sum(axis=1,keepdims=True)
        return p[:,x,:]
    p1=1/(1+np.exp(-a[:,None]*(grid[None,:]-b[:,None])))
    return p1 if x==1 else 1-p1
def item_arrays(it,draws):
    if draws:
        B=it['bayesian']; a=np.array(B['a_vector']); b=np.array(B['b_vector'])
        d=np.array(B['d_vector']).T if it['model']=='GPCM' and B.get('d_vector') and B['d_vector'][0] else (np.tile(it['d'],(a.size,1)) if it['model']=='GPCM' else None)
    else:
        a=np.array([it['a']]); b=np.array([it['b']]); d=np.array([it['d']]) if it['model']=='GPCM' else None
    return a,b,d
S={r['student']:r for r in csv.DictReader(open(R+'/scores.csv')) if r['domain_name']=='Overall'}
adm={}
for r in csv.DictReader(open(R+'/items_long.csv')):
    if int(r['itemScore'])>=0: adm.setdefault(r['simulee'],[]).append((r['itemID'],int(r['itemScore'])))
cache={}
def lik(it_id,x,draws):
    k=(it_id,x,draws)
    if k not in cache:
        a,b,d=item_arrays(bank[it_id],draws); cache[k]=probs(a,b,d,bank[it_id]['model'],x).mean(axis=0)
    return cache[k]
ids=sorted(S,key=int)
res=[]
for s in ids:
    lm=np.log(prior.copy()); lmg=np.log(prior.copy())
    for it,x in adm[s]:
        lm+=np.log(lik(it,x,False)); lmg+=np.log(lik(it,x,True))
    out=[]
    for l in (lm,lmg):
        w=np.exp(l-l.max()); w/=w.sum(); m=(w*grid).sum(); sd=math.sqrt((w*(grid-m)**2).sum()); out+= [m,sd]
    res.append([float(S[s]['true_theta']),float(S[s]['theta']),float(S[s]['se'])]+out)
res=np.array(res)
t=res[:,0]
for nm,(mc,sc) in (('engine',(1,2)),('exact modal',(3,4)),('exact marginal',(5,6))):
    e=res[:,mc]-t; print(nm,'RMSE %.4f meanSE %.4f cov %.4f'%(math.sqrt((e**2).mean()),res[:,sc].mean(),(abs(e)<=1.96*res[:,sc]).mean()))
r=res[:,6]/res[:,4]; print('marginal/modal SD ratio: mean %.4f, range %.4f-%.4f; by theta<-1.5 %.4f mid %.4f >1.5 %.4f'%(r.mean(),r.min(),r.max(),r[t<-1.5].mean(),r[(t>=-1.5)&(t<=1.5)].mean(),r[t>1.5].mean()))
print('max |mean diff| marginal-modal %.3f'%abs(res[:,5]-res[:,3]).max())
