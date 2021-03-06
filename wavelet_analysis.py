import pywt
from operator import itemgetter
from scipy.spatial import distance
from hcluster import pdist, linkage, dendrogram, fcluster
import pylab
from numpy import *

def partition(lst, n):
	division = len(lst) / float(n)
	return [ lst[int(round(division * i)): int(round(division * (i + 1)))] for i in xrange(n) ]

def wavelet_avg(Y,X=None,reshape=False,plot_xy=False,zero_thresh=1*10**-7,Names=None,Title=None):
	if not X:
		X = range(Y.shape[0])
	w = pywt.Wavelet('sym2')
	L = pywt.dwt_max_level(Y.shape[0],w)
	#L = pywt.swt_max_level(Y.shape[0])
	Zavg = zeros((L,len(X)))
	Zavgabs = zeros((L,len(X)))
	if reshape:
		e = [0,L,0,L]
	else:
		e = [X[0],X[-1],0,L]
	sp = 1
	if not Names:
		Names = ["Series "+str(y+1) for y in range(Y.shape[1])]
	if Title:
		pylab.suptitle(Title)
	if plot_xy:
		tp = str(Y.shape[1]+3)
		pylab.subplot(int(''.join([tp,'1',str(sp)])))
		P = pylab.plot(X,Y)
		pylab.legend(P,Names,prop={'size':6})
		pylab.xlim([0,Y.shape[0]-1])
		pylab.title("Relative Abundance Time Series")
		sp += 1
	else:
		tp = str(Y.shape[1]+2)
	for y in range(Y.shape[1]):
		C = pywt.wavedec(Y[:,y],w,level=L)
		#C = pywt.swt(Y[:,y],w,level=L)
		#C = ['dummy'] + [c[1] for c in C]
		Z = wavelet_matrix(C,L,len(X))
		pylab.subplot(int(''.join([tp,'1',str(sp)])))
		pylab.imshow(abs(Z),extent=e)
		pylab.title(Names[y])
		sp += 1
		Zavg += Z
		Zavgabs += abs(Z)
	Zavg /= Y.shape[1]
	Zavg[Zavg<zero_thresh] = 0
	pylab.subplot(int(''.join([tp,'1',str(sp)])))
	pylab.imshow(abs(Zavg),extent=e)
	pylab.title("Decomposition Avg")
	sp += 1
	Zavgabs /= Y.shape[1]
	Zavgabs[Zavgabs<zero_thresh] = 0
	pylab.subplot(int(''.join([tp,'1',str(sp)])))
	pylab.imshow(abs(Zavgabs),extent=e)
	pylab.title("abs(Decomposition) Sum")
	sp += 1
	pylab.show()
	return Zavg

def wavelet_matrix(C,r,c):
	Z = zeros((r,c))
	C1 = C[1:]
	for i in range(r):
		yp = partition(range(c),len(C1[i]))
		for yj in range(len(yp)):
			for j in yp[yj]:
				Z[i,j] += C1[i][yj]
	return Z

def nonlinear_ab(n,b1=0.02,b2=0.1,p=1):
	X = linspace(1,n,n)
	Y = zeros((len(X),2))
	Y[:p,0] = random.rand(1,p)
	Y[:p,1] = Y[0,0]
	for x in range(p,len(X)):
		Y[x,0] = Y[x-p,0]*(3.8 - 3.8*Y[x-p,0] - b1*Y[x-p,1])
		Y[x,1] = Y[x-p,1]*(3.5 - 3.5*Y[x-p,1] - b1*Y[x-p,0])
	return X,Y

def doppler(n,m,bmin=0,bmax=1,reflections=1):
	x = linspace(0,1,n/(reflections+1))
	y = sqrt(x*(1-x))*sin((2.1*pi)/(x+.05))
	y += -y.min()
	y /= y.max()
	if reflections%2 == 1:
		y = array((list(y) + list(y)[::-1])*((reflections+1)/2))
	elif reflections%2 == 0:
		y = list(y) + (list(y)[::-1] + list(y))*(reflections/2)
		y += list(y[:n-len(y)])
		y = array(y)
	ym = zeros((n,m))
	ym[:,0] = y
	for j in range(1,m):
		b = random.rand()*(bmax - bmin) + bmin
		ym[:,j] = ym[:,0]*b + random.rand(1,n)*(1-b)
	return range(n),ym

def rand_abundances(y0,n):
	Y = zeros((len(y0),n+1))
	Y[:,0] = y0
	for i in range(Y.shape[0]):
		r = random.rand(1,n)
		Y[i,1:] = r*(1-Y[i,0])/r.sum()
	return Y

def normalize_to_abundances(Y,amin=0,amax=1):
	for i in range(Y.shape[0]):
		Y[i,:] = Y[i,:]/Y[i,:].sum()
	return Y

def c_dists(Y,use_swt=True,level_weights=False):
	w = pywt.Wavelet('sym2')
	if use_swt:
		L = pywt.swt_max_level(Y.shape[0])
		C = [pywt.swt(Y[:,i],w,level=L) for i in range(Y.shape[1])]
		C = [[list(reshape(l[0],-1)) + list(reshape(l[1],-1)) for l in c] for c in C]
	else:
		L = pywt.dwt_max_level(Y.shape[0],w)
		C = [pywt.wavedec(Y[:,i],w,level=L) for i in range(Y.shape[1])]
	if level_weights:
		if use_swt:
			raise NameError('No level weights with SWT')
		Wc = [1. for x in range(1,L+1)]
		D = zeros((len(C),len(C)))
		for i in range(len(C)):
			for j in range(i+1,len(C)):
				d = sum([distance.cosine(C[i][x],C[j][x])*Wc[x] for x in range(L)])/sum(Wc)
				D[i,j] = d
				D[j,i] = d
		return D
	else:
		Cn = []
		for c in C:
			cn = []
			for l in c:
				cn += list(l)
			Cn.append(cn)
		return abs(pdist(Cn,'cosine'))

def wavelet_clusters(Y,ct=0.5,weights=False,return_clusters=False,swt=False):
	if weights:
		D = abs(c_dists(Y,level_weights=True,use_swt=False))
		Dr = []
		for i in range(D.shape[0]-1):
			Dr += list(D[i,i+1:])
	else:
		Dr = c_dists(Y,use_swt=swt)
	if return_clusters:
		L = linkage(Dr,method='single',metric='cosine')
		C = fcluster(L,ct,criterion='distance')
		return cluster_sets(C)
	plot_clusters(Dr,ct)

def time_series_clusters(Y,ct=0.5,return_clusters=False):
	D = pdist(transpose(Y),'correlation')
	D = abs(D)
	if return_clusters:
		L = linkage(D,method='single',metric='cosine')
		C = fcluster(L,ct,criterion='distance')
		return cluster_sets(C)
	plot_clusters(D,ct)

def plot_clusters(Dr,ct):
	L = linkage(Dr,method='single',metric='cosine')
	dendrogram(L,color_threshold=ct)
	pylab.show()

def create_population(types,n,rand_periods=1):
	Y = []
	for t in types:
		if t[0] == 'nonlinear':
			x,y = nonlinear_ab(n,b1=random.rand()*.2,b2=random.rand()*.75,p=random.randint(1,1+rand_periods))
		elif t[0] == 'doppler':
			x,y = doppler(n,t[1],reflections=random.randint(0,rand_periods))
		elif t[0] == 'rand':
			y = random.rand(n,t[1])
		if len(Y) > 0:
			Y = concatenate((Y,y),axis=1)
		else:
			Y = y
	return normalize_to_abundances(Y*random.power(.5,(1,Y.shape[1]))[:,newaxis][0])

def check_clusters(result_sets,answer_sets):
	s = 0
	matched_sets = []
	for rs in result_sets:
		for ans in answer_sets:
			if len(rs & ans)/float(len(rs)) > .7:
				s += len(rs)
				matched_sets.append(rs)
				break
	return s,matched_sets

def cluster_sets(cluster_array):
	cluster_array = sorted(enumerate(cluster_array),key=itemgetter(1))
	rs = [[cluster_array[0][0]]]
	for i in range(1,len(cluster_array)):
		if cluster_array[i][1] != cluster_array[i-1][1]:
			rs.append([])
		rs[-1].append(cluster_array[i][0])
	return [set(l) for l in rs if len(l)>1]

def iter_rand_clusters(types,m,n=10,r=3,swt=True):
	A = {'nonlinear': [],'doppler': []}
	i = 0
	for t in types:
		if t[0] != 'rand':
			A[t[0]] += range(i,i+t[1])
			i += t[1]
	A0 = [set(A['nonlinear'] + A['doppler'])]
	A = [set(v) for v in A.values() if len(v)>1]
	ta = sum([len(s) for s in A])
	print ta
	tsum = 0
	wsum = 0
	if swt:
		wct = 0.1
	else:
		wct = 0.08
	for i in range(n):
		Y = create_population(types,m,rand_periods=r)
		t = time_series_clusters(Y,ct=0.5,return_clusters=True)
		t0 = check_clusters(t,A)
		t = wavelet_clusters(Y,ct=wct,return_clusters=True,swt=swt)
		t1 = check_clusters(t,A0)
		if (t0[0] < ta) and (t1[0] < ta):
			print t0[0],t1[0]
			tsum += t0[0]
			wsum += t1[0]
	print tsum,wsum