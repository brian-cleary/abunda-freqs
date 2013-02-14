from numpy import *
from time import strptime
from collections import defaultdict
from operator import itemgetter
import ftplib

def download_stuff(file_type,domain='ftp.metagenomics.anl.gov',project_dir='~/projects/93/',output_path='data/',sample_subdir='/processed/'):
	ftp = ftplib.FTP(domain)
	ftp.login()
	ftp.cwd(project_dir)
	sample_dirs = ftp.nlst()
	for sd in sample_dirs:
		try:
			ftp.cwd(project_dir+sd+sample_subdir)
			f = open(output_path+sd+file_type,'wb')
			ftp.retrbinary("RETR "+file_type,f.write)
			f.close()
		except Exception,err:
			print str(err),sd

def parse_time_point(fp,min_count=0):
	f = open(fp)
	Lines = f.readlines()
	f.close()
	abundance_counts = {}
	for line in Lines:
		l = line.strip().split('\t')
		if int(l[1]) > min_count:
			abundance_counts[l[0]] = int(l[1])
	return abundance_counts

meta_fields = ['common_sample_site','collection_date','host_individual']
def parse_metadata(fp):
	f = open(fp)
	Lines = f.readlines()
	f.close()
	Meta = {}
	for line in Lines:
		l = line.strip().split()
		if l[1] in meta_fields:
			if l[1] == 'collection_date':
				Meta[l[1]] = strptime(l[2],'%m-%d-%Y')
			else:
				Meta[l[1]] = l[2]
	return Meta
	

def relative_abundances(counts):
	total_count = float(sum(counts.values()))
	return dict([(k,v/total_count) for k,v in counts.items()])

def abundance_time_series(count_file_paths,min_appearances=1):
	RA = []
	M = []
	all_ids = []
	for fp in count_file_paths:
		ac = parse_time_point(fp)
		ra = relative_abundances(ac)
		all_ids += ra.keys()
		try:
			M.append(parse_metadata(fp[:fp.index('.')]+'.3.meta'))
			RA.append(ra)
		except:
			print fp
	Host_Sites = defaultdict(list)
	for i in range(len(M)):
		Host_Sites[(M[i]['common_sample_site'],M[i]['host_individual'])].append((M[i]['collection_date'],i))
	for k,v in Host_Sites.items():
		Host_Sites[k] = sorted(v,key=itemgetter(0))
	all_ids = [_ for _ in enumerate(set(all_ids))]
	A = dict([(a[1],a[0]) for a in all_ids])
	all_ids = dict(all_ids)
	Host_Site_Time_Series = {}
	for k,v in Host_Sites.items():
		T = zeros((len(v),len(all_ids)))
		for i in range(len(v)):
			ra = RA[v[i][1]]
			for a in ra.items():
				T[i,A[a[0]]] = a[1]
		appearances = [i[0] for i in enumerate((T != 0).sum(0)) if i[1]>min_appearances]
		ai = [all_ids[i] for i in appearances]
		T = T[T.any(1),:]
		Host_Site_Time_Series[k] = (dict(enumerate(ai)),T[:,appearances])
	return Host_Site_Time_Series

def smooth_matrix(Y,max_change=10):
	Ys = Y
	for j in range(Y.shape[1]):
		Ys[:,j] = smooth_series(Y[:,j],max_change)
	return Ys

def smooth_series(y,mc):
	y_smooth = y
	ys = sorted(enumerate(y),key=itemgetter(1),reverse=True)
	ynzs = sorted(enumerate(y[y.nonzero()]),key=itemgetter(1))
	median_nz = ynzs[len(ynzs)/2][1]
	for yi in ys:
		if yi[0] == 0:
			neighbor_avg = y[1]
		elif yi[0] == len(y)-1:
			neighbor_avg = y[-2]
		else:
			neighbor_avg = (y[yi[0]-1] + y[yi[0]+1])/2.
		if neighbor_avg == 0:
			y_smooth[yi[0]] = median_nz
		elif yi[1]/neighbor_avg >mc:
			y_smooth[yi[0]] = neighbor_avg*mc
		else:
			break
	return y_smooth