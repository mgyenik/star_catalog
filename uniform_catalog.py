import math
import pandas
import numpy as np
import pylab
from pylab import *

# N is the number of points that should be uniformly distributed on the spiral.
N = 3000

INITIAL_MAG = 4.0
DEGREES = .017453 # 1 degree in radians
RAD_THRESH = 3*DEGREES

# uses HYG-Database from github
datafile = '/home/m/git/HYG-Database/hygxyz.csv'

def spiral_angles(k):
	m = math.sqrt(N*math.pi)
	theta = np.array([np.arccos(1 - float((2*j-1))/k) for j in xrange(1, k)])
	phi = m*theta
	return theta, phi

def point_plot_3d(x, y, z):
	from mpl_toolkits.mplot3d import Axes3D
	fig = pylab.figure()
	ax = fig.gca(projection='3d')
	#ax.plot(x,y,z)
	ax.scatter(x,y,z)
	show()

# Get unit vectors given spiral angles.
def xyz_from_sphere(th, ph):
	x = np.multiply(np.sin(th), np.cos(ph))
	y = np.multiply(np.sin(th), np.sin(ph))
	z = np.cos(th)
	return np.array((x, y, z)).T

def hipparcos_unit_vectors():
	df = pandas.DataFrame.from_csv(datafile)
	mag = np.sqrt(df['X']**2 + df['Y']**2 + df['Z']**2)
	df['x'] = df['X']/mag
	df['y'] = df['Y']/mag
	df['z'] = df['Z']/mag
	#return df[['x', 'y', 'z', 'Mag']]
	return df

# Let s be the full star catalog, and u the set of nearly uniformly distributed unit vectors.
def build_uniform_catalog(S, U):
	# c is the star catalog, which we initialize using all stars above some
	# arbitrary magnitude threshold (we want the brightest stars to be part
	# of our catalog no matter what!). Remember though, "magnitude" is like
	# golf scores; the lower the magnitude, the brighter the star.
	c = S[S.Mag <= INITIAL_MAG]

	# We also remove the stars we just added from s, transferring them to the catalog,
	# so they are not added twice when we process more stars below.
	S = S[S.Mag > INITIAL_MAG]

	# Now for each of the unit vectors (that point to one patch of sky), we
	# add all stars to our catalog that are within some small angular threshold
	# into the catalog, and remove them from further consideration, thus collecting
	# a few fainter stars within the patch to make sure it is covered.
	for u in U:
		v = S[np.arccos(np.dot(np.array((S.x, S.y, S.z)).T, u)) < RAD_THRESH]

		# Select the brightest star (lowest magnitude) in the set
		# of stars near the unit vector, and add it to our catalog!
		c = c.append(v[v.Mag == v.Mag.min()], ignore_index=True)
		v = v[v.Mag != v.Mag.min()]

	return c

def point_plot_catalog(c):
	x, y, z = c['x'].values, c['y'].values, c['z'].values
	point_plot_3d(x,y,z)
	show()

if __name__ == '__main__':
	print 'Finding hipparcos unit vectors...'
	hipparcos = hipparcos_unit_vectors()
	print 'Done finding hipparcos unit vectors: '
	#point_plot_3d(*xyz_from_sphere(*spiral_angles(N)))
	print 'Computing nearly uniform unit vectors for ', N, ' total patches...'
	t, p = spiral_angles(N)
	#x = np.arange(1,N)
	#plot(x, t)
	#plot(x, p)
	#show()
	u = xyz_from_sphere(t, p)
	print 'Done computing ', len(u), ' total unit vectors.'
	print 'Building uniform catalog...'
	catalog = build_uniform_catalog(hipparcos, u)
	print 'Done building uniform catalog!'
	print 'Total number of stars: ', catalog['x'].count()
	print 'Mean magnitude:        ', catalog.Mag.mean()

	magthresh = hipparcos[hipparcos.Mag <= 6.0]
	mt_ra, mt_dec = magthresh['RA'].values, magthresh['Dec'].values
	c_ra, c_dec = catalog['RA'].values, catalog['Dec'].values

	subplot(2,1,1)
	scatter(mt_ra, mt_dec)
	ylabel('Declination (degrees)')
	title('Magnitude 6.0 thresholded catalog ({} stars)'.format(len(mt_ra)))

	subplot(2,1,2)
	scatter(c_ra, c_dec)
	xlabel('Right Ascention (hours)')
	ylabel('Declination (degrees)')
	title('Nearly uniform catalog ({} stars)'.format(len(c_ra)))
	show()

