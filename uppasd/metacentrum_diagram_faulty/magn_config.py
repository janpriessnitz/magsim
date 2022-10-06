import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as hcluster
from scipy.spatial import Delaunay
import itertools
from matplotlib import cm
from sklearn.cluster import DBSCAN
from pathlib import Path
import tarfile
import scipy.ndimage as ndimage

################################################################################

def get_config_dataframe(path, tag, num_ens=5):
    """
    Read magnetic configuration from a data folder with results.
    
    Parameters:
    path (str): path to the data folder
    tag (str): UppASD simulation tag (8 characters)
    num_ens (int): number of ensambles (default 5)
    
    Returns:
    df: dataframe with magnetic configuration (spatial coordinates + local magnetization components)
    """
    
    # path to data files
    data_dir = Path(path)
    if not data_dir.is_dir():
        print("Directory {:s} does not exists.".format(str(data_dir)))
        return -1
    
    # extract data from archive
    archive = data_dir / "results.tar.gz"  # archive with data
    if archive.is_file():
        try:
            with tarfile.open(archive, "r:gz") as tar:
                def is_within_directory(directory, target):
                    
                    abs_directory = os.path.abspath(directory)
                    abs_target = os.path.abspath(target)
                
                    prefix = os.path.commonprefix([abs_directory, abs_target])
                    
                    return prefix == abs_directory
                
                def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
                
                    for member in tar.getmembers():
                        member_path = os.path.join(path, member.name)
                        if not is_within_directory(path, member_path):
                            raise Exception("Attempted Path Traversal in Tar File")
                
                    tar.extractall(path, members, numeric_owner=numeric_owner) 
                    
                
                safe_extract(tar, path=data_dir)
        except tarfile.ReadError:
            print("File {:s} is corrupt.".format(str(archive)))
            return -1
    else:
        print("Archive {:s} does not exist.".format(str(archive)))
        return -1
    
    # paths to data files
    coord  = data_dir / "coord.{:s}.out".format(tag)    # file with coordinates
    config = data_dir / "restart.{:s}.out".format(tag)  # file with magnetic moments

    # read data to dataframes
    # coordinates
    if coord.is_file():
        df_coord  = pd.read_csv(coord, delim_whitespace=True, header=None).iloc[:,1:3]
    else:
        print("File {:s} does not exist.".format(str(coord)))
        return -1
    
    # magnetization components
    if config.is_file():
        df_config = pd.read_csv(config, delim_whitespace=True, header=None, 
                                skiprows=7).iloc[:,[1,4,5,6]]
    else:
        print("File {:s} does not exist.".format(str(config)))
        return -1
    
    # remove output files from the directory
    OutFiles = list(data_dir.glob("**/*.out"))  # list of all *.out files
    for file in OutFiles:
        try:
            os.remove(file)
        except OSError:
            print("Error while deleting file {:s}.".format(file))

    # concatenate dataframes
    df_coord_ens = pd.concat([df_coord]*num_ens, ignore_index=True)
    df = pd.concat([df_coord_ens, df_config], axis=1)
    df.columns = ["x", "y", "nens", "mx", "my", "mz"]
    
    # return dataframe sorted by coordinates   
    return(df.sort_values(by=["nens", "y", "x"]))

##############################

def store_config(df, filename, **kwargs):
    """
    Store magnetic configuration in dataframe to a hdf5 file including metadata.
    """
    
    store = pd.HDFStore(filename, "w")  # new hdf5 store
    store.put("data", df)  # save dataframe
    store.get_storer("data").attrs.metadata = kwargs  # save dataframe
    store.close()

##############################
    
def read_config(filename):
    """
    Read magnetic configuration and metadata from a hdf5 file.
    """
    
    store = pd.HDFStore(filename)  # open hdf5
    if "/data" in store.keys():
        data = store["data"]  # get dataframe
        metadata = store.get_storer("data").attrs.metadata  # get metadata
        store.close()
    
        return data, metadata
    else:
        return None, None

##############################
    
def add_metadata_to_config(filename, **kwargs):
    """
    Add metadata to existig configuration in a hdf5 file.
    """
    
    # read existing metadata
    store = pd.HDFStore(filename, "r")  # open hdf5
    metadata = store.get_storer("data").attrs.metadata  # get metadata
    store.close()
    
    # add new metadata
    for key in kwargs:
        metadata[key] = kwargs[key]
    
    # add to hdf5 store
    store = pd.HDFStore(filename, "a")  # new hdf5 store
    store.get_storer("data").attrs.metadata = metadata  # save dataframe
    store.close()

##############################

def get_config_matrices(df, sizex=200, sizey=200):
    """
    Returns list of 2D matrices of 3D vectors of magnetic moments.
    """
    
    return [df[df["nens"] == nens].iloc[:,3:6].to_numpy().reshape(sizex,sizey,3) for nens in df["nens"].unique()]

##############################

def get_magnetization(M):
    """
    Return total magnetization vector
    """

    return M.mean(axis=(0,1))

##############################

def get_local_topological_charge(Si, Sj, Sk):
    """
    Calculate local topological charge (in an elementary triangle on a square lattice).
    """
    
    # calculate dot products
    dij = np.dot(Si, Sj)
    djk = np.dot(Sj, Sk)
    dki = np.dot(Sk, Si)
    
    # local topological charge
    q = (1. + dij + djk + dki) / (np.sqrt(2. * (1. + dij) * (1. + djk) * (1. + dki)))
    sgn = np.sign(np.dot(Si, np.cross(Sj, Sk)))
    
    # treat inputs
    if q < -1:
        q = -1.
    elif q > 1:
        q = 1
    
    return 2. * sgn * np.arccos(q)
    
##############################

def get_topological_charge(M):
    """
    Calculate total topological charge on a 2D square lattice.
    """
    
    # shape of the matrix
    sx, sy, sv = M.shape
    
    # reshape matrix
    M = M.reshape(sx*sy, sv)
    
    # matrix of indices
    I = np.arange(sx*sy).reshape(sx, sy)
    
    # extend matrix of indices
    I = np.vstack((I, I[0]))
    I = np.hstack((I, I[:,0].reshape(-1,1)))
        
    # list of indices
    C1 = np.concatenate((I[:-1, :-1].ravel(), I[1:, :-1].ravel()))
    C2 = np.concatenate((I[:-1, 1:].ravel(), I[:-1, 1:].ravel()))
    C3 = np.concatenate((I[1:, :-1].ravel(), I[1:, 1:].ravel()))
    
    Q = 0.  # total topological charge
    for i, j, k in zip(C1, C2, C3):
        Q += get_local_topological_charge(M[i], M[j], M[k])  # local topological charge
    
    return Q / (4. * np.pi)

##############################

def get_spin_struct_factors(M):
    """
    Calculate logitudinal and transverse spin structur factors.
    """
    
    # system size
    sizex, sizey, sv = M.shape
    num_spins = sizex*sizey  # total number of spins
    
    # calculate fft for each component
    MFFT = [np.fft.fftshift(np.fft.fft2(M[:,:,i])) for i in range(3)]
    
    # out-of-plne spin structure factor
    Xoop = np.power(np.abs(MFFT[2]), 2) / num_spins
    
    # in-plane spin structure factor
    Xip = (np.power(np.abs(MFFT[0]), 2) + np.power(np.abs(MFFT[1]), 2)) / num_spins
    
    return Xoop, Xip

##############################

def get_peaks_number(X, eta=6.):
    """
    Get list of peaks values in the FFT array.
    """
    
    # apply maximum filter 
    XM = ndimage.maximum_filter(X, size=(2, 2))
    
    # find areas greater than given threshold
    thresh = XM.mean() + XM.std() * eta
    labels, num_labels = ndimage.label(XM > thresh)
    
    values = ndimage.measurements.maximum(X, labels=labels, index=np.arange(1, num_labels + 1))
        
    return values

##############################

def get_clusters_centers(Data, Clusters):
    """
    Return list of centers of the clusters in 2D layer.
    """

    IndClusters  = set(Clusters)     # indices of clusters
    num_clusters = len(IndClusters)  # number of clusters

    # create lists of cluster poinst
    ClustDict = {cind: [] for cind in IndClusters}  # dictionary with lists of clusters points

    for p, cind in zip(Data, Clusters):
        ClustDict[cind].append(p)

    # calculate centers of clusters
    Centers = []
    for Points in ClustDict.values():
        Centers.append(sum(Points) / len(Points))

    return np.array(Centers)

##############################

def get_skyrmion_positions(Mag, mz_dn=-.9, thresh=5., eta = .1):
    """
    Using hierarchical clustering estimate positions of skyrmions in the dataframe.
    
    Parameters:
    Mag: numpy (sizey, sizex, 3) array with magnetic configuration
    mz_df = -.97 (float): maximum value of magnetization z-component which is assumed as -1 ("down")
    thresh = 5 (float): threshold for the hierarchical clustering algorithm
    
    Returns:
    SkyrPos (list(numpy array)): list of skyrmion positions
    """
    
    # lattice parameters
    Mz = Mag[:,:,2]  # z-components
    sizex, sizey = Mz.shape  # system size
    
    # enlarge system
    eta = .1  # ratio of border to the lattice size
    nx, ny = int(eta * sizex), int(eta * sizey)  # width of the borders
    # bottom and top border
    Mz = np.vstack([Mz[-nx:,:], Mz, Mz[:nx,:]])
    # left and right border
    Mz = np.hstack([Mz[:,-ny:], Mz, Mz[:,:ny]])
    
    # positions with spin-down
    Iy, Ix = np.where(Mz < mz_dn)  # indices
    DataDn = np.transpose([Ix - nx, Iy - ny])  # coordinates
    
    if len(DataDn) == 0:
        return []

    # clustering
    Clusters = hcluster.fclusterdata(DataDn, thresh, criterion="distance", metric="euclidean", method="average")
    
    if len(set(Clusters)) == 0:
        return []

    # calculate positions of the skyrmion centers in the enlarged supercell
    SkyrPos = get_clusters_centers(DataDn, Clusters)
    
    # select point in the original frame
    SkyrPosIn = []
    for p in SkyrPos:
        if (0 <= p[0] < sizex and 0 <= p[1] < sizey):
            SkyrPosIn.append(p)
            
    return np.array(SkyrPos)

##############################

def plot_mz_map(Mag, xlim=None, ylim=None, Positions=None, Triangulation=None, 
                title=None, cmap="bwr", pcol="tab:green", lcol="tab:green", figname=None):
    """
    Plot magnetic configuration (x, y, mz)
    
    Parameters:
    df (DataFrame): dataframe with magnetic configuration
    Attr (Dictionary): dictionary of attributes
    xlim (float, float): x subrange
    ylim (float, float): y subrange
    Positions: list of points positions
    Triangulation: triangulation between the points
    title (str): plot title
    cmap (str): color map
    pcol (str): points color
    lcol (str): lines color
    figname (str): name of output figure
        
    Returns:
    plot
    """
        
    fig, ax = plt.subplots(figsize=(12,12))
    
    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)
    if title:
        ax.set_title(title, size=18)
   
    # plot map of the magnetization z-component
    I = plt.imshow(Mag[:,:,2], cmap=cmap, vmin=-1., vmax=1.)

    fig.colorbar(I)

    # plot triangulation
    if Triangulation is not None:
        ax.triplot(*np.transpose(Positions), Triangulation.simplices, color=lcol)

    # plot centers of skyrmions
    if (Positions is not None) and (len(Positions) != 0):
        A = np.transpose(Positions)
        ax.scatter(A[0], A[1], marker="d", s=50, c=pcol)

    # export figure
    if figname:
        plt.savefig(figname)

    plt.show()
