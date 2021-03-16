"""A module to describe information coming from ensemble averaging

"""


class Population():
    max_moments = 2

    @classmethod
    def load(cls, filename):
        from numpy import load
        try:
            pop_dict = load(filename)
        except ValueError:
            pop_dict = load(filename, allow_pickle=True, encoding="latin1")
        items = pop_dict.item()
        pop = cls(labels=items['feature_labels'])
        for d, w, l in zip(items['datas'], items['weights'],
                           items['sample_labels']):
            pop.append(d, weight=w, label=l)
        return pop

    def __init__(self, labels=None):
        self.cumsums = None
        self.datas = []
        self.sample_labels = []
        self.feature_labels = labels
        self.weights = []

    def append(self, data, weight=1.0, label=None):
        """
        Insert new data to the population.

        Args:
            data (array-like, pandas.DataFrame, dict): the data to append.
                In case of a dictionary, the data is internally converted into
                a dataframe.
            weight (float): the weight that the sample has in the population
            label (str): the label of the sample
        """
        from pandas import DataFrame
        if isinstance(data, dict):
            data_ = DataFrame(data)
            if self.feature_labels is None:
                self.feature_labels = data_.columns
        else:
            data_ = data
            if self.feature_labels is None:
                self.feature_labels = [str(i) for i in range(len(data))]

        sums = safe_moments(data, self.max_moments, weight)
        # [weight*data_, weight*data_**2, weight]
        self.weights.append(weight)
        if self.cumsums is None:
            self.cumsums = sums
        else:
            for isum in range(self.max_moments+1):
                self.cumsums[isum] = safe_add(self.cumsums[isum],
                                              safe_multiply_and_pow(sums[isum],
                                                                    1, 1))
                # self.cumsums[isum] += sums[isum]
        self.datas.append(data_)
        self.sample_labels.append(
            str(len(self.datas)-1) if label is None else label)

    def to_dataframe(self, sample_labels=None, feature_labels=None):
        """
        Convert the population into a pandas dataframe

        Args:
            sample_labels (list): overrides the sample labels
            feature_labels (list): overrides the feature labels

        Returns:
            pandas.Dataframe: the dataframe of the population
        """
        import pandas as pd
        if sample_labels is None:
            sample_labels = self.sample_labels
        if feature_labels is None:
            feature_labels = self.feature_labels
        dd = {sample: data for sample, data in zip(sample_labels, self.datas)}
        return pd.DataFrame.from_dict(dd, orient='index',
                                      columns=feature_labels)

    def to_dict(self):
        """
        Convert the population to a dictionary
        """
        return {att: getattr(self, att)
                for att in ['datas', 'feature_labels', 'sample_labels',
                            'weights']}

    def to_file(self, filename):
        """
        Dump the populations to a numpy file.

        Args:
            filename (str): the file for the dumping
        """
        from numpy import save
        pop_dict = self.to_dict()
        save(filename, pop_dict)

    def to_excel(self, filename=None, writer=None, sheets=['mean', 'std'],
                 prefix='', mappable=None):
        """
        Dump the population to an excel file.

        Args:
            filename (str): Name of the file to write the data to
            writer (pandas.ExcelWriter): the instance of the writer class.
                Useful to append more sheets to the same excel file
            sheets (list, str): list of the mehtods that will be written per
                each of the sheet of the file. It can be also a string, like
                "all" to write all the data of the population in a separate
                file
            prefix (str): prefix of the sheets to be employed
            mappable (func): a function to be applied to the data before
                writing them to the file
        Returns:
            pandas.ExcelWriter: the instance of the writer needed to create the
                file.
        """
        import pandas as pd
        if writer is None:
            writer = pd.ExcelWriter(filename, engine='xlsxwriter')
        if sheets == 'all':
            datas = self.datas
            shn = [prefix + '-' + str(lb) for lb in self.sample_labels]
        else:
            datas = [getattr(self, sht) for sht in sheets]
            shn = [prefix + '-' + sht for sht in sheets]
        for d, n in zip(datas, shn):
            df = d if mappable is None else mappable(d)
            df.to_excel(writer, sheet_name=n)
        return writer

    @property
    def full_weight(self):
        return self.cumsums[self.max_moments]

    @property
    def mean(self):
        return safe_multiply_and_pow(self.cumsums[0],
                                     1.0/self.full_weight, 1)

    @property
    def std(self):
        from numpy import sqrt, abs  # to avoid rundoff problems
        if len(self.datas) > 1:
            tmp1 = safe_sub(safe_multiply_and_pow(self.cumsums[1],
                                                  1.0/self.full_weight, 1),
                            safe_multiply_and_pow(self.mean, 1.0, 2))
            tmp1 = safe_unary_op(tmp1, abs)
            tmp1 = safe_unary_op(tmp1, sqrt)
            return tmp1
            # sqrt(abs(self.cumsums[1]/self.full_weight-self.mean**2))
        else:
            return 0.0


def remove_nan(df, transpose=False):
    df1 = df.dropna(how='all')
    dft = df1.transpose()
    df1 = dft.dropna(how='all')
    if transpose:
        return df1
    else:
        return df1.transpose()


def symmetrize_df(df1):
    """
    From a dataframe that should be asymmetrix matrix,
    construct the symmetrized dataframe
    """
    import numpy as np
    from pandas import DataFrame
    A = df1.values
    W = np.tril(A) + np.triu(A.T, 1)
    return DataFrame(W, columns=df1.columns, index=df1.index)


def reorder(dft, transpose=False):
    from BigDFT.IO import reorder_fragments
    dft1 = dft.reindex(index=reorder_fragments(dft.index))
    df1 = dft1.transpose()
    dft1 = df1.reindex(index=reorder_fragments(df1.index))
    if transpose:
        return dft1
    else:
        return dft1.transpose()


def clean_dataframe(df, symmetrize=True):
    """
    Symmetrize a dataframe and remove the NaN rows and columns

    Args:
       df (Dataframe)
       symmetrize (bool): symmetrize the dataframe if applicable
    Returns:
       Dataframe: the cleaned dataframe
    """
    dft = remove_nan(df, transpose=True)
    if symmetrize:
        dft = symmetrize_df(dft)
    return reorder(dft, transpose=True)


def stacked_dataframe(pop):
    """
    Construct a stacked dataframe with all the data of the population

    Warning:
        Weights are ignored, therefore the average value of such stacked
        dataframe may be different from the population mean.
    """
    from pandas import DataFrame as DF
    df = DF()
    for dd, name in zip(pop.datas, pop.sample_labels):
        df[name] = dd.stack()
    return df.T


def _convert_into_percent(weights):
    ntot = sum(weights)
    pcs = [(float(w)/ntot)*100.0 for w in weights]
    assert abs(sum(pcs) - 100.0) < 1.e-3
    return pcs


def weighted_dataframe(dfs, wgts):
    """
    Construct a single dataframe that include the provided weigths.
    Useful for all the situations where one wants to have a single view
    of a population which is weighted

    Args:
        dfs (list): list of dataframes to be weighted
        wgts(list):  weights to be included, should be of same length of dfs
    """
    from pandas import DataFrame
    from numpy import rint
    pcs = _convert_into_percent(wgts)
    df_tot = DataFrame()
    for df, pc in zip(dfs, pcs):
        for p in range(int(rint(pc))):
            df_tot = df_tot.append(df)
    return df_tot


def transverse_dataframe(population, target_features, target_samples):
    """
    Contruct a transverse dataframe that gathers the data from some target
    features and samples of a population. This may be useful to show
    the variability of some data for particular entries
    """
    df = stacked_dataframe(population)
    lookup = [(x, res) for x in target_samples for res in target_features]
    return df[lookup]


def concatenate_populations(populations, extra_labels={}):
    """
    Write a file containing the set of the populations we want to serialized

    Args:
         populations (dict): dictionary of the populations, labeled by the key
         extra_labels (dict): dictionary of the feature_labels
             and sample_labels of the populations
    Returns:
         pandas.Dataframe: the concatenated populations
    """
    from pandas import concat
    dfs = []
    for label, pop in populations.items():
        extra_lb = extra_labels.get(label, {})
        df = pop.to_dataframe(**extra_lb)
        df['label'] = label
        dfs.append(df)
    return concat(dfs)


def dump_populations(populations):
    """
    Dump a dictionary of populations in a archive

    Args:
         populations(dict): the dictionary of the populations

    Returns:
         set: the set of filenames (extension '.npy')
             which have been produced
    """
    files = []
    for label, pop in populations.items():
        filename = str(label) + '.npy'
        pop.to_file(filename)
        files.append(filename)
    return files


def safe_moments(data, order, weight):
    """
    Calculate a list of the moments of the distribution
    """
    moments = []
    for pw in range(order):
        moments.append(safe_multiply_and_pow(data, weight, pw+1))
    moments.append(weight)
    return moments


def safe_add(a, b):
    from pandas import DataFrame as DF
    try:
        return a - safe_multiply_and_pow(b, -1, 1)
    except Exception:
        return DF(a) + DF(b)


def safe_sub(a, b):
    return safe_add(a, safe_multiply_and_pow(b, -1, 1))


def safe_unary_op(data, op):
    """
    Apply the operation to the data in a dataframe-compatible way
    """
    from pandas import DataFrame as DF

    def _safe_op(x):
        try:
            return op(x)
        except Exception:
            return x

    try:
        return op(data)
    except Exception:
        try:
            return data.apply(_safe_op)
        except Exception:
            df = DF(data)
            return df.apply(op)


def safe_multiply_and_pow(data, a, pw):
    """
    Perform a multiplication and power expansion that is
    resilient to dataframes that contain sequences that are not
    floating-point number compatible
    """
    return safe_unary_op(data, _transform_pw(float(a), pw))


def _atimesx(a, x):
    if a == 1.0:
        return x
    else:
        return a*x


def _axpw(a, pw, x):
    if pw == 1:
        return _atimesx(a, x)
    else:
        return _atimesx(a, x**pw)


def _transform_pw(a=1.0, pw=1):
    from functools import partial
    return partial(_axpw, a, pw)


class ClusterGrammer():
    """
    A class that facilitates the use of the clustergrammer objects
    Args:
        df (pandas.DataFrame): the dataframe to represent the clustergrammer
    """
    def __init__(self, df):
        import clustergrammer_widget as cw
        self.df = df
        self.net = cw.Network(cw.clustergrammer_widget)
        self.net.load_df(df)

    def categorize(self, axis, cats):
        """
        Define categories to be displayed next to the axis elements

        Args:
            axis (str): 'row' or 'col'
            cats (dict): label: cats dictionary where cats is a dictionary
                where each key contains a list of fragments
        """
        recats = [{'title': tl, 'cats': c} for tl, c in cats.items()]
        self.net.add_cats(axis, recats)

    def represent_only(self, axis, elements):
        """
        Represent only the elements that are indicated on the given axis

        Args:
            axis (str): 'row' of 'col'
            elements (list): list of the elements to be represented on the
                given axis
        """
        from tempfile import NamedTemporaryFile as tmp
        # from os import remove
        import sys
        # redirect stdout
        old_stdout = sys.stdout
        ftmp = tmp(mode='w+')
        sys.stdout = ftmp
        self.net.filter_names(axis, elements)
        sys.stdout = old_stdout
        # remove(ftmp.name)

    def show(self):
        """
        Display the ClusterGrammer
        """
        self.net.cluster(enrichrgram=True)
        return self.net.widget()

    def publish(self, name):
        """
        Produce a link of the clustergrammer object that is
        in the public domain

        Args:
            name (str): name of the file to be created remotely
        """
        import requests
        filename = name+'.tsv'
        self.net.write_matrix_to_tsv(filename=filename)
        upload_url = 'http://amp.pharm.mssm.edu/clustergrammer/matrix_upload/'
        self.req = requests.post(upload_url,
                                 files={'file': open(filename, 'rb')})

    def publish_link(self):
        """
        Returns the URL of the published ClusterGrammer

        Returns:
            str: URL of the link
        """
        return self.req.text


class CData():
    """
    A class that facilitates the calls to data clustering algorithms.

    Args:
        datas (array-like): the data samples to perform the clustering to
    """
    def __init__(self, datas):
        self.datas = datas

    @property
    def lookup_nonnan(self):
        from numpy import isnan, where
        return where([not isnan(d) for d in self.datas[0]])

    @property
    def scaled_data(self):
        if not hasattr(self, '_scaled_data'):
            from sklearn.preprocessing import MinMaxScaler
            mms = MinMaxScaler()
            arr = self.datas.copy()
            mms.fit(arr)
            self._scaled_data = mms.transform(arr)
        return self._scaled_data

    def cluster_data(self, n_clusters):
        from sklearn.cluster import KMeans
        from sklearn.metrics import silhouette_samples, silhouette_score
        km = KMeans(n_clusters=n_clusters)
        km = km.fit(self.scaled_data)
        self.n_clusters = n_clusters
        self.cluster_labels = km.fit_predict(self.scaled_data)
        # Labeling the clusters
        self.cluster_centers = km.cluster_centers_
        # total distances of the clusters
        self.squared_distances = km.inertia_
        if self.n_clusters >= 2:
            # silhouette method (valid if nclusters >= 2)
            self.silhouette_avg = silhouette_score(self.scaled_data,
                                                   self.cluster_labels)
            # Compute the silhouette scores for each sample
            self.silhouette_values = silhouette_samples(self.scaled_data,
                                                        self.cluster_labels)
        return km

    def silhouette_plot(self, ax1=None):
        import numpy as np
        import matplotlib.cm as cm
        import matplotlib.pyplot as plt
        if ax1 is None:
            ax1 = plt.axes()
        # The 1st subplot is the silhouette plot
        # The silhouette coefficient can range from -1, 1
        ax1.set_xlim([-1, 1])
        # The (n_clusters+1)*10 is for inserting blank space between silhouette
        # plots of individual clusters, to demarcate them clearly.
        ax1.set_ylim([0, len(self.scaled_data) + (self.n_clusters + 1) * 10])
        y_lower = 10
        for i in range(self.n_clusters):
            # Aggregate the silhouette scores for samples belonging to
            # cluster i, and sort them
            ith_cluster_silhouette_values = \
                self.silhouette_values[self.cluster_labels == i]

            ith_cluster_silhouette_values.sort()

            size_cluster_i = ith_cluster_silhouette_values.shape[0]
            y_upper = y_lower + size_cluster_i

            color = cm.nipy_spectral(float(i) / self.n_clusters)
            ax1.fill_betweenx(np.arange(y_lower, y_upper),
                              0, ith_cluster_silhouette_values,
                              facecolor=color, edgecolor=color, alpha=0.7)

            # Label the silhouette plots with cluster numbers at the middle
            ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

            # Compute the new y_lower for next plot
            y_lower = y_upper + 10  # 10 for the 0 samples

        ax1.set_title("The silhouette plot for the various clusters.")
        ax1.set_xlabel("The silhouette coefficient values")
        ax1.set_ylabel("Cluster label")

        # The vertical line for average silhouette score of all the values
        ax1.axvline(x=self.silhouette_avg, color="red", linestyle="--")

        ax1.set_yticks([])  # Clear the yaxis labels / ticks
        ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])
        return ax1

    # 2nd Plot showing the actual clusters formed
    def plot_clusters(self, feature_one=None, feature_two=None, axis=None):
        import matplotlib.cm as cm
        import matplotlib.pyplot as plt
        if axis is None:
            axis = plt.axes()
        one = list(self.datas.columns).index(feature_one)
        two = list(self.datas.columns).index(feature_one)
        name_one = str(feature_one)
        name_two = str(feature_two)
        # name_one, one, name_two, two = self.largest_features_ids(feature_one,
        #                                                          feature_two)
        colors = cm.nipy_spectral(
                    self.cluster_labels.astype(float) / self.n_clusters)
        axis.scatter(self.scaled_data[:, one], self._scaled_data[:, two],
                     marker='.', s=30, lw=0, alpha=0.7, c=colors,
                     edgecolor='k')
        # Draw white circles at cluster centers
        axis.scatter(self.cluster_centers[:, one],
                     self.cluster_centers[:, two], marker='o', c="white",
                     alpha=1, s=200, edgecolor='k')
        for i, c in enumerate(self.cluster_centers):
            axis.scatter(c[one], c[two], marker='$%d$' % i, alpha=1, s=50,
                         edgecolor='k')
        axis.set_title("The visualization of the clustered data.")
        axis.set_xlabel("Feature space for "+name_one)
        axis.set_ylabel("Feature space for "+name_two)
        return axis

    def plot_cluster_info(self, feature_one=None, feature_two=None):
        import matplotlib.pyplot as plt
        # Create a subplot with 1 row and 2 columns
        fig, (ax1, ax2) = plt.subplots(1, 2)
        fig.set_size_inches(18, 7)
        if self.n_clusters >= 2:
            self.silhouette_plot(ax1)
        self.plot_clusters(feature_one, feature_two, ax2)
        return plt

    def largest_features_ids(self, feature_one, feature_two):
        # name_one, name_two = self.get_largest_variances_names()
        if feature_one is not None:
            name_one = feature_one
        if feature_two is not None:
            name_two = feature_two
        labels = [self.feature_labels[i] for i, t in enumerate(
                                                   self.lookup_nonnan[0]) if t]
        one = labels.index(name_one)
        two = labels.index(name_two)
        return [name_one, one, name_two, two]

    def get_largest_variances_names(self):
        import numpy as np
        ordered_variances = np.argsort(self.std).tolist()
        ordered_variances.reverse()
        name_one = None
        name_two = None
        for iname in ordered_variances:
            if np.isnan(self.mean[iname]):
                continue
            if name_one is None:
                name_one = self.feature_labels[iname]
                continue
            if name_two is None:
                name_two = self.feature_labels[iname]
                break
        return name_one, name_two
