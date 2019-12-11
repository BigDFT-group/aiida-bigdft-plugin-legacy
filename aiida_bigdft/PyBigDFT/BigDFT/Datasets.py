"""Calculation datasets.

This module deals with the handling of series of calculations.
Classes and functions of this module are meant to simplify the approach to
ensemble calculations with the code, and to deal with parallel executions of
multiple instances of the code.

"""

from BigDFT.Calculators import Runner


def name_from_id(id):
    """Hash the id into a run name
    Construct the name of the run from the id dictionary

    Args:
        id (dict): id associated to the run

    Returns:
       str: name of the run associated to the dictionary ``id``
    """
    keys = list(id.keys())
    keys.sort()
    name = ''
    for k in keys:
        name += k + ':' + str(id[k]) + ','
    return name.rstrip(',')


def names_from_id(id):
    """Hash the id into a list of run names to search with the function id_in_names
       and add the separator ',' to have the proper value of a key
       (to avoid 0.3 in 0.39)
    """
    if id is None:
        return ['']
    else:
        return [name_from_id({k: v})+',' for k, v in id.items()]


class Dataset(Runner):
    """A set of calculations.

    Such class contains the various instances of a set of calculations with the
    code.
    The different calculations are labelled by parameter values and information
    that might then be retrieved for inspection and plotting.

    Args:
      label (str): The label of the dataset. It will be needed to identify the
          instance for example in plot titles or in the running directory.
      run_dir (str): path of the directory where the runs will be performed.
      input (dict): Inputfile to be used for the runs as default,
             can be overridden by the specific inputs of the run

    """

    def __init__(self, label='BigDFT dataset', run_dir='runs', **kwargs):
        """
        Set the dataset ready for appending new runs
        """
        from copy import deepcopy
        newkwargs = deepcopy(kwargs)
        Runner.__init__(self, label=label, run_dir=run_dir, **newkwargs)
        self.runs = []
        """
        List of the runs which have to be treated by the dataset these runs
        contain the input parameter to be passed to the various runners.
        """
        self.calculators = []
        """
        Calculators which will be used by the run method, useful to gather the
        inputs in the case of a multiple run.
        """

        self.results = {}
        """
        Set of the results of each of the runs. The set is not ordered as the
        runs may be executed asynchronously.
        """

        self.ids = []
        """
        List of run ids, to be used in order to classify and fetch the results
        """

        self.names = []
        """
        List of run names, needed for distinguishing the logfiles and
        input files. Each name should be unique to correctly identify a run.
        """

        self._post_processing_function = None

    def append_run(self, id, runner, **kwargs):
        """Add a run into the dataset.

        Append to the list of runs to be performed the corresponding runner and
           the arguments which are associated to it.

        Args:
          id (dict): the id of the run, useful to identify the run in the
             dataset. It has to be a dictionary as it may contain
             different keyword. For example a run might be classified as
             ``id = {'hgrid':0.35, 'crmult': 5}``.
          runner (Runner): the runner class to which the remaining keyword
             arguments will be passed at the input.

        Raises:
          ValueError: if the provided id is identical to another previously
             appended run.

        Todo:
           include id in the runs specification

        """
        from copy import deepcopy
        name = name_from_id(id)
        if name in self.names:
            raise ValueError('The run id', name,
                             ' is already provided, modify the run id.')
        self.names.append(name)
        # create the input file for the run, combining run_dict and input
        inp_to_append = deepcopy(self._global_options)
        inp_to_append.update(deepcopy(kwargs))
        # get the number of this run
        irun = len(self.runs)
        # append it to the runs list
        self.runs.append(inp_to_append)
        # append id and name
        self.ids.append(id)
        # search if the calculator already exists
        found = False
        for calc in self.calculators:
            if calc['calc'] == runner:
                calc['runs'].append(irun)
                found = True
                break
        if not found:
            self.calculators.append({'calc': runner, 'runs': [irun]})

    def process_run(self):
        """
        Run the dataset, by performing explicit run of each of the item of the
           runs_list.
        """
        self._run_the_calculations()
        return {}

    def _run_the_calculations(self, selection=None):
        for c in self.calculators:
            calc = c['calc']
            # we must here differentiate between a taskgroup run and a
            # separate run
            for r in c['runs']:
                if selection is not None and r not in selection:
                    continue
                inp = self.runs[r]
                name = self.names[r]
                self.results[r] = calc.run(name=name, **inp)

    def set_postprocessing_function(self, func):
        """Set the callback of run.

        Calls the function ``func`` after having performed the appended runs.

        Args:
           func (func): function that process the `inputs` `results` and
               returns the value of the `run` method of the dataset.
               The function is called as ``func(self)``.

        """
        self._post_processing_function = func

    def post_processing(self, **kwargs):
        """
        Calls the Dataset function with the results of the runs as arguments
        """
        if self._post_processing_function is not None:
            return self._post_processing_function(self)
        else:
            return self.results

    def fetch_results(self, id=None, attribute=None, run_if_not_present=True):
        """Retrieve some attribute from some of the results.

        Selects out of the results the objects which have in their ``id``
        at least the dictionary specified as input. May return an attribute
        of each result if needed.

        Args:
           id (dict): dictionary of the retrieved id. Return a list of the runs
               that have the ``id`` argument inside the provided ``id`` in the
               order provided by :py:meth:`append_run`.
           attribute (str): if present, provide the attribute of each of the
               results instead of the result object
           run_if_not_present (bool): If the run has not yet been performed
               in the dataset then perform it.

        Example:
           >>> study=Dataset()
           >>> study.append_run(id={'cr': 3},input={'dft':{'rmult':[3,8]}})
           >>> study.append_run(id={'cr': 4},input={'dft':{'rmult':[4,8]}})
           >>> study.append_run(id={'cr': 3, 'h': 0.5},
           >>>                  input={'dft':{'hgrids': 0.5, 'rmult':[4,8]}})
           >>> #append other runs if needed
           >>> #run the calculations (optional if run_if_not_present=True)
           >>> study.run()
           >>> # returns a list of the energies of first and the third result
           >>> # in this example
           >>> data=study.fetch_results(id={'cr': 3},attribute='energy')
        """
        names = names_from_id(id)
        fetch_indices = []
        selection_to_run = []
        for irun, name in enumerate(self.names):
            # add the separator ',' to have the proper value of a key
            # (to avoid 0.3 in 0.39)
            if not all([(n in name+',') for n in names]):
                continue
            if run_if_not_present and irun not in self.results:
                selection_to_run.append(irun)
            fetch_indices.append(irun)
        if len(selection_to_run) > 0:
            self._run_the_calculations(selection=selection_to_run)

        data = []
        for irun in fetch_indices:
            r = self.results[irun]
            data.append(r if attribute is None else getattr(r, attribute))
        return data

    def seek_convergence(self, rtol=1.e-5, atol=1.e-8, selection=None,
                         **kwargs):
        """
        Search for the first result of the dataset which matches the provided
        tolerance parameter. The results are in dataset order
        (provided by the :py:meth:`append_run` method) if `selection` is not
        specified.
        Employs the numpy :py:meth:`allclose` method for comparison.

        Args:
          rtol (float): relative tolerance parameter
          atol (float): absolute tolerance parameter
          selection (list): list of the id of the runs in which to perform the
               convergence search. Each id should be unique in the dataset.
          **kwargs: arguments to be passed to the :py:meth:`fetch_results`
               method.

        Returns:
          id,result (tuple): the id of the last run which matches the
                convergence, together with the result, if convergence is
                reached.

        Raises:
           LookupError: if the parameter for convergence were not found.
               The dataset has to be enriched or the convergence parameters
               loosened.
        """
        from numpy import allclose
        from futile.Utils import write
        to_get = self.ids if selection is None else selection

        id_ref = to_get[0]
        write('Fetching results for id "', id_ref, '"')
        ref = self.fetch_results(id=id_ref, **kwargs)
        ref = ref[0]
        for id in to_get[1:]:
            write('Fetching results for id "', id, '"')
            val = self.fetch_results(id=id, **kwargs)
            val = val[0]
            if allclose(ref, val, rtol=rtol, atol=atol):
                res = self.fetch_results(id=id_ref)
                label = self.get_global_option('label')
                write('Convergence reached in Dataset "' +
                      label+'" for id "', id_ref, '"')
                return (id_ref, res[0])
            ref = val
            id_ref = id
        raise LookupError('Convergence not reached, enlarge the dataset'
                          ' or change tolerance values')

    def get_times(self):
        from os import path as p
        from futile import Time as T
        time_files = {}
        run_dir = self.get_global_option('run_dir')
        time_files = [p.join(run_dir,
                             self.fetch_results(id=self.ids[c],
                                                attribute='data_directory')[0],
                             'time-' + self.names[c] + '.yaml'
                             ) for c in self.results
                      ]
        return T.TimeData(*time_files)

    def wait(self):
        from IPython.display import display, clear_output
        from aiida.orm import load_node
        running = len(self.results)
        while(running != 0):

            import time
            time.sleep(1)
            running = len(self.results)
            for c in self.results:
                pk = self.results[c]['node'].pk
                node = load_node(pk)
                if(node.is_finished):
                    running -= 1
                    # print(node.is_finished_ok)
            clear_output(wait=True)
            display(str(running)+" processes still running")

    def get_logfiles(self):
        logfiles = {}
        for c in self.results:
            try:
                logfiles[c] = self.calculators[0]['calc'].get_logs(
                    self.results[c]['node'].pk, self.names[c])
            except ValueError:
                logfiles[c] = self.results[c]
                print("no logfile for " + str(c))
        return logfiles


def combine_datasets(*args):
    """
    Define a new instance of the dataset class that should provide
    as a result a list of the runs of the datasets
    """
    full = Dataset(label='combined_dataset')
    # append the runs or each dataset
    for dt in args:
        for irun, runs in enumerate(dt.runs):
            calc = dt.get_runner(irun)
            id, dt.get_id(irun)
            full.append_run(id, calc, **runs)

    full.set_postprocessing_function(_combined_postprocessing_functions)


def _combined_postprocessing_functions(runs, results, **kwargs):
    pass
