"""Calculation datasets.

This module deals with the handling of series of calculations.
Classes and functions of this module are meant to simplify the approach to ensemble calculations
with the code, and to deal with parallel executions of multiple instances of the code.

"""

from Calculators import Runner

def name_from_id(id):
    """Hash the id into a run name
    Construct the name of the run from the id dictionary

    Args:
        id (dict): id associated to the run

    Returns:
       str: name of the run associated to the dictionary ``id``
    """
    keys=id.keys()
    keys.sort()
    name=''
    for k in keys:
        name += k+':'+str(id[k])+','
    return name.rstrip(',')


class Dataset(Runner):
    """A set of calculations.

    Such class contains the various instances of a set of calculations with the code.
    The different calculations are labelled by parameter values and information that might then be
    retrieved for inspection and plotting.

    Args:
      label (str): The label of the dataset. It will be needed to identify the instance for example
          in plot titles or in the running directory.
      run_dir (str): path of the directory where the runs will be performed.
      input (dict): Inputfile to be used for the runs as default,
             can be overridden by the specific inputs of the run

    """
    def __init__(self,label='BigDFT dataset',run_dir='runs',**kwargs):
        """
        Set the dataset ready for appending new runs
        """
        from copy import deepcopy
        from futile.Utils import make_dict
        newkwargs=deepcopy(kwargs)
        Runner.__init__(self,label=label,run_dir=run_dir,**newkwargs)
        self.runs=[]
        """List of the runs which have to be treated by the dataset these runs contain the input parameter to be passed to the various runners.
        """
        self.calculators=[]
        """
        Calculators which will be used by the run method, useful to gather the inputs in the case of a multiple run.
        """

        self.results={}
        """
        Set of the results of each of the runs. The set is not ordered as the runs may be executed asynchronously.
        """

        self.ids=[]
        """
        List of run ids, to be used in order to classify and fetch the results
        """

        self.names=[]
        """
        List of run names, needed for distinguishing the logfiles and input files.
        Eah name should be unique to correctly identify a run.
        """

        self._post_processing_function=None

    def append_run(self,id,runner,**kwargs):
        """Add a run into the dataset.

        Append to the list of runs to be performed the corresponding runner and the arguments which are associated to it.

        Args:
          id (dict): the id of the run, useful to identify the run in the dataset. It has to be a dictionary as it may contain
             different keyword. For example a run might be classified as ``id = {'hgrid':0.35, 'crmult': 5}``.
          runner (Runner): the runner class to which the remaining keyword arguments will be passed at the input.

        Raises:
          ValueError: if the provided id is identical to another previously appended run.

        Todo:
           include id in the runs spcification

        """
        from copy import deepcopy
        name=name_from_id(id)
        if name in self.names:
            raise ValueError('The run id',name,' is already provided, modify the run id.')
        self.names.append(name)
        #create the input file for the run, combining run_dict and input
        inp_to_append=deepcopy(self._global_options)
        inp_to_append.update(deepcopy(kwargs))
        #get the number of this run
        irun=len(self.runs)
        #append it to the runs list
        self.runs.append(inp_to_append)
        #append id and name
        self.ids.append(id)
        #search if the calculator already exists
        found = False
        for calc in self.calculators:
            if calc['calc'] == runner:
                calc['runs'].append(irun)
                found=True
                break
        if not found:
            self.calculators.append({'calc': runner, 'runs':[irun]})

    def process_run(self):
        """
        Run the dataset, by performing explicit run of each of the item of the runs_list.
        """
        for c in self.calculators:
            calc=c['calc']
            #we must here differentiate between a taskgroup run and a separate run
            for r in c['runs']:
                inp=self.runs[r]
                name=self.names[r]
                self.results[r]=calc.run(name=name,**inp)
        return {}

    def set_postprocessing_function(self,func):
        """Set the callback of run.

        Calls the function ``func`` after having performed the appended runs.

        Args:
           func (func): function that process the `inputs` `results` and returns the
               value of the `run` method of the dataset.
               The function is called as ``func(self)``.

        """
        self._post_processing_function=func

    def post_processing(self,**kwargs):
        """
        Calls the Dataset function with the results of the runs as arguments
        """
        if self._post_processing_function is not None:
            return self._post_processing_function(self)
        else:
            return self.results

    def fetch_results(self,id=None,attribute=None):
        """Retrieve some attribute from some of the results.

        Selects out of the results the objects which have in their ``id``
        at least the dictionary specified as input. May return an attribute
        of each result if needed.

        Args:
           id (dict): dictionary of the retrieved id. Return a list of the runs that
               have the ``id`` argument inside the provided ``id`` in the order provided by :py:meth:`append_run`.
           attribute (str): if present, provide the attribute of each of the results instead of the result object

        Example:
           >>> study=Dataset()
           >>> study.append_run(id={'cr': 3},input={'dft':{'rmult':[3,8]}})
           >>> study.append_run(id={'cr': 4},input={'dft':{'rmult':[4,8]}})
           >>> study.append_run(id={'cr': 3, 'h': 0.5},
           >>>                  input={'dft':{'hgrids': 0.5, 'rmult':[4,8]}})
           >>> #append other runs if needed
           >>> study.run()  #run the calculations
           >>> # returns a list of the energies of first and the third result in this example
           >>> data=study.fetch_results(id={'cr': 3},attribute='energy')
        """
        name='' if id is None else name_from_id(id)
        data=[]
        for irun,n in enumerate(self.names):
            if name not in n: continue
            r=self.results[irun]
            data.append(r if attribute is None else getattr(r,attribute))
        return data



def combine_datasets(*args):
    """
    Define a new instance of the dataset class that should provide
    as a result a list of the runs of the datasets
    """
    full=Dataset(label='combined_dataset')
    #append the runs or each dataset
    for dt in args:
        for irun,runs in enumerate(dt.runs):
            calc=dt.get_runner(irun)
            id,dt.get_id(irun)
            full.append_run(id,calc,**runs)

    full.set_postprocessing_function(_combined_postprocessing_functions)


def _combined_postprocessing_functions(runs,results,**kwargs):
    pass
