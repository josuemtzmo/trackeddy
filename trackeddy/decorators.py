import pandas as pd
from functools import wraps


def eddy_identifier(f):
    """
    eddy_identifier Wrapper of eddy store and discarded functions

    Parameters
    ----------
    f : function
        Function that stores the data of identified and discarded eddies.
    """
    def wrapped(*args, **kwargs):

        # If argument exit is passed, then reset the wrap function
        if 'exit' in kwargs:
            wrapped.calls = 0
            wrapped.df_store = pd.DataFrame({'' : []})
            return wrapped
        
        # Get table to store eddy data.
        table = f(*args, id = wrapped.calls, **kwargs)
        # if table within wrapper is empty, assign the output of the function
        # else concatenate the tables together. 
        if wrapped.df_store.empty:
            wrapped.df_store = table
        else:
            wrapped.df_store = pd.concat((wrapped.df_store, table))
        # Counter in wrapper that works as the eddy identifier. 
        wrapped.calls += 1
        return wrapped.df_store
    # Definition of the counter of the function call and the table to store data.
    wrapped.calls = 0
    wrapped.df_store = pd.DataFrame({'' : []})
    return wrapped


def filter_data(func):
    """
    filter_data Wrapper to decide filtering to apply to the date, the options are 'space', 'time', and 'both'

    Parameters
    ----------
    func : function
        Function to apply the filter (_filter_data_)

    Returns
    -------
    function
        Wrapper of function _filter_data_
    """
    @wraps(func)
    def wrapper(self, *args, **kwargs):
        filtered = self.treat_nan(*args)
        if kwargs.get('filter') == 'both':
            filtered = func(self, filtered, 'space')
            self.data2track = func(self, filtered, 'time')
        elif  kwargs.get('filter') == 'space' or  kwargs.get('filter') == 'time':
            self.data2track = func(self, *args, kwargs.get('filter'))
        else:
            self.data2track = filtered

    return wrapper