import pandas as pd
import pytz
from datetime import datetime
from datetime import timedelta
import os
import json
import clifpy
import numpy as np

#file paths
work_dir = os.path.abspath('..')
output_folder = os.path.join(work_dir,'output')

#Config
with open(os.path.join(work_dir,'config','config.json'), 'r') as file:
    config = json.load(file)
config['output_folder'] = output_folder

#My time_zone
my_tz = config['time_zone']

#Time Zone things
def ensure_datetime(s: pd.Series) -> pd.Series:
    """Parse to datetime; convert UTC into my_tz"""
    s = pd.to_datetime(s, errors='coerce')
    try:
        # Convert to EST if timezone-aware
        if s.dt.tz is not None:
            return s.dt.tz_convert(my_tz)
        else:
            return s.dt.tz_localize(my_tz)
    except (TypeError, AttributeError):
        return s

def convert_datetime_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Applies ensure_datetime_est to all datetime columns."""
    df = df.copy()
    for col in df.columns:
        # Only operate on datetime-like columns
        if pd.api.types.is_datetime64_any_dtype(df[col]):
            df[col] = ensure_datetime(df[col])
        
    return df

def load_data(data_set: str, name: str, folder: str = '', type: str = None) -> pd.DataFrame:
    """Loads a table from the CLIF/MIMIC/mobilization/output data sets. Input: string with table name. Output: DataFrame. Automatically converts to my_tz."""
    
    def _load_path(path: str, file_type: str) -> pd.DataFrame:
        if file_type == 'parquet':
            return pd.read_parquet(path)
        elif file_type in ('csv', 'csv.gz'):
            return pd.read_csv(path)
        else:
            raise ValueError('File type needs to be parquet, csv or csv.gz.')

    def _build_path(file_type: str) -> str:
        if folder:
            return os.path.join(config[data_set], folder, f"{name}.{file_type}")
        return os.path.join(config[data_set], f"{name}.{file_type}")

    if type is not None:
        path = _build_path(type)
        df = _load_path(path, type)
    else:
        for candidate in ('parquet', 'csv', 'csv.gz'):
            path = _build_path(candidate)
            if os.path.exists(path):
                df = _load_path(path, candidate)
                break
        else:
            raise FileNotFoundError(
                f"No file found for '{name}' in '{data_set}' "
                f"(tried parquet, csv, csv.gz)."
            )

    df = convert_datetime_columns(df)
    return df

def missing_summary(in_df:pd.DataFrame,f_name:str=None) -> pd.DataFrame:
    def missing(x): return round(np.mean(x.isna())*100,2)
    miss_df = in_df.agg(missing).reset_index()
    miss_df.columns = ['column_name','missing']
    miss_df.iloc[0] = {'column_name':'size', 'missing':in_df.shape[0]}
    if f_name:
        path = os.path.join(output_folder,'final',f'{f_name}_missing.csv')
        miss_df.to_csv(path)
        return path
    else:
        return miss_df

def table_summary(in_df:pd.DataFrame) -> pd.DataFrame:
    def missing(x): return round(np.mean(x.isna())*100,2)
    def Q1(x): return x.quantile(0.25)
    def Q3(x): return x.quantile(0.25)
    numeric_df = in_df.select_dtypes(include=['number','bool'])
    for col in numeric_df.columns:
        if pd.api.types.is_bool_dtype(numeric_df[col]):
            numeric_df[col] = numeric_df[col].astype(int)
    sum_df = numeric_df.agg(['min',Q1,'median','mean',Q3,'max', missing]).reset_index().T
    cols = sum_df.iloc[0]
    sum_df = sum_df.iloc[1:]
    sum_df.columns = cols
    return sum_df

#time_bins object
class time_bins:

    def __init__(self, in_name:str=None, in_df:pd.DataFrame=None, in_eb:pd.DataFrame=None, resort:bool=True):
        '''
        INIT
        inputs:
        in_name: name of the file where the data frame is initially stored. It will look for it in "output/intermediate", assumes parquet file.
        in_df: a DataFrame with time_bins already built in.
        in_eb: a DataFrame of encounter blocks to build the time_bins.
        NOTES: If given a DF rather than building it, it must macth expected structure.
        '''
        self.bin_array = np.arange(0, config['time_end'] + config['time_bin_size'], config['time_bin_size'])
        self.bins_df = pd.DataFrame({
            'bin_start': self.bin_array[:-1],
            'bin_end': self.bin_array[1:]
        })
        self.bins_df['time_bin'] = self.bins_df['bin_start'].astype(str)
        
        if in_name:
            self.df = load_data('output_folder',in_name, folder='intermediate')
            extra_bins = sum(~self.df['time_bin'].isin(self.bins_df['time_bin'])) > 0
        elif in_df is not None:
            self.df = in_df.copy()
            extra_bins = sum(~self.df['time_bin'].isin(self.bins_df['time_bin'])) > 0
        elif in_eb is not None:
            self.df = pd.merge(in_eb, self.bins_df, how='cross').reset_index()
            extra_bins = False
        else:
            raise ValueError("Input to initialize hourly_blocks must include in_path=file name or in_df=DataFrame.")

        self.required_cols =  ['encounter_block','bin_start','bin_end','time_bin']
        missing = [col for col in self.required_cols if col not in self.df.columns]
        if missing:
            raise ValueError(f"Missing columns: {missing}")

        if resort:
            self.df = self.df.sort_values(by=['encounter_block','bin_start'])
        if extra_bins:
            print("WARNING: There are time_bin values in the provided DataFrame which are unexpected.")
    
    def classify_time_bin(self, t_diff:pd.Series) -> pd.Series:
        '''
        INPUT
        t_diff = Pandas DataSeries with a time_delta time element or a numerical element assumed to be hours
        OUTPUT
        DataSeries with time_bins[time_bin] as the value.
        '''
        #If time_delta format convert to hours.
        if pd.api.types.is_timedelta64_dtype(t_diff):
            t_diff = t_diff.dt.total_seconds()/3600
    
        return pd.cut(t_diff, bins=self.bin_array, labels=self.bins_df['time_bin'].tolist(), right=True, include_lowest=True)

    def bin_sort_fill(self,
                      val_col:str,
                      fill_with):
        '''
        Sorts time_bins and fill them in with specified request.
        INPUT
        val_col = val_col string name in the self.df.
        fill_with = either forward fill or 1
        '''
        
        if fill_with == 'ffill':
            self.df = self.df.sort_values(by=['encounter_block','bin_start']) #A bit repetitive but hard to have do this other way
            self.df[val_col] = self.df.groupby('encounter_block')[val_col].ffill()
        else:
            # NB: chained `self.df[col].fillna(..., inplace=True)` is a silent no-op
            # under pandas 3.0 Copy-on-Write. Assign the result back instead.
            self.df[val_col] = self.df[val_col].fillna(fill_with)
    
    def gather_time_bins(self,
                         input_df:pd.DataFrame,
                         val_col:str,
                         agg_func:str = 'mean',
                         fill_with = np.nan):
        '''
        INPUTS:
        df = DataFrame required to have:
            Either a column called 'time_diff' which is a date time difference value or any number (assumed hours) or 'time_bin'.
            val_col, the data itself to be aggregated.
            encounter_block
            time_bin = unique time string can be created with function above.
        agg = aggregation function to aggregate
        '''
        #Time bins
        if 'time_bin' not in input_df.columns:
            input_df['time_bin'] = self.classify_time_bin(input_df['time_diff'])
    
        #Copy data frame
        in_df = input_df[['encounter_block','time_bin',val_col]].copy()
        #Remove time windows outside of the time bins
        in_df = in_df[in_df['time_bin'].notna()]
        
        #Name of new variable
        new_name = f"{val_col}_{agg_func}"
    
        #Flag for agg_func = flag or all
        or_flag = agg_func == 'flag'
        if or_flag:
            agg_func = 'max'
            in_df[val_col] = in_df[val_col].astype(bool)
        and_flag = agg_func == 'all'
        if and_flag:
            agg_func = 'mean'
            in_df[val_col] = in_df[val_col].astype(bool)
    
        #Aggregation step
        print(f'Aggregation step for {new_name}')
        agg_df = in_df.groupby(by=['encounter_block','time_bin'])[val_col].agg(agg_func).reset_index()

        #Rename column
        agg_df.rename(columns={val_col:new_name},inplace=True)
        
        #Spcial situation for "or" and "and" to convert to boolean and back to int to 0 or 1.
        if or_flag or and_flag:
            agg_df[new_name] = agg_df[new_name] == 1
            agg_df[new_name] = agg_df[new_name].astype(int)
        
        #Merge onto time bins
        print(f'Merging step for {new_name}')
        self.df = self.df.merge(agg_df, on=['encounter_block','time_bin'], how='left', sort=False)

        #Sort and fill
        print(f'Sort/Fill step for {new_name}')
        self.bin_sort_fill(new_name, fill_with)

    def add_event(self,
                  input_df:pd.DataFrame,
                  new_name:str):
        '''
        INPUT
        in_df: A data frame which requires 'encounter_block', and either 'time_diff' or 'time_bin'.
        new_name: a string with a name for the new columns.
        NOTE: It will be back fill with 0 and forward fill with 1.
        '''
        #Time bins
        if ~('time_bin' in input_df.columns):
            input_df['time_bin'] = self.classify_time_bin(input_df['time_diff'])
    
        #Copy data frame
        in_df = input_df[['encounter_block','time_bin']].copy().drop_duplicates()
        #Remove time windows outside of the time bins
        in_df = in_df[in_df['time_bin'].notna()]
        #Dummy true variable
        in_df[new_name] = 1
        #Put into bin architecture
        self.df = pd.merge(self.df, in_df, on=['encounter_block','time_bin'], how='left', sort=False)
        #Fill ahead with 1 if present.
        self.bin_sort_fill(new_name, 'ffill')
        #Fill everything else with zero
        self.bin_sort_fill(new_name, 0)
        in_df[new_name] = in_df[new_name].astype(bool)

    def remove_based_on_censor(self, censor_col:str, keep_first:bool = False):
        '''
        Removes rows from the time bins based on some censoring variable.
        INPUT
        censor_col: String with the name of the column in time bins to use for this.
        keep_first: Boolean to either keep the first columns to be censored (ie for death where we keep the bin in which they died but no the subsequent bins).
        '''
        start_N = self.df.shape[0]
        if keep_first:
            rm = self.df.duplicated(subset=['encounter_block',censor_col], keep='first') & self.df[censor_col]
        else:
            rm = self.df[censor_col]
        self.df = self.df[~rm]
        end_N = self.df.shape[0]
        print(f"Removed {start_N - end_N} out of {start_N} from time_bin.df based on censor {censor_col} with keep_first = {keep_first}")
    
    def table_summary(self) -> pd.DataFrame:
        return table_summary(self.df)
        
    def save(self, suffix:str=''):
        path = os.path.join(output_folder,'intermediate',f'time_bin{suffix}.parquet')
        self.df.to_parquet(path)
        path = os.path.join(output_folder, 'final',f'time_bin{suffix}_summary.csv')
        _summary = self.table_summary()
        _summary.to_csv(path)

#Older aggregation by time function for a single time period
def aggregate_by_time(in_df:pd.DataFrame, val_col:str, min_time:int = 0, max_time:int = 24, agg_func:str = 'mean'):
    '''
    INPUTS:
    df = DataFrame required to have:
        column called 'time_diff' which is a date time difference value or any number (assumed hours)
        val_col
        encounter_block
    min_time, max_time = hours as an int, +/-999 representing infinity.
    agg = aggregation function to pivot table
    '''
    df = in_df[['encounter_block','time_diff',val_col]].copy()
    #If time_delta format convert to hours.
    if pd.api.types.is_timedelta64_dtype(df['time_diff']):
        df['time_diff'] = df['time_diff'].dt.total_seconds()/3600
    
    #Time mask, note that -999 and +999 are treated as infinate
    if min_time == -999:
        time_mask = df['time_diff'] < max_time
    elif max_time == 999:
        time_mask = df['time_diff'] > min_time
    else:
        time_mask = df['time_diff'].between(min_time, max_time)
    df = df[time_mask]
    
    #Rename the variable for usefulness
    new_name = f"{val_col}_{min_time}_{max_time}h_{agg_func}"
    
    #Flag for agg_func = flag or all
    or_flag = agg_func == 'flag'
    if or_flag:
        agg_func = 'max'
        df[val_col] = df[val_col].astype(bool)
    and_flag = agg_func == 'all'
    if and_flag:
        agg_func = 'min'
        df[val_col] = df[val_col].astype(bool)

    print(f'Aggregation step for {new_name}')
    agg_df = df.groupby('encounter_block')[val_col].agg(agg_func).reset_index()
    agg_df.rename(columns={val_col:new_name},inplace=True)

    #Spcial situation for "true" and "all" to convert to boolean
    if or_flag or and_flag: agg_df[new_name] = agg_df[new_name] == 1
    
    return agg_df

class hourly_blocks:

    def __init__(self, in_name:str=None, in_df:pd.DataFrame=None):
        '''
        INIT
        This class does not built the hourly data frame but just holds it and adds to it.
        inputs:
        in_name: name of the file where the data frame is initially stored. It will look for it in "output/intermediate", assumes parquet file.
        in_df: a DataFrame with hour blocks already built in.
        NOTES: The data frame must have columns 'encounter_block':int and 'time_from_vent':int.
        '''
        if in_name:
            self.df = load_data('output_folder',in_name, folder='intermediate')
        elif in_df is not None:
            self.df = in_df.copy()
        else:
            raise ValueError("Input to initialize hourly_blocks must include in_path=file name or in_df=DataFrame.")

        self.required_cols = ['encounter_block', 'time_from_vent']
        missing = [col for col in self.required_cols if col not in self.df.columns]
        if missing:
            raise ValueError(f"Missing columns: {missing}")

        self.df['encounter_block'] = self.df['encounter_block'].astype(int)
        self.df['time_from_vent'] = self.df['time_from_vent'].astype(int)

        self.df = self.df.sort_values(by=self.required_cols)

    def calc_time_from_vent(self, t_diff:pd.Series) -> pd.Series:
        '''
        INPUT
        t_diff = Pandas DataSeries with a time_delta time element or a numerical element assumed to be hours
        OUTPUT
        DataSeries with integers rounded down as the value.
        '''
        #If time_delta format convert to hours.
        if pd.api.types.is_timedelta64_dtype(t_diff):
            t_diff = t_diff.dt.total_seconds()/3600

        _res = np.floor(t_diff)
        return _res.astype(int)

    def hourly_fill(self,
                      val_col:str,
                      fill_with):
        '''
        Sorts time_from_vent and fill them in with specified request.
        INPUT
        val_col = val_col string name in the self.df.
        fill_with = either forward fill or back and forward fill or a value.
        NOTE: Assumes data frame remains sorted.
        '''
        
        if fill_with == 'ffill':
            self.df[val_col] = self.df.groupby('encounter_block')[val_col].ffill()
        elif fill_with == 'bffill':
            #Forward fill first then back fill
            self.df[val_col] = self.df.groupby('encounter_block')[val_col].ffill()
            self.df[val_col] = self.df.groupby('encounter_block')[val_col].bfill()
        else:
            # NB: chained `self.df[col].fillna(..., inplace=True)` is a silent no-op
            # under pandas 3.0 Copy-on-Write. Assign the result back instead.
            self.df[val_col] = self.df[val_col].fillna(fill_with)
    
    def addto_blocks(self,
                     input_df:pd.DataFrame,
                     val_col:str,
                     agg_func:str = 'mean',
                     fill_with = np.nan,
                     new_name:str = None,
                     reorder:bool = False):
        '''
        INPUTS:
        df = DataFrame required to have:
            Either a column called 'time_diff' which is a date time difference value or any number (assumed hours) or 'time_from_vent'.
            val_col, the data itself to be aggregated.
            encounter_block
            'time_from_vent' or 'time_diff' = unique time string can be created with function above.
        agg = aggregation function to aggregate
        new_name = for the column in hourly_blocks
        '''
        #Time bins
        if 'time_from_vent' not in input_df.columns:
            input_df['time_from_vent'] = self.calc_time_from_vent(input_df['time_diff'])
    
        #Copy data frame
        in_df = input_df[['encounter_block','time_from_vent',val_col]].copy()
        
        #Name of new variable
        if not new_name:
            new_name = f"{val_col}_{agg_func}"
    
        #Flag for agg_func = flag or all
        or_flag = agg_func == 'flag'
        if or_flag:
            agg_func = 'max'
            in_df[val_col] = in_df[val_col].astype(bool)
        and_flag = agg_func == 'all'
        if and_flag:
            agg_func = 'min'
            in_df[val_col] = in_df[val_col].astype(bool)
    
        #Aggregation step
        print(f'Aggregation step for {new_name}')
        agg_df = in_df.groupby(by=['encounter_block','time_from_vent'])[val_col].agg(agg_func).reset_index()

        #Rename column
        agg_df.rename(columns={val_col:new_name},inplace=True)

        if reorder:
            agg_df = agg_df.sort_values(self.required_cols)
        
        #Spcial situation for "or" and "and" to convert to boolean and back to int to 0 or 1.
        if or_flag or and_flag:
            agg_df[new_name] = agg_df[new_name] == 1
            agg_df[new_name] = agg_df[new_name].astype(int)
        
        #Merge onto hourly blocks
        print(f'Merging step for {new_name}')
        self.df = pd.merge(self.df,
                           agg_df,
                           on=['encounter_block','time_from_vent'],
                           how='left')

        #Sort and fill
        if fill_with:
            print(f'Sort/Fill step for {new_name}')
            self.hourly_fill(new_name, fill_with)

    def table_summary(self) -> pd.DataFrame:
        return table_summary(self.df)

    def save(self, suffix:str=''):
        path = os.path.join(output_folder,'intermediate',f'hourly_df{suffix}.parquet')
        self.df.to_parquet(path)
        path = os.path.join(output_folder, 'final',f'hourly_df{suffix}_summary.csv')
        _summary = self.table_summary()
        _summary.to_csv(path)