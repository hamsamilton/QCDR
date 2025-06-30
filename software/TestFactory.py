import os as os
import sys
import copy
import cProfile
import pstats
from abc import ABC, abstractmethod
from concurrent.futures import ProcessPoolExecutor, as_completed
import traceback

class TestFactory(ABC):
  """
  This class provides a framework for testing other functions in parallel.
  It defines a dictionary containing dictionaries, where each inner dict has the arguments to
  send to the test functions
  """
  @abstractmethod
  def __init__(self,SaveDir = '',MaxWorkers = 10):
    """RunInfoDict: A dictionary to be filled with tests with their information on how to run
    SaveDir    : Specify the super folder where output information should be saved
    MaxWorkers : How many distinct jobs do you want to run at once, too many and you'll crash the job"""
    RunInfoDir   = {}
    self.SaveDir    = SaveDir
    self.MaxWorkers = MaxWorkers

    """Subclasses should have an __init__ statement with a dictionary containing the parameters required to run
    the test compatible with **kwargs"""
    self.RunInfoDict = {}
    pass

  @staticmethod
  def MakeTestReplicates(TestDict,NReps):
    # produce replicates of given tests
    NewTestDict = {}

    for TestName, Test in TestDict.items():

      NewTestDict[TestName] = copy.deepcopy(Test)

      for rep in range(1, NReps + 1):

        RepTestName = f"{TestName}_{rep}"
        RepTest = copy.deepcopy(Test)
        RepTest["SaveDir"] += f"_{rep}"

        NewTestDict[RepTestName] = RepTest

    return NewTestDict

  def RunTests(self,TestFun):
    """ Provides a framework for running multiple tests
      RunInfoDict: A dictionary containing dictionaries which stores the parameters
                   for each test. THESE DICTIONARIES MUST HAVE A SAVEDIR KEY SO THAT
                   PLOTS CAN BE SAVED PROPERLY
      SaveDir    : Where should the TestOutputs be stored?
      TestFun   : What test (function) should be run?
    """

    print('List of tests',self.RunInfoDict)
    # Make the external folder 
    os.makedirs(self.SaveDir,
                exist_ok = True)

    print(f'CPUs available for parallelization :{os.cpu_count()}')
    # for each directory path for all dictionaries  in RunInfoDict, add the parent directory
    with ProcessPoolExecutor(max_workers = self.MaxWorkers) as executor:
      futures = []

      # Make a save location for each tests output 
      for RunName, RunInfo in self.RunInfoDict.items():

        logfilename = os.path.join(self.SaveDir,
                                   f"{RunName}.txt")

        future = executor.submit(self.SendJob,
                                 TestFun,
                                 RunInfo,
                                 logfilename)
        futures.append(future)
      for future in as_completed(futures):
        try:
          future.result()
        except Exception as e:
          TraceBackStr = ''.join(traceback.format_exception(None,
                                                            e,
                                                            e.__traceback__))
          print(f"A task failed with an exception: {TraceBackStr}")

  def SendJob(self,TestFun, RunInfo,logfilename):
    # Redirects output to handle the logging per the test run
    with open(logfilename, 'w') as logfile:
      original_stdout = sys.stdout
      original_stderr = sys.stderr
      try:
        sys.stdout = logfile
        sys.stderr = logfile

        #Setup profiler
        profiler = cProfile.Profile()
        profiler.enable()

        print(f'Running Test {TestFun} with parameters {RunInfo}')
        TestFun(**RunInfo)

      except Exception as e:
        TraceBackStr = ''.join(traceback.format_exception(None,
                                                          e,
                                                          e.__traceback__))
        print(f"A task failed with an exception: {TraceBackStr}")

      finally:
        profiler.disable()
        ps = pstats.Stats(profiler,
                          stream = sys.stdout).sort_stats('cumulative')
        ps.print_stats(40)
        #Restore original stdout and stderr
        sys.stdout = original_stdout
        sys.stderr = original_stderr
        print(f'Finished Test{TestFun} with parameters {RunInfo}')

  def ToListOfLists(InputList):
    # A simple utility
    return [[item] for item in InputList]
