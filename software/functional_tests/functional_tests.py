import traceback
import cProfile
import pstats
import os
import subprocess
import concurrent.futures

class PyTestSuite:
    "Run tests to test QCDR"

    def __init__(self):
        self.path_to_data = "../../data/"
        self.OutputDir = "testoutputs/"
        if not os.path.exists(self.OutputDir):
            os.makedirs(self.OutputDir)
        self.SCRIPT_B11 = "SCRIPT/SCRIPT_B11stats.csv"
        self.SCRIPT_AllBatches = "SCRIPT/SCRIPT_stats_allbatches.csv"
        self.SCRIPT_B11_GCInfo = "SCRIPT/SCRIPT_B11_GC_info.csv"
        self.SCRIPT_man_cutoff_table = "SCRIPT/manual_cutoff_table.xlsx"
        self.SCRIPT_CountTable = "SCRIPT/SCRIPT_CountTable.xlsx"
        self.SCRIPT_B11HistData= "SCRIPT/SCRIPT_B11histdata.5.csv"
        self.SCRIPT_AllGCInfo = "SCRIPT/SCRIPT_GC_info.csv"
        self.LungTransplantStats= "LungTransplant/LungTransplantStats.csv"
        self.LungTransplantGBC  = "LungTransplant/LungTransplant_GBC.csv"
        self.LungTransplantGeneHist = "LungTransplant/LungTransplant_GeneHist.5.csv"
        self.prepend_path_to_data()
        self.commands = self.ConstructCommands()
        print('about to run tests')
        self.RunQCDRTests()
        print('testsrun') 
    def prepend_path_to_data(self):
            attributes = [attr for attr in dir(self) if not callable(getattr(self, attr)) and not attr.startswith("__")]
            for attr in attributes:
                if isinstance(getattr(self, attr), str):
                    setattr(self, attr, self.path_to_data + getattr(self, attr))

    def ConstructCommands(self):
        "Construct the commands for QCDR and run them"
        SCRIPTB11_BaseCase =  {"ip"  :self.SCRIPT_B11,
                     "out" : "testoutputs/SCRIPTB11_BAseCase.pdf",
                     "bgd" : self.SCRIPT_AllBatches,
                     "gc"  : self.SCRIPT_B11_GCInfo,
                     "hist": self.SCRIPT_B11HistData}

        SCRIPTB11_noGCHIST = {"ip"  :self.SCRIPT_B11,
                     "out" :"testoutputs/SCRIPTB11_BaseCasenoGCHIST.pdf",
                     "bgd" : self.SCRIPT_AllBatches}

        SCRIPTAll          = {"ip"  :self.SCRIPT_AllBatches,
                     "out" :"testoutputs/SCRIPTAll.pdf",
                     "bgd" : self.SCRIPT_AllBatches,
                     "gc"  : self.SCRIPT_AllGCInfo}

        LungTransplant = {"ip" : self.LungTransplantStats,
                          "out": "testoutputs/LungTransplantSelfTest.pdf",
                          "bgd": self.LungTransplantStats,
                          "gc" : self.LungTransplantGBC,
                          "hist":self.LungTransplantGeneHist}

        LungTransplantSCRIPTbgd = {"ip" : self.LungTransplantStats,
                          "out": "testoutputs/LungTransplantSCRIPTbgd.pdf",
                          "bgd": self.SCRIPT_AllBatches,
                          "gc" : self.LungTransplantGBC,
                          "hist":self.LungTransplantGeneHist}

        LungTransplantNoGCHist = {"ip" : self.LungTransplantStats,
                          "out": self.OutputDir + "LungTransplantSCRIPTbgd.pdf",
                          "bgd": self.LungTransplantStats}

        commands = [SCRIPTB11_BaseCase]

        return commands

    def ExecuteCommand(self,test):

        TerminalCommand = ["python3","../QCDR_main.py",
                       "-ip", test["ip"],
                       "-out",test["out"],
                       "-bgd",test["bgd"]
                      ]

        if "gc" in test:
            TerminalCommand += ["-gc", test["gc"]]
        if "hist" in test:
            TerminalCommand += ["-hist", test["hist"]]
        if "ctf" in test:
            TerminalCommand += ["-ctf", test["ctf"]]

        with open((test["out"] + ".txt"), 'w') as log_file:
            subprocess.run(TerminalCommand, stdout=log_file, stderr=subprocess.STDOUT)

    def RunQCDRTests(self):

        print(self.commands)
        with concurrent.futures.ProcessPoolExecutor() as executor:
            futures = [executor.submit(self.ExecuteCommand,test) for test in self.commands]

            for future in concurrent.futures.as_completed(futures):
                try:
                    future.result()
                    print('Test is finished')
                except Exception as e:
                    TraceBackStr = ''.join(traceback.format_exception(None,
                                                                      e,
                                                                      e.__traceback__))
                    print(f'A task failed with an exception: {TraceBackStr}')

if __name__ == '__main__':
    TestQCDR = PyTestSuite()
