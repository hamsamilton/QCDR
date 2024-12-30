from TestFactory import TestFactory
from QCDR_main import *

class QCDRTestFactory(TestFactory):

    def __init__(self,SaveDir):

        super().__init__(SaveDir)# Initializes self.SaveDir and self.RunInfoDir
        self.DataLoc = "../data/"
        self.fail_alpha = .05
        self.warn_alpha = .1

        self.SCRIPTDataLoc = self.DataLoc + 'SCRIPT/'
        self.SCRIPTCutoffTable = self.SCRIPTDataLoc + "manual_cutoff_table.xslx"
        #Full SCRIPT data paths
        self.SCRIPTAllData = self.SCRIPTDataLoc + "SCRIPT_stats_allbatches.csv"
        self.SCRIPTCountTable = self.SCRIPTDataLoc + "SCRIPT_CountTable.xlsx"
        self.SCRIPT_GC_info = self.SCRIPTDataLoc + "SCRIPT_GC_info.csv"
        # B11 SCRIPT data paths for faster testing and iterating
        self.SCRIPTB11 = self.SCRIPTDataLoc + "SCRIPT_B11stats.csv"
        self.SCRIPTB11CountTable = self.SCRIPTDataLoc + "SCRIPT_CountTable_B11.xlsx"
        self.SCRIPT_B11_GC_info = self.SCRIPTDataLoc + "SCRIPT_B11_GC_info.csv"
        # Lung Transplant files
        self.LungDataLoc = self.DataLoc + 'LungTransplant/'
        self.LungTransplant = self.LungDataLoc + "LungTransplantStats.csv"
        self.LungTransplantGBC = self.LungDataLoc + "LungTransplant_GBC.csv"
        self.LungTransplantCountTable = self.LungDataLoc + "LungTransplantRawCounts.xlsx"
        self.CutoffTable= self.LungDataLoc + "LTcutoffs.xlsx"

        self.BaseTest = {"qry_filename"    : self.SCRIPTB11,
                         "op_folder"       : self.SaveDir + "SCRIPTB11Test",
                         "_bgd_file"       : self.SCRIPTAllData,
                         "_gc_file"        : self.SCRIPT_B11_GC_info,
                         "_hist_file"      : self.SCRIPTB11CountTable,
                         "cutoff_filename" : None,
                         "fail_alpha"      : self.fail_alpha,
                         "warn_alpha"      : self.warn_alpha}

    def MakeAllSCRIPTTest(self):
        # DO a test over all SCRIPT data
        TestName = "SCRIPTWholeDataTest"
        AllSCRIPTTest = self.BaseTest.copy()
        AllSCRIPTTest["qry_filename"] = self.SCRIPTAllData
        AllSCRIPTTest["op_folder"] = self.SaveDir + TestName
        AllSCRIPTTest['_gc_file'] = self.SCRIPT_GC_info
        AllSCRIPTTest['_hist_file'] = self.SCRIPTCountTable

        self.RunInfoDict[TestName] = AllSCRIPTTest

        return self

    def MakeB11SCRIPTTest(self):
        # Do a test over just B11 which can be useful because it runs much faster.
        TestName = "SCRIPTB11Test"

        SCRIPTB11Test = self.BaseTest.copy()
        SCRIPTB11Test["op_folder"] = self.SaveDir + TestName

        self.RunInfoDict[TestName] = SCRIPTB11Test

        return self

    def NoGBCTest(self):
        # Perform a test when the GBC is not added
        TestName =  'NoGBCTest'

        NoGBCTest = self.BaseTest.copy()
        NoGBCTest['op_folder'] = self.SaveDir + TestName
        NoGBCTest['_gc_file'] = None

        self.RunInfoDict[TestName] = NoGBCTest

        return self

    def NoHistTest(self):
        # Perform a test when the Hist is not added
        TestName = 'NoHistTest'

        NoHistTest = self.BaseTest.copy()
        NoHistTest['op_folder'] = self.SaveDir + TestName
        NoHistTest['_hist_file'] = None

        self.RunInfoDict[TestName] = NoHistTest

        return self


    def NoHistGBCTest(self):
        # Perform a test when the Hist and GBC are not added
        TestName = 'NoHistGBCTest'

        NoHistGBCTest = self.BaseTest.copy()
        NoHistGBCTest['op_folder'] = self.SaveDir + TestName
        NoHistGBCTest['_hist_file'] = None

        self.RunInfoDict[TestName] = NoHistGBCTest

        return self

    def ChangeFailandWarn(self):
        # Test how the test performs when the fail and warn are changed
        TestName = 'ChangeWarnAndFail'
        ChangedWarnAndFail = self.BaseTest.copy()
        ChangedWarnAndFail['op_folder'] = self.SaveDir + TestName
        ChangedWarnAndFail['fail_alpha'] = .01
        ChangedWarnAndFail['fail_alpha'] = .3

        self.RunInfoDict[TestName] = ChangedWarnAndFail

        return self

    def LungTransplantTest(self):
        # Perform a test on the lung transplant dataset

        TestName = "LungTransplantTest"
        LungTransplantTest = self.BaseTest.copy()
        LungTransplantTest["qry_filename"] = self.LungTransplant
        LungTransplantTest["op_folder"] = self.SaveDir + TestName
        LungTransplantTest["_bgd_file"] = self.LungTransplant
        LungTransplantTest["_gc_file"] = self.LungTransplantGBC
        LungTransplantTest["_hist_file"] = self.LungTransplantCountTable

        self.RunInfoDict[TestName] = LungTransplantTest

        return self


    def LungTransplantTestwCutoff(self):
        # Perform a test on the lung transplant dataset

        TestName = "LungTransplantwCutoffTest"
        LungTransplantTest = self.BaseTest.copy()
        LungTransplantTest["qry_filename"] = self.LungTransplant
        LungTransplantTest["op_folder"] = self.SaveDir + TestName
        LungTransplantTest["_bgd_file"] = self.LungTransplant
        LungTransplantTest["_gc_file"] = self.LungTransplantGBC
        LungTransplantTest["_hist_file"] = self.LungTransplantCountTable
        LungTransplantTest["cutoff_filename"] = self.CutoffTable


        self.RunInfoDict[TestName] = LungTransplantTest

        return self

    def LungTransplantTestSCRIPTbgd(self):
        # Perform a test on the LungTransplant dataset but with the SCRIPT background

        TestName = "LungTestSCRIPTbgd"
        self.LungTransplantTest()
        LungTestSCRIPTBgd = self.BaseTest.copy()
        LungTestSCRIPTBgd['qry_filename'] = self.LungTransplant
        LungTestSCRIPTBgd['op_folder'] = self.SaveDir + TestName
        LungTestSCRIPTBgd['_bgd_file'] = self.SCRIPTAllData
        LungTestSCRIPTBgd['_gc_file'] = self.LungTransplantGBC
        LungTestSCRIPTBgd['_hist_file'] = self.LungTransplantCountTable

        self.RunInfoDict[TestName] = LungTestSCRIPTBgd

    def ComprehensiveTesting(self):

        # Run a comprehensive list of the above tests for thorough testing
        self.LungTransplantTestSCRIPTbgd()
        self.LungTransplantTest()
        self.MakeAllSCRIPTTest()
        self.MakeB11SCRIPTTest()
        self.NoGBCTest()
        self.NoHistTest()
        self.NoHistGBCTest()
        self.ChangeFailandWarn()
        self.LungTransplantTestwCutoff()

        return self

    def QuickTest(self):
        # Run a limited number of fast tests to iterate quickly and bugfix

        #self.MakeB11SCRIPTTest()
        #self.LungTransplantTest()
        #self.LungTransplantTestSCRIPTbgd()
        self.LungTransplantTestwCutoff()

        return self





if __name__ == '__main__':

    print('Running tests')
    [QCDRTestFactory(SaveDir = "QCDRTestOutputs/").
     ComprehensiveTesting().
     RunTests(TestFun = QCDR_main)]
    print("Tests Finished Running")











