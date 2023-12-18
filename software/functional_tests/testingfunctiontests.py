import subprocess

class PyTestSuite:
    "Run tests to test QCDR"

    def __init__(self):
        self.path_to_data = "../../data/"
        self.SCRIPT_B11 = "SCRIPT/SCRIPT_B11stats.csv"
        self.SCRIPT_AllBatches = "SCRIPT/SCRIPT_stats_allbatches.csv"
        self.SCRIPT_B11_GCInfo = "SCRIPT/SCRIPT_B11_GC_info.csv"
        self.SCRIPT_man_cutoff_table = "SCRIPT/manual_cutoff_table.xlsx"
        self.SCRIPT_CountTable = "SCRIPT/SCRIPT_CountTable.xlsx"
        self.SCRIPT_B11HistData= "SCRIPT/SCRIPT_B11histdata.5.csv"
        self.LungTransplantStats= "LungTransplant/LungTransplantStats"
        self.LungTransplantGBC  = "LungTransplant/LungTransplant_GBC.csv"
        self.LungTransplantGeneHist = "LungTransplant/LungTransplant_GeneHist.5.csv"
        self.prepend_path_to_data()
        self.commands = self.MakeQCDR_TestSuite()
        self.TerminalCommand = self.construct_sh()
        self.write_commands()
        self.ExecuteCommands()
    def prepend_path_to_data(self):
            attributes = [attr for attr in dir(self) if not callable(getattr(self, attr)) and not attr.startswith("__")]
            for attr in attributes:
                if isinstance(getattr(self, attr), str):
                    setattr(self, attr, self.path_to_data + getattr(self, attr))

    def construct_command(self, ip = "",out = "",gc = "",hist = "",bgd = "",ctf = ""):
        
        
        command = "python3 ../QCDR_main.py -ip {} -out {} -gc {} -hist {} -bgd {} -ctf {}".format(ip, out, gc, hist, bgd, ctf)
        return command

        for command in self.commands:
            sh_command = sh_command + command + "\n"

        return sh_command

    def MakeQCDR_TestSuite(self):
        "Construct the commands for QCDR and run them"
        BaseCase = self.construct_command(ip  = self.SCRIPT_B11,
                                          out = "BaseCase.pdf",
                                          bgd = self.SCRIPT_AllBatches,
                                          gc  = self.SCRIPT_B11_GCInfo,
                                          hist= self.SCRIPT_B11HistData)

        commands = [BaseCase]
        
        return commands        

    def write_commands(self):
        "write commands to a text file for record keeping and bugfixing"
        
        with open("TestRecords.txt","w") as file:
            file.write(self.TerminalCommand)

    def ExecuteCommands(self):
        "send commands to terminal"
        result = subprocess.run(self.TerminalCommand)
        print(result.stdout)


