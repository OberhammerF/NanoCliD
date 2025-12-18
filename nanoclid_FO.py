import glob
import os
import subprocess
from argparse import ArgumentParser
import pandas as pd
import re
import yaml
import socket
import traceback
try:
   from utils.utils import getOutputs
except:
   pass
try:
    from utils.curie.utils import *
    curieFunctions = curieUtils()
    curieNetwork = True
except:
    curieNetwork = False

class NanoClid:
    _DEFAULT_P_DIR = {"abacus" : "/mnt/beegfs/EH/pipelines/prod/.p", "calcsub" : "/data/bioinfo-clinique-public/prod/.p", "standalone" : "/data/bioinfo-clinique-public/prod/.p"}

    def __init__(self, inputFolder=None, bedDir=None, bedFile=None, run=None, outDir=None, dryRun=None, genomeVersion=None, until=None, samples="", outTemplate=None, snpEffDir=None, copyToTransverse=None, copyToWorkspace=None, runOnCluster=None, transverseFolder = None, sampleSheet = None, refDir = None, hostName = None, snakemakeBin = None, email = "", profile = None, queue = None, containersFolder = None, doradoModelsFolder = None, basecaller=None):
        self.inputFolder = inputFolder
        self.bedDir = bedDir
        self.bedFile = bedFile
        self.run = run
        self.outDir = outDir
        self.dryRun = dryRun
        self.genomeVersion = genomeVersion
        self.until = until
        self.samplesToRun = samples.split(",")
        self.outTemplate = outTemplate
        self.gitDir = os.path.dirname(os.path.realpath(__file__))
        self.snpEffDir = snpEffDir
        self.doradoModelsFolder = doradoModelsFolder
        self.basecaller = basecaller if basecaller else "dorado"
        # containersFolder is now required - must be passed as argument
        if not containersFolder:
            raise ValueError("containersFolder is required. Please provide the path to the singularity containers folder using -C/--containersFolder argument.")
        self.containersFolder = os.path.abspath(containersFolder)
        if not os.path.isdir(self.containersFolder):
            raise ValueError(f"containersFolder does not exist: {self.containersFolder}")
        if snakemakeBin is None and profile is None and curieNetwork:
            self.snakemakeBin, profile = curieFunctions._setSnakemakeBinAndProfile()
        else:
            # If snakemakeBin not provided, try to find it in PATH or venv
            if snakemakeBin is None:
                import shutil
                snakemakeBin = shutil.which("snakemake")
                if snakemakeBin is None:
                    # Try venv in gitDir
                    venv_snakemake = os.path.join(self.gitDir, "venv", "bin", "snakemake")
                    if os.path.exists(venv_snakemake):
                        snakemakeBin = venv_snakemake
                    else:
                        raise ValueError("snakemake not found in PATH. Please provide -B/--snakemakeBin argument or activate the venv.")
            self.snakemakeBin = snakemakeBin
        self.configTemplate = self.__setConfigTemplate(self.gitDir, self.run, curieNetwork, profile)
        self.profile = os.path.join(self.gitDir, "profiles", "curie", profile) if curieNetwork else os.path.join(self.gitDir, "profiles", "externe", profile)
        if refDir == "" and curieNetwork:
            self.refDir = curieFunctions._setRefDir(refDir, genomeVersion, profile)
        else:
            self.refDir = refDir
        self.copyToTransverse = copyToTransverse
        self.copyToWorkspace = copyToWorkspace
        self.runOnCluster = runOnCluster
        self.transverseFolder = transverseFolder
        self.sampleSheet = sampleSheet
        self.analysis = None
        self.combinaison = "''"
        self.fastqConcatenated = False
        self.fromBlow5 = "'no'"
        self.fromFast5 = "'no'"
        self.fromPod5 = "'yes'"
        self.hostName = hostName if hostName else ""
        self.email = email
        self.queue = queue
        if curieNetwork and not self.email:
            self.email = "bioinfo-clinique@curie.fr"

    def __setConfigTemplate(self, gitDir, run, curieNetwork, profile):
        folder = "curie" if curieNetwork else "externe"
        if run == "A000":
            return os.path.join(gitDir, "data", folder, f"test_{profile}.yaml")
        return os.path.join(gitDir, "config", folder, profile, "config.yaml")

    def __checkBedIntegrity(self, bed):
        df = pd.read_csv(bed, sep = "\t", header = None)
        if not "chr" in df[0].iloc[0]:
            df = df.sort_values(by = [0,1])
            df[0] = df[0].astype(str)
            df[0] = "chr" + df[0]
            df.to_csv(bed, sep = "\t", index = None, header = False)

    def parseSampleSheet(self, sampleSheet, bedFile):
        subprocess.call(f'mkdir -p {os.path.join(self.inputFolder, self.run, "archive")}', shell = True)
        if bedFile:
            subprocess.call(f"cp {bedFile} {self.inputFolder}/{self.run}/archive/", shell = True)
            self.bed = f"{self.inputFolder}/{self.run}/archive/{os.path.basename(bedFile)}"
        elif len(glob.glob(f"{self.inputFolder}/{self.run}/archive/*bed")) != 0:
            self.bed = glob.glob(f"{self.inputFolder}/{self.run}/archive/*bed")[0]
        else:
            code = subprocess.call(f'grep TargetBED {sampleSheet}', shell = True)
            if code == 1:
                subprocess.call(f"touch {self.inputFolder}/{self.run}/archive/empty.bed", shell = True)
                self.bed = f"{self.inputFolder}/{self.run}/archive/empty.bed"
            else:
                if "A000" in os.path.basename(sampleSheet):
                    bed = os.path.join(self.gitDir, "data", "test.bed")
                else:
                    bed = subprocess.check_output(f'grep TargetBED {sampleSheet}', shell=True).decode('utf-8').rstrip().split(',')[1]
                    if curieNetwork:
                        bed = curieFunctions._getBed(bed, os.getenv("USER"), self.genomeVersion, self.run)
                    else:
                        bed = os.path.join(self.bedDir, bed)
                subprocess.call(f"cp {bed} {self.inputFolder}/{self.run}/archive/", shell = True)
                self.bed = f"{self.inputFolder}/{self.run}/archive/{os.path.basename(bed)}"
                if not ".bed" in self.bed:
                    self.bed = f"{self.bed}.bed"
        if self.bed != f"{self.inputFolder}/{self.run}/archive/empty.bed":
            self.__checkBedIntegrity(self.bed)
        self.sequencer = ""
        self.flowCellType = subprocess.check_output(f"grep FlowcellType {sampleSheet}", shell = True).decode("utf-8").rstrip().split(",")[1].lower()
        self.barcodingKits = subprocess.check_output(f"grep Assay {sampleSheet}", shell = True).decode("utf-8").rstrip().split(",")[1]
        self.barcodingKits = self.barcodingKits.replace('.', '-')
        skipRows = subprocess.check_output(f'grep -n "\[Data\]" {sampleSheet}', shell = True).decode("utf-8").split(":")[0]
        sampleSheet = pd.read_csv(sampleSheet, sep  = ",", skiprows = int(skipRows), keep_default_na = False)
        sampleSheet = sampleSheet[sampleSheet["Index_ID"] != ""]
        self.mergedSamples = sampleSheet["Index_ID"].tolist()
        sampleSheet["Index_ID"] = sampleSheet["Index_ID"].str[2:]
        sampleSheet["Index_ID"] = "barcode" + sampleSheet["Index_ID"]
        self.samplesPath = dict(zip(self.mergedSamples, sampleSheet.Index_ID))
        self.__setTypeAnalysis(self.inputFolder, self.run, self.mergedSamples)
        self.__setDataPath()

    def __setWildcardsCombinaison(self, injections, fast5Dir):
        self.combinaison = []
        injectionsForMerge = []
        samples = sorted(fast5Dir.keys())
        for injection in injections:
            injectionsForMerge.append(injection)
            for sample in samples:
                if injection in fast5Dir[sample] and "_" in sample:
                    self.combinaison.append({"samples" : sample, "injections" : injection})

    '''
    Sometimes, demultiplexing is not performed, so we need to perform demultiplexing
    '''
    def __isDemultiplexingToDo(self, fast5Paths):
        self.demultiplexing = False
        fast5 = []
        for fast5Pass in fast5Paths.values():
            if isinstance(fast5Pass, list):
                fast5.append(",".join(fast5Pass))
            else:
                fast5.append(fast5Pass)
        fast5 = set(",".join(fast5).split(",")) #get all fast5 paths
        if self.multiplexing and len(fast5) <= len(self.injections):
            self.demultiplexing = True

    '''
    Sometimes, basecalling can be performed on minknown. So if we have the fastq folder, just concat them in order to avoid to redo basecalling
    '''
    def __setFastqPath(self, folders):
        r = re.compile("fastq_pass")
        if len(list(filter(r.search, folders))) != 0:
            return list(filter(r.search, folders))[0]
        return ""

    def __concatFastq(self, fastqPaths, combinaison):
        for sample in fastqPaths.keys():
            if fastqPaths[sample] != "":
                resFolder = f"{self.outDir}"
                if self.analysis == "multipleInjections":
                    resFolder = os.path.join(f"{self.outDir}/{self.run}")
                injection = [dico["injections"] for dico in combinaison if dico["samples"] == sample][0]
                subprocess.call(f"mkdir -p {resFolder}", shell = True)
                subprocess.call(f"mkdir -p {resFolder}/{injection}", shell = True)
                subprocess.call(f"mkdir -p {resFolder}/{injection}/{sample}", shell = True)
                subprocess.call(f"mkdir -p {resFolder}/{injection}/{sample}/FASTQ", shell = True)
                subprocess.call(f"cat {fastqPaths[sample]}/*fastq.gz > {resFolder}/{injection}/{sample}/FASTQ/{sample}.fastq.gz", shell = True)
                self.fastqConcatenated = True

    def __setFast5Path(self, extension):
        self.fast5Paths = {}
        self.fastqPaths = {}
        self.samplesToMerge = {}
        self.samples = []
        samples = sorted(list(self.samplesPath.keys()))
        samples = set(self.samplesToRun).intersection(samples) if self.samplesToRun != [""] else samples
        runDir = f"{self.inputFolder}" if self.analysis == "simpleInjection" else os.path.join(self.inputFolder ,self.run)
        if self.analysis == "simpleInjection":
            runDir = f"{self.inputFolder}"
        for sample in samples:
            if self.multiplexing or subprocess.call(f'find {runDir}/ -type d -regextype posix-egrep -regex ".*(pod5|fast5).*" | grep barcode', shell = True) == 0:
                pattern = self.samplesPath[sample]
            else:
                pattern = f"{extension}*"
            self.samplesToMerge[sample] = []
            self.fast5Paths[sample] = []
            for i in range(len(self.injections)):
                fast5Dir = self.getFiles(os.path.join(runDir, (self.injections[i])), 'd', pattern, True)
                fastqDir = [f for f in fast5Dir if "fastq" in f and not "fail" in f and not "skip" in f]
                fast5Dir = [f for f in fast5Dir if extension in f and not "fail" in f and not "skip" in f]
                if len(fast5Dir) == 0 and pattern == self.samplesPath[sample] and self.multiplexing: #if multiplexing and no fast5 with barcodeNb pattern found look for fast5 folder supposing demultiplexing is needed
                    fast5Dir = self.getFiles(os.path.join(runDir, (self.injections[i])), 'd', f'{extension}*', True)
                    fastqDir = [f for f in fast5Dir if "fastq" in f and not "fail" in f and not "skip" in f]
                    fast5Dir = [f for f in fast5Dir if extension in f and not "fail" in f and not "skip" in f]
                if len(fast5Dir) > 0:
                    self.fastqPaths[f"{sample}_{i+1}"] = self.__setFastqPath(fastqDir)
                    self.fast5Paths[sample].append(fast5Dir[0])
                    self.fast5Paths[f"{sample}_{i+1}"] = fast5Dir[0]
                    self.samples.append(f"{sample}_{i+1}")
                    self.samplesToMerge[sample].append(f"{sample}_{i+1}")
        samplesToRemove = []
        for sample in self.fast5Paths.keys():
            if len(self.fast5Paths[sample]) == 0:
                samplesToRemove.append(sample)
        for sample in samplesToRemove:
            del self.fast5Paths[sample]
            if not "_" in sample:
                del self.samplesToMerge[sample]
        self.mergedSamples = [sample for sample in self.fast5Paths.keys() if "_" not in sample]

    def __setReportFilesPath(self, fast5Paths, injections):
        self.reportFiles = {}
        for sample in fast5Paths.keys():
            if self.analysis == "simpleInjection":
                self.reportFiles[sample] = self.getFiles(os.path.join(self.inputFolder, self.run), "f", "report*.md", False)
            else:
                for injection in injections:
                    if injection in fast5Paths[sample]:
                        reportFile = self.getFiles(os.path.join(self.inputFolder, self.run, injection), "f", "report*.md", False)
                        if reportFile == "":
                            self.reportFiles[sample] = '""'
                        else:
                            self.reportFiles[sample] = reportFile

    def __setDataPath(self):
        self.__setFast5Path("pod5")
        if self.fast5Paths == {} and self.fastqPaths == {}:
            self.__setFast5Path("fast5")
            if self.fast5Paths != {}:
                self.fromFast5 = "'yes'"
            if self.fast5Paths == {}:
                self.__setFast5Path("blow5")
                self.fromBlow5 = "'yes'"
        self.__setReportFilesPath(self.fast5Paths, self.injections)
        self.__setWildcardsCombinaison(sorted(self.injections), self.fast5Paths)
        self.__isDemultiplexingToDo(self.fast5Paths)
        self.__concatFastq(self.fastqPaths, self.combinaison)

    def __setTypeAnalysis(self, inputFolder, run, samples):
        self.multiplexing = False
        if len(samples) > 1:
            self.multiplexing = True
        folders = [folder for folder in os.listdir(f"{inputFolder}/{run}/") if os.path.isdir(f"{inputFolder}/{run}/{folder}")]
        injections = [folder for folder in folders if re.match(f"{run}_[0-9][0-9]*", folder)]
        self.analysis = "simpleInjection"
        if len(injections) > 0:
            self.analysis = "multipleInjections"
            self.injections = sorted([f"{injection}" for injection in injections])
        else:
            self.injections = [run]


    def parseLineConfig(self, line):
        if line.rstrip() != "":
            key, value = line.rstrip().split(":", 1)[0], line.rstrip().split(":", 1)[1]
            value = value.replace(" ", "", 1)
            if key[-1] == " ":
                key = key.replace(" ", "", 1)
            return key, value
        return "", ""

    def loadConfig(self, path):
        with open(path, "r") as file:
            configYaml = yaml.safe_load(file)
        return configYaml

    def writeConfig(self, dico, path):
        if os.path.basename(path) == "config.yaml":
            # This is a profile config - use yaml.dump for proper formatting
            with open(path, "w") as file:
                yaml.dump(dico, file, sort_keys=False)
        else:
            # This is a run config - use custom format
            f = open(path, "w")
            for key in dico.keys():
                if type(dico[key]) != dict:
                    f.write(f"{key}: {dico[key]}\n")
                else:
                    f.write(f'{key}: \n')
                    for subKey in dico[key].keys():
                        if dico[key][subKey] == "":
                            f.write(f" {subKey}: ''\n")
                        else:
                            f.write(f" {subKey}: {dico[key][subKey]}\n")
                f.write("\n")
            f.close()

    def __updateConfig(self, config):
        
        # Ensure all tool sections are dicts (template may have them as strings or missing)
        tool_sections = [
            "dorado", "guppy", "minimap2", "slow5tools", "slow5_merge", "bluecrab",
            "f5c", "f5c_call_methylation", "samtools_sort", "samtools_merge",
            "bedtools", "pod5tools", "nanoplot", "mosdepth", "samtools_stats",
            "clair3", "pepper", "nanocaller", "sniffles", "cuteSV", "svim", "nanovar",
            "annotSV", "combineVariants", "bioInfoCliTools", "R", "cnv", "circos",
            "computeQC", "concatSV", "snpEff", "wildcards"
        ]
        for section in tool_sections:
            if not isinstance(config.get(section), dict):
                config[section] = {}

        # Set minimal defaults for all tools
        config["dorado"].setdefault("parameters", "")
        config["dorado"].setdefault("model", config.get("dorado", {}).get("model", ""))
        config["dorado"].setdefault("sif", "dorado.sif")
        config["dorado"].setdefault("demux", "")
        # Set dorado models path from command line argument
        config["dorado"]["models_path"] = self.doradoModelsFolder if self.doradoModelsFolder else ""
        
        config["guppy"].setdefault("parameters", config.get("guppy", {}).get("parameters", ""))
        config["guppy"].setdefault("parameters_standalone", config.get("guppy", {}).get("parameters_standalone", ""))
        config["guppy"].setdefault("sif", "guppy.sif")
        
        config["minimap2"].setdefault("mmi", "")
        config["minimap2"].setdefault("parameters", "")
        config["minimap2"].setdefault("sif", "minimap2.sif")
        config["minimap2"].setdefault("cn", "")
        
        config["wildcards"].setdefault("samples", getattr(self, "samples", []))
        config["wildcards"].setdefault("injections", getattr(self, "injections", []))
        config["wildcards"].setdefault("mergedSamples", getattr(self, "mergedSamples", []))
        config["wildcards"].setdefault("run", [getattr(self, "run", "")])

        config["slow5tools"].setdefault("parameters", "")
        config["slow5tools"].setdefault("sif", "slow5tools.sif")
        
        config["slow5_merge"].setdefault("parameters", "")
        
        config["bluecrab"].setdefault("parameters", "")
        config["bluecrab"].setdefault("sif", "bluecrab.sif")
        
        config["f5c"].setdefault("bin", "f5c")
        config["f5c"].setdefault("parameters", "")
        config["f5c"].setdefault("sif", "f5c.sif")
        
        config["f5c_call_methylation"].setdefault("parameters", "")
        
        config["samtools_sort"].setdefault("parameters", "")
        config["samtools_sort"].setdefault("sif", "samtools.sif")
        
        config["samtools_merge"].setdefault("parameters", "")
        config["samtools_merge"].setdefault("sif", "samtools.sif")
        
        config["bedtools"].setdefault("sif", "bedtools.sif")
        
        config["pod5tools"].setdefault("sif", "pod5tools.sif")

        config["nanoplot"].setdefault("parameters", "")
        config["nanoplot"].setdefault("script", "")
        config["nanoplot"].setdefault("sif", "nanoplot.sif")

        config["mosdepth"].setdefault("parameters", "")
        config["mosdepth"].setdefault("sif", "mosdepth.sif")

        config["samtools_stats"].setdefault("parameters", "")

        if "calculate_methylation_frequency" not in config:
            config["calculate_methylation_frequency"] = "-c 2.5 -s"

        config["clair3"].setdefault("parameters", "")
        config["clair3"].setdefault("sif", "clair3.sif")
        config["clair3"].setdefault("threads", 1)

        config["pepper"].setdefault("sif", "pepper.sif")
        config["pepper"].setdefault("threads", 1)

        config["nanocaller"].setdefault("parameters", "")
        config["nanocaller"].setdefault("script", "")
        config["nanocaller"].setdefault("sif", "nanocaller.sif")
        config["nanocaller"].setdefault("threads", 1)

        config["sniffles"].setdefault("parameters", "")
        config["sniffles"].setdefault("sif", "sniffles.sif")

        config["cuteSV"].setdefault("parameters", "")
        config["cuteSV"].setdefault("sif", "cuteSV.sif")
        config["cuteSV"].setdefault("threads", 1)

        config["svim"].setdefault("parameters", "")
        config["svim"].setdefault("qual", "")
        config["svim"].setdefault("sif", "svim.sif")

        config["nanovar"].setdefault("parameters", "")
        config["nanovar"].setdefault("sif", "nanovar.sif")
        config["nanovar"].setdefault("threads", 1)

        config["annotSV"].setdefault("parameters", "")
        config["annotSV"].setdefault("sif", "annotSV.sif")

        config["combineVariants"].setdefault("javaParameters", "")
        config["combineVariants"].setdefault("parameters", "")
        config["combineVariants"].setdefault("sif", "combineVariants.sif")

        config["bioInfoCliTools"].setdefault("sif", "bioInfoCliTools.sif")

        config["R"].setdefault("sif", "R.sif")

        config["cnv"].setdefault("bigwigsif", "")
        config["cnv"].setdefault("brain_bed", "")
        config["cnv"].setdefault("brain_bed_5mb", "")
        config["cnv"].setdefault("cnvfrombamsif", "")
        config["cnv"].setdefault("deeptoolssif", "")
        config["cnv"].setdefault("effectiveGenomeSize", "")
        config["cnv"].setdefault("genome_subsampling", "")
        config["cnv"].setdefault("script", "")
        config["cnv"].setdefault("threads", 1)

        config["circos"].setdefault("chromosome", "")
        config["circos"].setdefault("cytobande", "")
        config["circos"].setdefault("geneList", "")
        config["circos"].setdefault("gtf", "")
        config["circos"].setdefault("script", "")

        config["computeQC"].setdefault("script", "")

        config["concatSV"].setdefault("script", "")

        config["snpEff"].setdefault("javaParameters", "-Xmx8G")
        config["snpEff"].setdefault("parameters", "")
        config["snpEff"].setdefault("sif", "snpEff.sif")

        # dataDir must be a dict with genome versions as keys (required by snv_calling.snk line 185)
        # The Snakemake file does: config["snpEff"]["dataDir"]["hg19"]
        if "dataDir" not in config["snpEff"]:
            snpeff_path = getattr(self, "snpEffDir", "")
            config["snpEff"]["dataDir"] = {
                "hg19": snpeff_path,
                "hg38": snpeff_path
            }
        elif isinstance(config["snpEff"]["dataDir"], str):
            # Convert string to dict format
            path = config["snpEff"]["dataDir"]
            config["snpEff"]["dataDir"] = {
                "hg19": path,
                "hg38": path
            }
        elif isinstance(config["snpEff"]["dataDir"], dict):
            # Already a dict - ensure required genome versions exist
            snpeff_path = getattr(self, "snpEffDir", "")
            config["snpEff"]["dataDir"].setdefault("hg19", snpeff_path)
            config["snpEff"]["dataDir"].setdefault("hg38", snpeff_path)

        # Basic flags
        config["analysis"] = getattr(self, "analysis", "")
        config["demultiplexing"] = getattr(self, "demultiplexing", False)

        # Flowcell/model strings (only if flowCellType available)
        flow = getattr(self, "flowCellType", "")
        seq = getattr(self, "sequencer", "")
        if config["dorado"].get("model"):
            if "r10" in flow:
                config["dorado"]["model"] = config["dorado"]["model"].replace("FLOWCELL", f"dna_{flow}_e8.2_400bps_hac@v5.0.0")
            else:
                config["dorado"]["model"] = config["dorado"]["model"].replace("FLOWCELL", f"dna_{flow}_e8_hac@v3.3")
        else:
            config["dorado"]["model"] = f"dna_{flow}_e8_hac@v3.3" if flow else ""

        # Set guppy model config file based on flow cell type
        # Model files are inside container at /opt/ont-guppy/data/
        if flow:
            if "r10" in flow:
                guppy_model_cfg = f"dna_{flow}_e8.2_400bps_hac{seq}.cfg"
            else:
                guppy_model_cfg = f"dna_{flow}_450bps_hac{seq}.cfg"
        else:
            guppy_model_cfg = "dna_r9.4.1_450bps_hac.cfg"  # default fallback
        
        # Ensure guppy section exists and set model
        if "guppy" not in config:
            config["guppy"] = {}
        config["guppy"]["model"] = guppy_model_cfg

        # Set default guppy standalone parameters if not set (for GPU basecalling)
        if not config["guppy"].get("parameters_standalone"):
            config["guppy"]["parameters_standalone"] = "--device cuda:0"

        if config["guppy"].get("parameters_standalone"):
            if "r10" in flow:
                cfg = f"dna_{flow}_e8.2_400bps_hac{seq}.cfg"
            else:
                cfg = f"dna_{flow}_450bps_hac{seq}.cfg"
            config["guppy"]["parameters_standalone"] = config["guppy"]["parameters_standalone"].replace("FLOWCELL", cfg)

        # Replace FLOWCELL placeholder in guppy.parameters if present
        try:
            if "FLOWCELL" in config["guppy"].get("parameters", ""):
                if 'cfg' not in locals():
                    if "r10" in flow:
                        cfg = f"dna_{flow}_e8.2_400bps_hac{seq}.cfg"
                    else:
                        cfg = f"dna_{flow}_450bps_hac{seq}.cfg"
                config["guppy"]["parameters"] = config["guppy"]["parameters"].replace("FLOWCELL", cfg)
        except Exception:
            pass

        # If a guppy .sif exists in the containers folder, set its absolute path in the run config
        try:
            if getattr(self, "containersFolder", None):
                import glob
                guppy_candidates = glob.glob(os.path.join(self.containersFolder, "guppy*.sif"))
                if guppy_candidates:
                    sel = next((c for c in guppy_candidates if os.path.basename(c) == "guppy.sif"), guppy_candidates[0])
                    config["guppy"]["sif"] = sel
        except Exception:
            pass



        # Demultiplexing params
        if getattr(self, "demultiplexing", False):
            bk = getattr(self, "barcodingKits", "")
            config["dorado"]["demux"] = f'--emit-fastq --kit-name "{bk}"'
            config["guppy"]["parameters_standalone"] = f'{config["guppy"]["parameters_standalone"]} --barcode_kits "{bk}"'

        # Basic run/ref settings
        config["flowCellType"] = flow.split(".")[0] if flow else ""
        config["bed"] = getattr(self, "bed", "")
        config["sampleSheet"] = getattr(self, "sampleSheet", "")
        config["run"] = getattr(self, "run", "")
        config["genome_version"] = getattr(self, "genomeVersion", "")

        # Set containers_path - use the required containersFolder argument directly
        config["containers_path"] = self.containersFolder
        if not config["containers_path"].endswith("/"):
            config["containers_path"] = config["containers_path"] + "/"

        # Genome / minimap2 paths
        refdir = getattr(self, "refDir", "")
        gv = getattr(self, "genomeVersion", "")
        
          # DEBUG: Print what we have
        print(f"DEBUG __updateConfig: self.refDir = {self.refDir}")
        print(f"DEBUG __updateConfig: refdir = {refdir}")
        print(f"DEBUG __updateConfig: gv = {gv}")
        
        
        if refdir and gv:
            # Convert to absolute path so Snakemake can find files from any working directory
            refdir_abs = os.path.abspath(refdir)
            
            # Check which FASTA format exists
            fasta_f = os.path.join(refdir_abs, f"{gv}.fasta")
            fasta_a = os.path.join(refdir_abs, f"{gv}.fa")
            if os.path.exists(fasta_f):
                config["genome"] = fasta_f
            elif os.path.exists(fasta_a):
                config["genome"] = fasta_a
            else:
                config["genome"] = ""
            
            # Set genome file and minimap2 index with absolute paths
            config["genomeFile"] = os.path.join(refdir_abs, f"{gv}.genome")
            config["minimap2"]["mmi"] = os.path.join(refdir_abs, f"{gv}.mmi")
        else:
            config["genome"] = ""
            config["genomeFile"] = ""
            config["minimap2"]["mmi"] = ""

        config["samplesToMerge"] = getattr(self, "samplesToMerge", {})
        config["combinaison"] = getattr(self, "combinaison", "")
        config["git_dir"] = self.gitDir
        config["input_dir"] = self.inputFolder
        config["output_dir"] = f"{self.outDir}"
        config["template"] = self.outTemplate
        config["reportFiles"] = getattr(self, "reportFiles", {})
        config["fast5Dir"] = getattr(self, "fast5Paths", {})

        # Apply snpEffDir for non-Curie networks (update dict values, don't replace dict)
        if not curieNetwork and hasattr(self, "snpEffDir") and self.snpEffDir:
            config["snpEff"]["dataDir"]["hg19"] = self.snpEffDir
            config["snpEff"]["dataDir"]["hg38"] = self.snpEffDir

        config["fromBlow5"] = getattr(self, "fromBlow5", "'no'")
        config["fromFast5"] = getattr(self, "fromFast5", "'no'")
        config["fromPod5"] = getattr(self, "fromPod5", "'yes'")

        if self.email:
            config["email"] = self.email

        # --- set preferred basecaller so Snakemake config reflects choice ---
        config["basecaller"] = getattr(self, "basecaller", "dorado")
        if config["basecaller"] == "guppy":
            # avoid binding dorado models and ensure guppy section present
            config["dorado"]["models_path"] = ""
            config["dorado"]["sif"] = ""
            config["guppy"].setdefault("sif", "guppy.sif")


        if curieNetwork:
            config = curieFunctions._addSpecificCurieInfoToConfig(config, self.run, self.hostName, self.gitDir, self.email, self.transverseFolder)

        if self.runOnCluster:
            config["runAllAnalysisOnCluster"] = "yes"

        return config

    def __createConfig(self, template, folder, run):
        
        os.makedirs(folder, exist_ok=True)
        configFile = f"{run}.yaml"
        config = self.loadConfig(template)
        
        
        # Remove any tool sections from template that are not dicts (force __updateConfig to set them)
        tool_sections = [
            "dorado", "guppy", "minimap2", "slow5tools", "slow5_merge", "bluecrab",
            "f5c", "f5c_call_methylation", "samtools_sort", "samtools_merge",
            "bedtools", "pod5tools", "nanoplot", "mosdepth", "samtools_stats",
            "clair3", "pepper", "nanocaller", "sniffles", "cuteSV", "svim", "nanovar",
            "annotSV", "combineVariants", "bioInfoCliTools", "R", "cnv", "circos",
            "computeQC", "concatSV", "snpEff"
        ]
        for section in tool_sections:
            if section in config and not isinstance(config[section], dict):
                del config[section]  # remove bad entries, __updateConfig will recreate them
                config[section] = {}
        config = self.__updateConfig(config)
        
        
        self.writeConfig(config, os.path.join(folder, configFile))
        
        
        return os.path.join(folder, configFile)

    def __updateProfile(self, profile, outDir, run, queue):
        env = "prod" if "prod" in self.gitDir else "dev"
        profileType = profile.split("/")[-1]
        subprocess.call(f"cp -r {profile} {outDir}/{run}", shell=True)
        profileDico = self.loadConfig(f"{outDir}/{run}/{profileType}/config.yaml")
        
        # Update singularity-args with actual paths
        profileDico["singularity-args"] = profileDico["singularity-args"].replace("GIT_DIR", self.gitDir)
        if os.path.exists(self.refDir):
            profileDico["singularity-args"] = profileDico["singularity-args"].replace("REF_DIR", self.refDir)
        else:
            profileDico["singularity-args"] = profileDico["singularity-args"].replace("REF_DIR,", "")
        profileDico["singularity-args"] = profileDico["singularity-args"].replace("PDIR_VAR", NanoClid._DEFAULT_P_DIR[profileType])
        
        # Add --nv flag for GPU support if not already present
        if "--nv" not in profileDico["singularity-args"]:
            profileDico["singularity-args"] = "--nv " + profileDico["singularity-args"]
        # Add dorado models bind mount if provided (only once)
        try:
            dorado_bind = getattr(self, "doradoModelsFolder", None)
            if dorado_bind and os.path.exists(dorado_bind):
                # Only append if not already present in the singularity-args string
                if dorado_bind not in profileDico["singularity-args"]:
                    profileDico["singularity-args"] += f",{dorado_bind}"
        except Exception:
            pass

        # Set singularity-prefix directly from containersFolder argument
        profileDico["singularity-prefix"] = self.containersFolder
        if not profileDico["singularity-prefix"].endswith("/"):
            profileDico["singularity-prefix"] = profileDico["singularity-prefix"] + "/"
        
        if queue:
            idxQueue = [i for i in range(len(profileDico['default-resources'])) if 'partition' in profileDico['default-resources'][i]][0]
            profileDico['default-resources'][idxQueue] = f'slurm_partition={queue}'
        else:
            profileDico["default-resources"] = ",".join(profileDico["default-resources"]).replace("ENV", env).split(",")
        if "cluster" in profileDico.keys():
            profileDico["cluster"] = profileDico["cluster"].replace("logs_cluster", f"{self.outDir}/{self.run}/logs_cluster")
        
        self.writeConfig(profileDico, f"{outDir}/{run}/{profileType}/config.yaml")
        if "standalone" in profile and not self.until:
            #standalone and no until means you want to run the pipe in standalone configuration, so we remove run_on_cluster.txt from expected output
            newTemplate = os.path.join(self.outDir, self.run, os.path.basename(self.outTemplate))
            subprocess.call(f"mkdir -p {os.path.join(self.outDir, self.run)}", shell = True)
            subprocess.call(f"cp {self.outTemplate} {newTemplate}", shell = True)
            subprocess.call(f"sed -i '/run_on_cluster.txt/d' {newTemplate}", shell = True)
            self.outTemplate = newTemplate
        profile = f"{outDir}/{run}/{profileType}"
        return profile

    def runSnakemake(self, snakemakeBin, snakefile, configFile, profile, dryRun):
        cmd = f"{snakemakeBin} -s {snakefile} --profile {profile} --configfile {configFile} -d {os.path.join(self.outDir, self.run)}"
        demultiplexingRule = os.path.join(self.gitDir, "workflow/rules/demultiplexing.snk")
        if self.demultiplexing:
            cmdDemultiplexing = f"{snakemakeBin} -s {demultiplexingRule} --profile {profile} --configfile {configFile} -d {os.path.join(self.outDir, self.run)}"
            cmd = f"{cmdDemultiplexing} -t && {cmd} --rerun-triggers mtime"
        if self.until and "standalone" in profile:
            cmd = " ".join((cmd, f"-U {self.until}"))
        if dryRun:
            cmds = cmd.split(" && ")
            cmds = [f"{cmd} -n" for cmd in cmds]
            cmd = cmds[0]
            if len(cmds) > 1:
                cmd = " && ".join(cmds)
            print(cmd)
            subprocess.call(f"{cmd}", shell=True)
            return
        cmdNanoclidFinished = f"touch {self.outDir}/{self.run}/nanoclid_done.txt"
        cmd = " && ".join((cmd, cmdNanoclidFinished))
        print(cmd)
        subprocess.call(f"{cmd}", shell=True)

    @staticmethod
    def getFiles(path, kind, pattern, files=False):
        find = subprocess.check_output(f"find -L {path} -type {kind} -name '{pattern}' 2> /dev/null", shell=True).decode("utf-8").rstrip().split("\n")
        if len(find) == 0:
            print(f"No files found in path {path} with pattern {pattern}")
            return ""
        elif files:
            if pattern == "fast5*":
                return [file for file in find if not "fail" in file]
            return sorted(find)
        else:
            return find[0]

    def getSpecificPath(self, runDir, pattern, kind, injections):
        summaryFiles = self.getFiles(runDir, kind, pattern, files=True)
        summaryFilesPath = {}
        for injection in injections:
            for summaryFile in summaryFiles:
                if injection in summaryFile:
                    summaryFilesPath[injection] = summaryFile
                    break
        if pattern == "report*.md":
            injectionWithNoReport = set(injection).difference(set(summaryFilesPath.keys()))
            for injection in injectionWithNoReport:
                summaryFilesPath[injection] = ""
        if summaryFilesPath == {}:
            for injection in injection:
                summaryFilesPath[injection] = ""
        if pattern != "report*.md":
            self.injections = list(set(injections).intersection(set(summaryFilesPath.keys())))
        return summaryFilesPath

    def _runNanoClid(self):
        
        if "abacus" in self.profile and curieNetwork:
            try:
                curieFunctions._deleteOldRun(self.inputFolder, self.run)
            except:
                pass #cas ou on lance avec un autre user qui n'a pas les droits pour supprimer un run deja existant
        if not self.sampleSheet and curieNetwork:
            self.sampleSheet = curieFunctions._getSampleSheet(self.run, self.inputFolder, self.loadConfig(self.configTemplate), self.loadConfig(os.path.join(self.profile, "config.yaml")), self.profile)
        self.parseSampleSheet(self.sampleSheet, self.bedFile)
        configFile = self.__createConfig(self.configTemplate, f"{self.outDir}/{self.run}", self.run)
        
        
        profile = self.__updateProfile(self.profile, self.outDir, self.run, self.queue)
        
        
        if curieNetwork:
            curieFunctions.runSnakemake(self.snakemakeBin, os.path.join(self.gitDir, "workflow/Snakefile"), configFile, profile, self.dryRun, self.outDir, self.run, self.gitDir, self.fastqConcatenated, self.until, self.demultiplexing, self.copyToTransverse, self.copyToWorkspace, self.runOnCluster, self.fromBlow5, self.fromFast5)
        else:
            self.runSnakemake(self.snakemakeBin, os.path.join(self.gitDir, "workflow/Snakefile"), configFile, profile, self.dryRun)


    def _install(self, singularityFolder):
        """
        Install NanoCliD: download images, nanovar, set up venv, and retrieve snpEff annotations.
        Each step is validated before proceeding to the next.
        """
        print(f"Installing NanoCliD to: {singularityFolder}")
        gitDir = os.path.dirname(os.path.realpath(__file__))
        
        # URLs for downloads
        imagesPath = "http://xfer.curie.fr/get/moe5ZZ8PGPD/images.tar.gz"
        nanovarPath = "http://xfer.curie.fr/get/Al1MrbPdC0C/nanovar.tar.gz"
        snpEffPath = "http://downloads.sourceforge.net/project/snpeff/databases/v4_3/snpEff_v4_3_hg19.zip"
        
        # Step 1: Create singularity folder
        print("\n[1/7] Creating singularity folder...")
        code = subprocess.call(f"mkdir -p {singularityFolder}", shell=True)
        if code != 0:
            raise ValueError(f"Failed to create singularity folder: {singularityFolder}")
        if not os.path.isdir(singularityFolder):
            raise ValueError(f"Singularity folder does not exist after creation: {singularityFolder}")
        print("[1/7] Singularity folder created OK")
        
        # Step 2: Download images
        print("\n[2/7] Downloading images...")
        imagesArchive = os.path.join(singularityFolder, "images.tar.gz")
        code = subprocess.call(f"wget -P {singularityFolder}/ {imagesPath}", shell=True)
        if code != 0:
            raise ValueError(f"Download images failed (exit code {code}). URL: {imagesPath}")
        if not os.path.isfile(imagesArchive):
            raise ValueError(f"Images archive not found after download: {imagesArchive}")
        print("[2/7] Download images OK")
        
        # Step 3: Untar images
        print("\n[3/7] Extracting images archive...")
        code = subprocess.call(
            f"tar -xvf {imagesArchive} -C {singularityFolder} && "
            f"mv {singularityFolder}/images/* {singularityFolder} && "
            f"rm -rf {singularityFolder}/images {imagesArchive}",
            shell=True
        )
        if code != 0:
            raise ValueError(f"Untar images failed (exit code {code}). Check download integrity.")
        # Verify at least one .sif file exists
        sif_files = [f for f in os.listdir(singularityFolder) if f.endswith('.sif')]
        if not sif_files:
            raise ValueError(f"No .sif files found in {singularityFolder} after extraction. Extraction may have failed.")
        print(f"[3/7] Extracted {len(sif_files)} container images OK")
        
        # Step 4: Download nanovar archive
        print("\n[4/7] Downloading nanovar archive...")
        annotationsDir = os.path.join(gitDir, "annotations")
        nanovarArchive = os.path.join(annotationsDir, "nanovar.tar.gz")
        code = subprocess.call(f"mkdir -p {annotationsDir}", shell=True)
        if code != 0:
            raise ValueError(f"Failed to create annotations folder: {annotationsDir}")
        code = subprocess.call(f"wget -P {annotationsDir}/ {nanovarPath}", shell=True)
        if code != 0:
            raise ValueError(f"Download nanovar archive failed (exit code {code}). URL: {nanovarPath}")
        if not os.path.isfile(nanovarArchive):
            raise ValueError(f"Nanovar archive not found after download: {nanovarArchive}")
        print("[4/7] Download nanovar archive OK")
        
        # Step 5: Untar nanovar
        print("\n[5/7] Extracting nanovar archive...")
        code = subprocess.call(
            f"tar -xvf {nanovarArchive} -C {annotationsDir} && rm {nanovarArchive}",
            shell=True
        )
        if code != 0:
            raise ValueError(f"Untar nanovar folder failed (exit code {code}). Check git-lfs or download integrity.")
        nanovarDir = os.path.join(annotationsDir, "nanovar")
        if not os.path.isdir(nanovarDir):
            raise ValueError(f"Nanovar folder not found after extraction: {nanovarDir}")
        print("[5/7] Nanovar extraction OK")
        
        # Step 6: Set up Python virtual environment
        print("\n[6/7] Setting up Python virtual environment...")
        venvDir = os.path.join(gitDir, "venv")
        
        # Install virtualenv
        code = subprocess.call("pip3 install virtualenv", shell=True)
        if code != 0:
            raise ValueError(f"Failed to install virtualenv (exit code {code})")
        
        # Create venv
        code = subprocess.call(f"mkdir -p {venvDir} && python3 -m venv {venvDir}", shell=True)
        if code != 0:
            raise ValueError(f"Failed to create virtual environment at {venvDir} (exit code {code})")
        if not os.path.isdir(venvDir):
            raise ValueError(f"Virtual environment directory not created: {venvDir}")
        
        # Upgrade pip/wheel/setuptools
        try:
            subprocess.run(
                ['bash', '-c', f"source {venvDir}/bin/activate && pip3 install --upgrade pip wheel setuptools"],
                check=True
            )
        except subprocess.CalledProcessError as e:
            raise ValueError(f"Failed to upgrade pip/wheel/setuptools: {e}")
        
        # Install snakemake and pandas
        try:
            subprocess.run(
                ['bash', '-c', f"source {venvDir}/bin/activate && pip3 install snakemake pandas"],
                check=True
            )
        except subprocess.CalledProcessError as e:
            raise ValueError(f"Failed to install snakemake/pandas: {e}")
        
        # Verify snakemake is installed
        snakemakeBin = os.path.join(venvDir, "bin", "snakemake")
        if not os.path.isfile(snakemakeBin):
            raise ValueError(f"Snakemake not found after installation: {snakemakeBin}")
        print("[6/7] Python virtual environment setup OK")
        
        # Step 7: Retrieve snpEff annotations
        print("\n[7/7] Retrieving snpEff annotations...")
        snpEffZip = os.path.join(annotationsDir, "snpEff_v4_3_hg19.zip")
        code = subprocess.call(f"wget -P {annotationsDir}/ {snpEffPath}", shell=True)
        if code != 0:
            raise ValueError(f"Failed to download snpEff annotations (exit code {code}). URL: {snpEffPath}")
        if not os.path.isfile(snpEffZip):
            raise ValueError(f"snpEff zip not found after download: {snpEffZip}")
        
        code = subprocess.call(f"cd {annotationsDir} && unzip -o {snpEffZip} && rm {snpEffZip}", shell=True)
        if code != 0:
            raise ValueError(f"Failed to unzip snpEff annotations (exit code {code})")
        
        snpEffDataDir = os.path.join(annotationsDir, "data")
        if not os.path.isdir(snpEffDataDir):
            raise ValueError(f"snpEff data directory not found after extraction: {snpEffDataDir}")
        print("[7/7] snpEff annotations OK")
        
        # Installation complete
        print("\n" + "=" * 60)
        print("Installation completed successfully!")
        print("=" * 60)
        print(f"  Singularity images: {singularityFolder}")
        print(f"  Virtual environment: {venvDir}")
        print(f"  Annotations: {annotationsDir}")
        print("\nNanoCliD is ready to use. Please run the test command to verify.")
        print("=" * 60)
        
        
    def _runTest(self, refDir, profile, email, containersFolder, doradoModelsFolder):
        gitDir = os.path.dirname(os.path.realpath(__file__))
        print("Extracting data test folder...")
        code = subprocess.call(f"tar -xvf {os.path.join(gitDir, 'data', 'input', 'input.tar.gz')} -C {os.path.join(gitDir, 'data')}", shell = True)
        if code == 1:
            raise ValueError("Extracting data test folder OK")
        inputFolder = os.path.join(gitDir, "data")
        bedDir = os.path.join(gitDir, "data")
        run = "A000"
        outDir = os.path.join(gitDir, "data")
        dryRun = False
        genomeVersion = "hg19"
        until = False
        samples = ""
        outTemplate = os.path.join(gitDir, "data", "externe", "TEST.template")
        snpEffDir = os.path.join(gitDir, "annotations", "data")
        copyToTransverse = False
        copyToWorkspace = False
        sampleSheet = os.path.join(gitDir, "data", "A000_samplesheet.csv")
        hostName = None
        snakemakeBin = os.path.join(gitDir, "venv", "bin", "snakemake")
        nanoclid = NanoClid(
            inputFolder=inputFolder,
            bedDir=bedDir,
            run=run,
            outDir=outDir,
            dryRun=dryRun,
            genomeVersion=genomeVersion,
            until=until,
            samples=samples,
            outTemplate=outTemplate,
            snpEffDir=snpEffDir,
            copyToTransverse=copyToTransverse,
            copyToWorkspace=copyToWorkspace,
            sampleSheet=sampleSheet,
            refDir=refDir,
            hostName=hostName,
            snakemakeBin=snakemakeBin,
            email=email,
            profile=profile,
            containersFolder=containersFolder,
            doradoModelsFolder=doradoModelsFolder,
            basecaller=self.basecaller
        )
        nanoclid._runNanoClid()

if __name__ == "__main__":
    parser = ArgumentParser(
    description="Launch NanoClid pipeline")

    subs = parser.add_subparsers(dest = "command")

    install_parser = subs.add_parser("install", help='Install NanoClid images.')
    install_parser.add_argument("-S", "--singularityFolder", help="Where to download singularity images ? Default is GIT_DIR/singularity", default=os.path.join(os.path.dirname(os.path.realpath(__file__)), 'singularity'))

    test_parser = subs.add_parser("test", help="Run test")
    test_parser.add_argument("-e", "--email", help="Email adress to send run informations.")
    test_parser.add_argument("-p", "--profile", required=True, help="Profile to use to launch NanoCliD. Must be standalone|calcsub|abacus")
    test_parser.add_argument("-R", "--refDir", required=True, help="Path to folder containing genome files")
    test_parser.add_argument("-C", "--containersFolder", required=True, help="Path to singularity containers folder (required).")
    test_parser.add_argument("-M", "--doradoModelsFolder", help="Path to dorado models folder.", required=True)
    test_parser.add_argument("-X", "--basecaller", choices=["dorado","guppy"], default="dorado", help="Choose basecaller to use (dorado or guppy).")

    run_parser = subs.add_parser("run", help='Run NanoClid')
    run_parser.add_argument("-b", "--bedDir", help="Path to folder containing bed files.", default = "")
    run_parser.add_argument("--bedFile", help="Path to bed file.", default = "")
    run_parser.add_argument("-B", "--snakemakeBin", help="Path to snakemake bin.")
    run_parser.add_argument("-c", "--noCopyToTransverse", action='store_false', help="Do not copy results to transverse. Only for curie network.")
    run_parser.add_argument("-C", "--containersFolder", required=True, help="Path to singularity containers folder (required).")
    run_parser.add_argument("-D", "--snpEffDir", help="Folder containing snpEff annotation files.",  default=os.path.join(os.path.dirname(os.path.realpath(__file__)), 'annotations/data'))
    run_parser.add_argument("-e", "--email", help="Used by snakemake to send error/launching message.")
    run_parser.add_argument("-g", "--genomeVersion", required=True, help="Genome version for the analysis.")
    run_parser.add_argument("--hostName", help="Cluster on which launch the analysis. Only for curie network")
    run_parser.add_argument("-I", "--inputFolder", required=True, help="Path to the folder containing nanopore sequencing data.")
    run_parser.add_argument("-n", "--dryRun", action='store_true', help="Run nanoclid in dry run mode.")
    run_parser.add_argument("-O", "--outDir", help="Path to output folder.")
    run_parser.add_argument("-p", "--profile", help="Profile to use to launch NanoCliD. Must be standalone|calcsub|abacus")
    run_parser.add_argument("-q", "--queue", help="Queue to use for cluster launch")
    run_parser.add_argument("--noRunAllAnalysisOnCluster", action='store_false', help="Do not run on cluster. Only for curie network.")
    run_parser.add_argument("-r", "--runID", required = True, help="runID")
    run_parser.add_argument("-R", "--refDir", help="Ref dir containing genome annotations and reference files.", default="")
    run_parser.add_argument("-s", "--samples", help="Only run analysis on these samples. Must be comma separated.", default="")
    run_parser.add_argument("-S", "--samplesheet", help="Path to the samplesheet.")
    run_parser.add_argument("-T", "--outputTemplate", help="Path to the output template.", default=os.path.join(os.path.dirname(os.path.realpath(__file__)), 'templates', 'externe', 'ADAPTIVE.template'))
    run_parser.add_argument("-t", "--transverseFolder", help="Path to transverse folder.")
    run_parser.add_argument("-U", "--until", help="Specify until which rule you want to run the workflow", default="")
    run_parser.add_argument("-w", "--noCopyToWorkspace", action='store_false', help="Do not copy results to workspace. Only for curie network")
    run_parser.add_argument("-M", "--doradoModelsFolder", help="Path to dorado models folder.", required=True, default=None)
    run_parser.add_argument("-X", "--basecaller", choices=["dorado","guppy"], default="dorado", help="Choose basecaller to use (dorado or guppy).")

    args = parser.parse_args()

    if args.command == "install":
        NanoClid(profile = "standalone", containersFolder=args.singularityFolder)._install(args.singularityFolder)

    if args.command == "test":
        NanoClid(profile = "standalone", containersFolder=args.containersFolder, doradoModelsFolder=args.doradoModelsFolder, basecaller=args.basecaller)._runTest(args.refDir, args.profile, args.email, args.containersFolder, args.doradoModelsFolder)

    if args.command == "run":
        
        if not args.outDir:
            args.outDir = args.inputFolder
        try:
            nanoclid = NanoClid(args.inputFolder, 
                                args.bedDir,
                                args.bedFile,
                                args.runID,
                                args.outDir,
                                args.dryRun,
                                args.genomeVersion, 
                                args.until,
                                args.samples,
                                args.outputTemplate,
                                args.snpEffDir,
                                args.noCopyToTransverse,
                                args.noCopyToWorkspace,
                                args.noRunAllAnalysisOnCluster,
                                args.transverseFolder,
                                args.samplesheet,
                                args.refDir,
                                args.hostName,
                                args.snakemakeBin,
                                args.email,
                                args.profile,
                                args.queue,
                                args.containersFolder,
                                args.doradoModelsFolder,
                                args.basecaller
                                )
            
            
            nanoclid._runNanoClid()
        except Exception as e:
            with open(os.path.join(nanoclid.outDir, nanoclid.run, 'errorLaunching.txt'), 'w') as f:
                f.write(str(e))
                f.write(traceback.format_exc())
            profile = nanoclid.loadConfig(os.path.join(nanoclid.profile, "config.yaml"))
         #   containersPath = profile["singularity-prefix"]
            containersPath = args.containersFolder
            config = nanoclid.loadConfig(nanoclid.configTemplate)
            config["email"] = "f.oberhammer+debug@prinsesmaximacentrum.nl"
            config["errorMail"]["content"] = traceback.format_exc()
            config["errorMail"]["subject"] = config["errorMail"]["subject"][1:-1] #remove ''
            from utils.utils import sendMail
            sendMail(config, containersPath, curieNetwork, "onerror", os.path.join(nanoclid.outDir, nanoclid.run, 'errorLaunching.txt'))