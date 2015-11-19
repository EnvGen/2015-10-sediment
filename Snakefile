__author__ = "Johannes Alneberg"
__license__ = "MIT"


import os
import sys
import shutil
import glob
from subprocess import check_output

# Check that no submodule git repo is dirty
#submodules = ["BLUEPRINT_pipeline"]
#for submodule in submodules:
#    submodule_status = check_output(["git", "status", "--porcelain", submodule])
#    if not submodule_status == b"":
#        print(submodule_status)
#        raise Exception("Submodule {} is dirty. Commit changes before proceeding.".format(submodule))

# Check that the git repo is not dirty
#submodule_status = check_output(["git", "status", "--porcelain"])
#if not submodule_status == b"":
#    print(submodule_status)
#    raise Exception("Repo is dirty. Commit changes before proceeding.")

# Chose config file based on if we're on uppmax or not
configfile: "config.json"

config["fastqc_rules"]["reads"] = {}
config["cutadapt_rules"]["reads"] = {}
config["megahit_rules"]["samples"] = {}

with open("sample_indices.json") as si:
    sample_indices = json.load(si)

for read_file in glob.glob("samples/raw/*.fq.gz"):
    read_basename = os.path.basename(read_file)
    read_name = read_basename.replace(".fq.gz", "")
    config["fastqc_rules"]["reads"][read_name] = read_file
   
    # Add all steps to fastqc - this will cause fastqc to run after each step 
    # as well as on the raw reads 
    for trim_params_name, trim_params_dict in config["cutadapt_rules"]["trim_params"].items():
        config["fastqc_rules"]["reads"]["cutadapt_"+trim_params_name+"_"+read_name] = \
            "cutadapt/adapt_cutting/{trim_params}/{read}".format(
                trim_params=trim_params_name,
                read = read_basename
            ) 

        config["fastqc_rules"]["reads"]["fastuniq_"+trim_params_name+"_"+read_name] = \
            "fastuniq/{trim_params}/{read}".format(
                trim_params=trim_params_name,
                read = read_basename
            ) 
    # Hack to get read pairs in a list
    read_name = read_name.replace("_R1", "").replace("_R2", "")

    if read_name in config["cutadapt_rules"]["reads"]:
        config["cutadapt_rules"]["reads"][read_name].append(read_file)
        config["cutadapt_rules"]["reads"][read_name].sort()
    else:
        config["cutadapt_rules"]["reads"][read_name] = [read_file]
    
    # Add the variable barcode sequences for each sample to cutadapt config
    if read_name in sample_indices["R1_index"]:
        for trim_params_name, trim_params_config in config["cutadapt_rules"]["trim_params"].items():
            if "common_variables" in trim_params_config.keys():
                variables = config["cutadapt_rules"]["trim_params"][trim_params_name]["common_variables"].copy()
                variables["R1_index"] = sample_indices["R1_index"][read_name]
                if "R2_rev_index" in sample_indices:
                    variables["R2_rev_index"] = sample_indices["R2_rev_index"][read_name]
                else:
                    variables["R2_rev_index"] = None
                if "variables" not in config["cutadapt_rules"]["trim_params"][trim_params_name]:
                    config["cutadapt_rules"]["trim_params"][trim_params_name]["variables"] = {}
                config["cutadapt_rules"]["trim_params"][trim_params_name]["variables"][read_name] = variables

for read_file in glob.glob("finished_reads/*.fq.gz"):
    read_basename = os.path.basename(read_file)
    read_name = read_basename.replace(".fq.gz", "")
    sample_name = read_name.replace("_R1", "").replace("_R2", "")
    
    if sample_name in config["megahit_rules"]["samples"]:
        config["megahit_rules"]["samples"][sample_name].append(read_file)
        config["megahit_rules"]["samples"][sample_name].sort()
    else:
        config["megahit_rules"]["samples"][sample_name] = [read_file]

for contigs_f in glob.glob("assembly/megahit_coassembly/default/parts/contigs.*.fasta"):
    part = contigs_f.split('.')[-2]
    config["prokka_extended_rules"]["contigs"]['megahit_coassembly.{}'.format(part)] = contigs_f

WORKFLOW_DIR = "snakemake-workflows/"

#include: os.path.join(WORKFLOW_DIR, "rules/mapping/bowtie2.rules")
#include: os.path.join(WORKFLOW_DIR, "rules/mapping/samfiles.rules")
#include: os.path.join(WORKFLOW_DIR, "rules/quantification/rpkm.rules")
include: os.path.join(WORKFLOW_DIR, "bio/ngs/rules/trimming/cutadapt.rules")
include: os.path.join(WORKFLOW_DIR, "bio/ngs/rules/quality_control/fastqc.rules")
include: os.path.join(WORKFLOW_DIR, "bio/ngs/rules/duplicate_removal/fastuniq.rules")
include: os.path.join(WORKFLOW_DIR, "bio/ngs/rules/assembly/megahit.rules")
include: os.path.join(WORKFLOW_DIR, "bio/ngs/rules/annotation/prokka.rules")

rule preprocess_all:
    input:
        htmls=expand("fastqc/{reads}/{reads}_fastqc.html", reads=config["fastqc_rules"]["reads"]),
        zips=expand("fastqc/{reads}/{reads}_fastqc.zip", reads=config["fastqc_rules"]["reads"])


# Testing rules and configs
test_reads_orig = {
            "P2237_101_R1": "samples/raw/P2237_101_R1.fq.gz",
            "P2237_101_R2": "samples/raw/P2237_101_R2.fq.gz",
            "P2237_102_R1": "samples/raw/P2237_102_R1.fq.gz",
            "P2237_102_R2": "samples/raw/P2237_102_R2.fq.gz",
            "P2237_111_R1": "samples/raw/P2237_110_R1.fq.gz",
            "P2237_111_R2": "samples/raw/P2237_110_R2.fq.gz"
        }

test_reads = {}
for read_name, read_file in test_reads_orig.items():
    test_reads[read_name] = read_file
    read_basename = os.path.basename(read_file)
    for trim_params_name, trim_params_dict in config["cutadapt_rules"]["trim_params"].items():
        test_reads["cutadapt_" + trim_params_name + "_" + read_name] = \ 
            "cutadapt/adapt_cutting/{trim_params}/{read}".format(
                trim_params=trim_params_name,
                read = read_basename
            )

        test_reads["fastuniq_"+trim_params_name+"_"+read_name] = \
            "fastuniq/{trim_params}/{read}".format(
                trim_params=trim_params_name,
                read = read_basename
            ) 

rule fastqc_all_test:
    input:
        htmls=expand("fastqc/{reads}/{reads}_fastqc.html", reads=test_reads),
        zips=expand("fastqc/{reads}/{reads}_fastqc.zip", reads=test_reads)

