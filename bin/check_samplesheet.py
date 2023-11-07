#!/usr/bin/env python3

"""
adapted from https://github.com/nf-core/chipseq/blob/51eba00b32885c4d0bec60db3cb0a45eb61e34c5/bin/check_samplesheet.py
python check_samplesheet.py /data/CCBR/projects/ccbr1301/analysis/cruise_231106/assets/sample_manifest.csv /data/CCBR/projects/ccbr1301/analysis/cruise_231106/assets/contrast_manifest.csv 
"""

import os
import sys
import errno
import argparse


def parse_args(args=None):
    Description = "Reformat samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN_SAMPLESHEET> <FILE_IN_CONTRASTSHEET>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN_SAMPLESHEET", help="Input samplesheet file.")
    parser.add_argument("FILE_IN_CONTRASTSHEET", help="Input contrastsheet file.")
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)

def check_samplesheet(file_in_samplesheet, file_in_contrastsheet):
    """
    This function checks that the samplesheet follows the following structure:
    sample,fastq_1,fastq_2,treat_or_ctrl,groupID
    SPT5_T0_REP1,SRR1822153_1.fastq.gz,SRR1822153_2.fastq.gz,treatment,group1
    SPT5_T0_REP2,SRR1822154_1.fastq.gz,SRR1822154_2.fastq.gz,control,group1
    SPT6_T0_REP1,SRR1822155_1.fastq.gz,SRR1822155_2.fastq.gz,treatment,group2
    SPT6_T0_REP2,SRR1822156_1.fastq.gz,SRR1822156_2.fastq.gz,control,group2

    This function checks that the contrastsheet follows the following structure:
    groupID,libraryID
    group1,lib-01
    group2,lib-02
    """
    #####################################################################################
    # check samplesheet, create mapping dict
    sample_mapping_dict = {}
    with open(file_in_samplesheet, "r", encoding="utf-8-sig") as fin:
        ## Check header
        MIN_COLS = 2
        HEADER = ["sample", "fastq_1", "fastq_2", "treat_or_ctrl", "groupID"]
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        if header[: len(HEADER)] != HEADER:
            print(
                f"ERROR: Please check samplesheet header -> {','.join(header)} != {','.join(HEADER)}"
            )
            sys.exit(1)

        ## Check sample entries
        for line in fin:
            lspl = [x.strip().strip('"') for x in line.strip().split(",")]

             # If it's a blank line, next
            if len(line.strip()) == 0:
                print("Skipping blank line")
                continue

            # Check valid number of columns per row
            if len(lspl) < len(HEADER):
                print_error(
                    "Invalid number of columns (minimum = {})!".format(len(HEADER)),
                    "Line",
                    line,
                )
            num_cols = len([x for x in lspl if x])
            if num_cols < MIN_COLS:
                print_error(
                    "Invalid number of populated columns (minimum = {})!".format(
                        MIN_COLS
                    ),
                    "Line",
                    line,
                )

            ## Check sample name entries
            sample, fastq_1, fastq_2, treat_or_ctrl, groupID = lspl[: len(HEADER)]
            if sample.find(" ") != -1:
                print(
                    f"WARNING: Spaces have been replaced by underscores for sample: {sample}"
                )
                sample = sample.replace(" ", "_")
            if not sample:
                print_error("Sample entry has not been specified!", "Line", line)

            ## Check FastQ file extension, exists
            for fastq in [fastq_1, fastq_2]:
                if fastq:
                    if fastq.find(" ") != -1:
                        print_error("FastQ file contains spaces!", "Line", line)
                    if not fastq.endswith(".fastq.gz") and not fastq.endswith(".fq.gz"):
                        print_error(
                            "FastQ file does not have extension '.fastq.gz' or '.fq.gz'!",
                            "Line",
                            line,
                        )

            ## check treatment or control designation
            if treat_or_ctrl not in {"treatment", "control"}:
                print_error(
                    "treat_or_ctrl column can only contain values `treatment` or `control`.",
                    "Line",
                    line,
                )

            ## Auto-detect paired-end/single-end
            if sample and fastq_1 and fastq_2:  ## Paired-end short reads
                is_single = "0"
            elif sample and fastq_1 and not fastq_2:  ## Single-end short reads
                is_single = "1"
            else:
                print_error("Invalid combination of columns provided!", "Line", line)

            ## Create sample mapping dictionary = {group: [[ sample: [[ single_end, fastq_1, fastq_2,]] ]]}
            if groupID not in sample_mapping_dict:
                    sample_mapping_dict[groupID] = {}

            sample_info = [is_single, fastq_1, fastq_2, treat_or_ctrl]

            if sample not in sample_mapping_dict[groupID]:
                sample_mapping_dict[groupID][sample] = [sample_info]
            else:
                if sample_info in sample_mapping_dict[groupID][sample]:
                    print_error("Samplesheet contains duplicate rows!", "Line", line)
                else:
                    sample_mapping_dict[groupID][sample].append(sample_info)

    ###################################################################################################
    # check contrastsheet
    contrast_mapping_dict = {}
    with open(file_in_contrastsheet, "r") as fin:
        ## Check header
        MIN_COLS = 2
        HEADER = ["groupID", "libraryID"]
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        if header[: len(HEADER)] != HEADER:
            print(
                f"ERROR: Please check contrastsheet header -> {','.join(header)} != {','.join(HEADER)}"
            )
            sys.exit(1)

        ## Check sample entries
        for line in fin:
            lspl = [x.strip().strip('"') for x in line.strip().split(",")]

            # Check valid number of columns per row
            num_cols = len([x for x in lspl if x])
            if num_cols < MIN_COLS:
                print_error(
                    "Invalid number of populated columns (minimum = {})!".format(MIN_COLS),
                    "Line",
                    line,
                )
            
            # validate groupID is in samplesheet
            groupID, libraryID = lspl[: len(HEADER)]
            if groupID not in sample_mapping_dict:
                print_error("groupID entry found in the contrastsheet is not one listed in the samplesheet. Check groupIDs!", "Line", groupID)

    ###################################################################################################
    ## Write validated samplesheet(s) with appropriate columns
    if len(sample_mapping_dict) > 0:
        for groupID in sorted(sample_mapping_dict.keys()):
            file_out=groupID + ".csv"
            with open(file_out, "w") as fout:
                fout.write(
                    ",".join(
                        ["sample", "single_end", "fastq_1", "fastq_2", "treat_or_ctrl"]
                    )
                    + "\n"
                )
                for sample in sorted(sample_mapping_dict[groupID].keys()):
                    ## Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
                    if not all(
                        x[0] == sample_mapping_dict[groupID][sample][0][0]
                        for x in sample_mapping_dict[groupID][sample]
                    ):
                        print_error(
                            f"Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end!",
                            "Sample",
                            sample,
                        )

                    for idx, val in enumerate(sample_mapping_dict[groupID][sample]):
                        plus_T = (
                            f"_T{idx+1}" if len(sample_mapping_dict[groupID][sample]) > 1 else ""
                        )  # do not append _T{idx} if not needed
                        fout.write(",".join([f"{sample}{plus_T}"] + val) + "\n")
    else:
        print_error(f"No entries to process!", "Samplesheet: {file_in_samplesheet}")


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN_SAMPLESHEET, args.FILE_IN_CONTRASTSHEET)


if __name__ == "__main__":
    sys.exit(main())
