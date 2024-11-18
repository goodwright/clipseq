#!/usr/bin/env python

import os
import sys
import errno
import argparse
import platform

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


def dump_versions(process_name):
    with open("versions.yml", "w") as out_f:
        out_f.write(process_name + ":\n")
        out_f.write("    python: " + platform.python_version() + "\n")


def check_samplesheet(process_name, samplesheet, output):
    """
    This function checks that the samplesheet follows the following structure:

    group,replicate,control_group,control_replicate
    WT_BRD4_plus,1,WT_No_Ab,1
    WT_BRD4_plus,2,WT_No_Ab,2
    WT_FUS_plus,1,WT_No_Ab,1
    WT_FUS_plus,2,WT_No_Ab,2
    """

    # Dump version file
    dump_versions(process_name)

    rows = []
    samples = []

    with open(samplesheet, "r") as fin:
        ## Check header
        MIN_COLS = 4
        HEADER = ["group", "replicate", "control_group", "control_replicate"]
        # HEADER_LEN = len(HEADER)
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        ACTUAL_HEADER_LEN = len(header)

        if header[: len(HEADER)] != HEADER:
            print("ERROR: Please check samplesheet header -> {} != {}".format(",".join(header), ",".join(HEADER)))
            sys.exit(1)

        ## Check sample entries
        line_no = 1
        for line in fin:
            lspl = [x.strip().strip('"') for x in line.strip().split(",")]

            ## Check if its just a blank line so we dont error
            if line.strip() == "":
                continue

            ## Check valid number of columns per row
            if len(lspl) != ACTUAL_HEADER_LEN:
                print_error(
                    "Invalid number of columns (found {} should be {})! - line no. {}".format(
                        len(lspl), len(HEADER), line_no
                    ),
                    "Line",
                    line,
                )

            ## Check valid number of populated columns per row
            num_cols = len([x for x in lspl if x])
            if num_cols < MIN_COLS:
                print_error(
                    "Invalid number of populated columns (minimum = {})!".format(MIN_COLS),
                    "Line",
                    line,
                )

            ## Check sample name entries
            group, replicate, control_group, control_replicate = lspl[: len(HEADER)]
            if group:
                if group.find(" ") != -1:
                    print_error("Group entry contains spaces!", "Line", line)
            else:
                print_error("Group entry has not been specified!", "Line", line)

            if group:
                if group.find(".") != -1:
                    print_error("Group entry contains dots!", "Line", line)

            if control_group:
                if control_group.find(" ") != -1:
                    print_error("Control Group entry contains spaces!", "Line", line)
            else:
                print_error("Control Group entry has not been specified!", "Line", line)
            
            if control_group:
                if control_group.find(".") != -1:
                    print_error("Control Group entry contains dots!", "Line", line)

            ## Check replicate entry is integer
            if not replicate.isdigit():
                print_error("Replicate id not an integer", "Line", line)
            replicate = int(replicate)
            if replicate <= 0:
                print_error("Replicate must be > 0", "Line", line)
            
            ## Check control_replicate entry is integer
            if not control_replicate.isdigit():
                print_error("Control Replicate id not an integer", "Line", line)
            control_replicate = int(control_replicate)
            if control_replicate <= 0:
                print_error("Control Replicate must be > 0", "Line", line)
            
            ## Check replicate and control_replicate are not the same
            if replicate == control_replicate:
                print_error("Replicate and Control Replicate must be different", "Line", line)

            ## Check for duplicate entries
            sample_info = (group, replicate, control_group, control_replicate)
            if sample_info in rows:
                print_error("Samplesheet contains duplicate rows!", "Line", line)
            else:
                rows.append(sample_info)

            ## Create sample id and control_sample id
            ## This id matches with the final id in meta (without _Txxx)
            sample_id = "{}_R{}".format(group, replicate)
            control_sample_id = "{}_R{}".format(control_group, control_replicate)

            # Save the sample id and control sample id in front of the row
            lspl = [sample_id, control_sample_id] + lspl
            samples.append(lspl)

            line_no += 1
        
        # Write out the samplesheet in csv with header
        make_dir(os.path.dirname(output))
        with open(output, "w") as fout:
            output_header = ["id", "control_id"] + HEADER
            fout.write(",".join(output_header) + "\n")
            for sample in samples:
                fout.write(",".join(sample) + "\n")


if __name__ == "__main__":
    # Allows switching between nextflow templating and standalone python running using arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--process_name", default="!{process_name}")
    parser.add_argument("--samplesheet", default="!{samplesheet}")
    parser.add_argument("--output", default="!{output}")
    args = parser.parse_args()

    check_samplesheet(args.process_name, args.samplesheet, args.output)