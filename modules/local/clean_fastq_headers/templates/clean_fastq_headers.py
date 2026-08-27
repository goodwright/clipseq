#!/usr/bin/env python3

import platform
import subprocess


def process_fastq(input_file, output_file, threads):
    """Process FASTQ file to remove spaces from headers using pigz."""
    is_gzipped = input_file.endswith('.gz')

    in_file = None
    if is_gzipped:
        decompress = subprocess.Popen(['pigz', '-p', str(threads), '-dc', input_file],
                                      stdout=subprocess.PIPE)
    else:
        in_file = open(input_file, 'rb')
        decompress = subprocess.Popen(['cat'], stdin=in_file, stdout=subprocess.PIPE)

    with open(output_file, 'wb') as out_handle:
        compress = subprocess.Popen(['pigz', '-p', str(threads), '-c'],
                                    stdin=subprocess.PIPE, stdout=out_handle)
        try:
            line_count = 0
            for line in decompress.stdout:
                if line_count % 4 == 0:  # Header line
                    line = line.decode().strip().replace(' ', '_').encode() + b'\n'
                compress.stdin.write(line)
                line_count += 1
        finally:
            decompress.stdout.close()
            decompress.wait()
            if in_file is not None:
                in_file.close()
            compress.stdin.close()
            compress.wait()

    # Both ends of the pipe must be checked. pigz reports a corrupt or
    # truncated input by exiting non-zero after emitting a partial stream, so
    # without this the task "succeeds" and hands a silently truncated FASTQ to
    # the rest of the pipeline.
    if decompress.returncode != 0:
        raise RuntimeError(
            "reading %s failed (exit %d) - the file is likely truncated or corrupt"
            % (input_file, decompress.returncode))
    if compress.returncode != 0:
        raise RuntimeError(
            "writing %s failed (exit %d)" % (output_file, compress.returncode))

    if line_count % 4 != 0:
        raise RuntimeError(
            "%s ended mid-record after %d lines - the file is likely truncated"
            % (input_file, line_count))


def pigz_version():
    """pigz reports its version on stderr on some builds, stdout on others."""
    proc = subprocess.run(['pigz', '--version'], stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT)
    return proc.stdout.decode().strip().replace('pigz ', '')


# This pipeline is single-end only: single_end is parsed into meta but never
# read downstream, and each sample is one FASTQ. Fail loudly rather than
# quietly emit files the rest of the pipeline cannot use.
reads = "!{reads}".replace('[', '').replace(']', '').replace(',', ' ').split()
if len(reads) != 1:
    raise SystemExit(
        "CLEAN_FASTQ_HEADERS expects exactly one FASTQ per sample, got %d (%s)"
        % (len(reads), ", ".join(reads)))

process_fastq(reads[0], "!{prefix}.clean.fastq.gz", !{task.cpus})

with open("versions.yml", "w") as out_f:
    out_f.write("!{process_name}" + ":\n")
    out_f.write("    python: " + platform.python_version() + "\n")
    out_f.write("    pigz: " + pigz_version() + "\n")
