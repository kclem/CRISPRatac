"""
CRISPRatac aligns ATAC-sequencing reads to a genome and to regions produced by combination of CRISPR cut sites

Program requirements:
samtools
bowtie2
matplotlib
CRISPResso2

Parameters may be passed in via the command line, or a settings file.

To run, create the file 'settings.txt'
=== settings.txt ====
fastq_r1	read1.fq.gz
fastq_r2	read2.fq.gz
genome	/data/pinello/COMMON_DATA/REFERENCE_GENOMES/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome.fa
cuts	chr1:1111111 chr2:222222
=====================

Then run 'python CRISPRatac.py settings.txt'

"""
import argparse
from collections import defaultdict
import json
import logging
import gzip
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import multiprocessing as mp
import os
import re
import subprocess
import sys
from CRISPResso2 import CRISPRessoMultiProcessing
from CRISPResso2 import CRISPRessoShared

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

__version__ = "v0.0.1"

def main():
    settings, logger = parse_settings(sys.argv)

    #data structures for plots for report
    summary_plot_objects=[]  # list of PlotObjects for plotting

    assert_dependencies(
            samtools_command=settings['samtools_command'],
            bowtie2_command=settings['bowtie2_command'],
            crispresso_command=settings['crispresso_command'],
            )

    av_read_length = get_av_read_len(settings['fastq_r1'])
    num_reads_input = get_num_reads_fastq(settings['fastq_r1'])
    logger.info('%d reads in input'%num_reads_input)

    cut_sites = settings['cuts']
    cut_annotations = {}
    for cut_site in settings['cuts']:
        cut_annotations[cut_site] = 'Cut site'

    target_padding = settings['alignment_extension']
    target_length = av_read_length
    if target_padding < 0:
        target_length += target_padding
        target_padding = 0
        target_padding += av_read_length

    if len(cut_sites) > 0:
        (custom_index_fasta,target_names,target_info) = make_artificial_targets(
                root = settings['root']+'.customTargets',
                cuts=cut_sites,
                cut_annotations=cut_annotations,
                genome=settings['genome'],
                target_length=target_length,
                target_padding=target_padding,
                samtools_command=settings['samtools_command'],
                bowtie2_command=settings['bowtie2_command'],
                bowtie2_threads = settings['n_processes'],
                can_use_previous_analysis = settings['can_use_previous_analysis']
                )
    else:
        custom_index_fasta = None
        target_names = []
        target_info = {}


    curr_r1_file = settings['fastq_r1'] #if alignment to genome happens first, the input for artificial target mapping will be reads that don't align to the genome
    curr_r1_description = "input R1 reads"
    curr_r2_file = settings['fastq_r2']
    curr_r2_description = "input R2 reads"

    genome_aligned_count = 0 # number of reads aligned to genome
    custom_aligned_count = 0 # number of reads aligned to custom targets

    (genome_r1_assignments,genome_r2_assignments,genome_unmapped_r1, genome_unmapped_r2, genome_aligned_count, genome_mapped_bam_file,genome_chr_aln_plot_obj,genome_tlen_plot_object
        ) = align_reads(
                root = settings['root']+'.genomeAlignment',
                fastq_r1 = curr_r1_file,
                fastq_r2 = curr_r2_file,
                reads_name = curr_r1_description.replace(" R1",""),
                bowtie2_reference = settings['bowtie2_genome'],
                reference_name = 'Genome',
                target_info = target_info,
                bowtie2_command = settings['bowtie2_command'],
                bowtie2_threads = settings['n_processes'],
                samtools_command = settings['samtools_command'],
                ignore_n = settings['ignore_n'],
                keep_intermediate = settings['keep_intermediate'],
                can_use_previous_analysis = settings['can_use_previous_analysis']
                )

    if genome_tlen_plot_object is not None:
        genome_tlen_plot_object.order = 2
        summary_plot_objects.append(genome_tlen_plot_object)
    if genome_chr_aln_plot_obj is not None:
        genome_chr_aln_plot_obj.order = 5
        summary_plot_objects.append(genome_chr_aln_plot_obj)

    curr_r1_file = genome_unmapped_r1
    curr_r1_description = "R1 reads not aligned to genome"
    curr_r2_file = genome_unmapped_r2
    curr_r2_description = "R2 reads not aligned to genome"

    crispresso_run_names = []
    crispresso_sub_htmls = {}
    if len(target_names) > 0:
        #first, align the unaligned r1s
        (custom_r1_assignments,none_file,custom_unmapped_r1, none_file2, custom_r1_aligned_count, custom_r1_mapped_bam_file,custom_r1_chr_aln_plot_obj,none_custom_tlen_plot_object
            ) = align_reads(
                    root = settings['root']+'.r1.customTargetAlignment',
                    fastq_r1 = curr_r1_file,
                    fastq_r2 = None,
                    reads_name = curr_r1_description,
                    bowtie2_reference=custom_index_fasta,
                    reference_name='Custom targets',
                    target_info = target_info,
                    arm_min_seen_bases = settings['arm_min_seen_bases'],
                    arm_min_matched_start_bases = settings['arm_min_matched_start_bases'],
                    bowtie2_command=settings['bowtie2_command'],
                    bowtie2_threads=settings['n_processes'],
                    samtools_command=settings['samtools_command'],
                    ignore_n = settings['ignore_n'],
                    keep_intermediate = settings['keep_intermediate'],
                    can_use_previous_analysis = settings['can_use_previous_analysis']
                    )

        crispresso_run_names_r1, crispresso_sub_htmls_r1 = run_crispresso2(
                    root = settings['root']+'.CRISPResso_r1',
                    input_fastq_file = curr_r1_file,
                    assignment_file = custom_r1_assignments,
                    target_info = target_info,
                    crispresso_name_suffix = 'R1',
                    crispresso_min_count = settings['crispresso_min_count'],
                    crispresso_min_aln_score = settings['crispresso_min_aln_score'],
                    crispresso_quant_window_size = settings['crispresso_quant_window_size'],
                    crispresso_command=settings['crispresso_command'],
                    n_processes = settings['n_processes']
                    )
        for name in crispresso_run_names_r1:
            crispresso_run_names.append(name)
            crispresso_sub_htmls[name] = crispresso_sub_htmls_r1[name]

        custom_aligned_count += custom_r1_aligned_count
        curr_r1_file = custom_unmapped_r1
        curr_r1_description = "R1 reads not aligned to custom targets"

        if custom_r1_chr_aln_plot_obj is not None:
            custom_r1_chr_aln_plot_obj.order = 10
            summary_plot_objects.append(custom_r1_chr_aln_plot_obj)

        #if input is paired, align second reads to the custom targets (and genome as well if not performed previously)
        if curr_r2_file is not None:
            (custom_r2_assignments,none_file,custom_unmapped_r2, none_file2, custom_r2_aligned_count, custom_r2_mapped_bam_file,custom_r2_chr_aln_plot_obj,none_custom_tlen_plot_object
                ) = align_reads(
                        root = settings['root']+'.r2.customTargetAlignment',
                        fastq_r1 = curr_r2_file,
                        fastq_r2 = None,
                        reads_name = curr_r2_description,
                        bowtie2_reference=custom_index_fasta,
                        reference_name='Custom targets',
                        target_info = target_info,
                        arm_min_seen_bases = settings['arm_min_seen_bases'],
                        arm_min_matched_start_bases = settings['arm_min_matched_start_bases'],
                        bowtie2_command=settings['bowtie2_command'],
                        bowtie2_threads=settings['n_processes'],
                        samtools_command=settings['samtools_command'],
                        ignore_n = settings['ignore_n'],
                        keep_intermediate = settings['keep_intermediate'],
                        can_use_previous_analysis = settings['can_use_previous_analysis']
                        )
            custom_aligned_count += custom_r2_aligned_count
            crispresso_run_names_r2, crispresso_sub_htmls_r2 = run_crispresso2(
                        root = settings['root']+'.CRISPResso_r2',
                        input_fastq_file = curr_r2_file,
                        assignment_file = custom_r2_assignments,
                        target_info = target_info,
                        crispresso_name_suffix = 'R2',
                        crispresso_min_count = settings['crispresso_min_count'],
                        crispresso_min_aln_score = settings['crispresso_min_aln_score'],
                        crispresso_quant_window_size = settings['crispresso_quant_window_size'],
                        crispresso_command=settings['crispresso_command'],
                        n_processes = settings['n_processes']
                        )

        for name in crispresso_run_names_r2:
            crispresso_run_names.append(name)
            crispresso_sub_htmls[name] = crispresso_sub_htmls_r2[name]

            if custom_r2_chr_aln_plot_obj is not None:
                custom_r2_chr_aln_plot_obj.order = 11
                summary_plot_objects.append(custom_r2_chr_aln_plot_obj)

            curr_r2_file = custom_unmapped_r2
            curr_r2_description = "R2 reads not aligned to custom targets"

    
    labels = ["Input Reads","Aligned Custom Targets","Aligned Genome"]
    values = [num_reads_input,custom_aligned_count,genome_aligned_count]
    alignment_summary_root = settings['root']+".alignmentSummary"
    with open(alignment_summary_root+".txt",'w') as summary:
        summary.write("\t".join(labels)+"\n")
        summary.write("\t".join([str(x) for x in values])+"\n")
    fig = plt.figure(figsize=(12,12))
    ax = plt.subplot(111)
    pie_values = []
    pie_labels = []
    for i in range(1,len(labels)):
        if values[i] > 0:
            pie_values.append(values[i])
            pie_labels.append(labels[i]+"\n("+str(values[i])+")")
    ax.pie(pie_values, labels=pie_labels, autopct="%1.2f%%")
    plot_name = alignment_summary_root
    plt.savefig(plot_name+".pdf",pad_inches=1,bbox_inches='tight')
    plt.savefig(plot_name+".png",pad_inches=1,bbox_inches='tight')
    plot_count_str = "<br>".join(["%s N=%s"%x for x in zip(labels,values)])


    summary_plot_objects.append(
            PlotObject(plot_name = plot_name,
                    plot_title = 'Alignment Summary',
                    plot_label = 'Pie chart showing final assignment of reads.<br>' + plot_count_str,
                    plot_datas = [('Alignment summary',alignment_summary_root + ".txt")]
                    ))

    make_report(report_file=settings['root']+".html",
            report_name = 'Report',
            crispratac_folder = '',
            crispresso_run_names = crispresso_run_names,
            crispresso_sub_html_files = crispresso_sub_htmls,
            summary_plot_objects = summary_plot_objects,
            )


    logger.info('Successfully completed!')

    # FINISHED

def parse_settings(args):
    """
    Parses settings from the command line
        First parses from the settings file, then parses from command line

    param:
        args: command line arguments

    returns:
        settings: dict of parsed settings
    """
    parser = argparse.ArgumentParser(description='CRISPRatac: tools for aligning ATAC-seq reads to CRISPR-cleaved genomes', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--version', action='version', version="%(prog)s "+__version__)
    parser.add_argument('settings_file', nargs='*', help='Tab-separated settings file')

    parser.add_argument('--debug', action='store_true', help='Tab-separated settings file')
    parser.add_argument('--root', type=str, default=None, help='Output directory file root')
    parser.add_argument('--keep_intermediate',action='store_true',help='If true, intermediate files are not deleted')

    parser.add_argument('--cuts','--cut_sites', nargs='*', help='Cut sites in the form chr1:234 (multiple cuts are separated by spaces)', default=[])

    parser.add_argument('--genome', help='Genome sequence file for alignment. This should point to a file ending in ".fa", and the accompanying index file (".fai") should exist.', default=None)
    parser.add_argument('--bowtie2_genome', help='Bowtie2-indexed genome file.',default=None)

    parser.add_argument('--fastq_r1', help='Input fastq r1 file', default=None)
    parser.add_argument('--fastq_r2', help='Input fastq r2 file', default=None)

    custom_group = parser.add_argument_group('Custom target settings')
    custom_group.add_argument('--alignment_extension', type=int, help='Number of bp to extend beyond av read length around cut site for custom index', default=50)

    #min alignment cutoffs for alignment to each arm/side of read
    a_group = parser.add_argument_group('Alignment cutoff parameters')
    a_group.add_argument('--arm_min_seen_bases', type=int, help='Number of bases that are required to be seen on each "side" of translocated reads. E.g. if a artificial target represents a translocation between chr1 and chr2, arm_min_seen_bases would have to be seen on chr1 as well as on chr2 for the read to be counted.', default=15)
    a_group.add_argument('--arm_min_matched_start_bases', type=int, help='Number of bases that are required to be matching (no indels or mismatches) at the beginning of the read on each "side" of artifical targets. E.g. if a artificial target represents a translocation between chr1 and chr2, the first arm_min_matched_start_bases of the read would have to match exactly to chr1 and the last arm_min_matched_start_bases of the read would have to match exactly to chr2', default=10)
    a_group.add_argument('--ignore_n', type=bool, help='If set, "N" bases will be ignored. By default (False) N bases will count as mismatches in the number of bases required to match at each arm/side of the read', default=False)

    #CRISPResso settings
    c_group = parser.add_argument_group('CRISPResso settings')
    c_group.add_argument('--crispresso_min_count', type=int, help='Min number of reads required to be seen at a site for it to be analyzed by CRISPResso', default=50)
    c_group.add_argument('--crispresso_min_aln_score', type=int, help='Min alignment score to reference sequence for quantification by CRISPResso', default=20)
    c_group.add_argument('--crispresso_quant_window_size', type=int, help='Number of bp on each side of a cut to consider for edits', default=1)

    #sub-command parameters
    p_group = parser.add_argument_group('Pipeline parameters')
    p_group.add_argument('--samtools_command', help='Command to run samtools', default='samtools')
    p_group.add_argument('--bowtie2_command', help='Command to run bowtie2', default='bowtie2')
    p_group.add_argument('--crispresso_command', help='Command to run CRISPResso2', default='CRISPResso')
    p_group.add_argument('--n_processes', type=str, help='Number of processes to run on (may be set to "max")', default='1')

    cmd_args = parser.parse_args(args[1:])

    settings = {}
    settings_file_args = {}

    # try to read in args from the settings file
    if (len(cmd_args.settings_file) > 0):
        for s_file in cmd_args.settings_file:
            with open(s_file, 'r') as fin:
                for line in fin:
                    if line.startswith("#"):
                        continue
                    line_els = line.split("#")[0].rstrip('\n').split("\t")
                    if line_els == ['']:
                        continue
                    if (len(line_els) < 2):
                        raise Exception('Cannot parse line "' + line + '"\nA tab must separate the key and value in the settings file')
                    key = line_els[0].strip()
                    val = line_els[1].strip()
                    settings_file_args[key] = val

    settings['root'] = cmd_args.root
    if 'root' in settings_file_args:
        settings['root'] = settings_file_args['root']
        settings_file_args.pop('root')
    if cmd_args.root is None:
        if len(cmd_args.settings_file) > 0:
            settings['root'] = cmd_args.settings_file[0] + ".CRISPRatac"
        else:
            settings['root'] = "CRISPRatac"

    settings['debug'] = cmd_args.debug
    if 'debug' in settings_file_args:
        settings['debug'] = (settings_file_args['debug'].lower() == 'true')
        settings_file_args.pop('debug')

    logger = logging.getLogger('CRISPRatac')
    logging_level = logging.INFO
    if settings['debug']:
        logging_level=logging.DEBUG

    logger.setLevel(logging.DEBUG)

    log_formatter = logging.Formatter("%(asctime)s:%(levelname)s: %(message)s","%Y-%m-%d %H:%M:%S")

    sh = logging.StreamHandler()
    sh.setFormatter(log_formatter)
    sh.setLevel(logging_level)
    logger.addHandler(sh)

    fh = logging.FileHandler(settings['root']+".log")
    fh.setFormatter(log_formatter)
    logger.addHandler(fh)


    logger.info('CRISPRatac ' + __version__)

    logger.info('Parsing settings file')

    settings['keep_intermediate'] = cmd_args.keep_intermediate
    if 'keep_intermediate' in settings_file_args:
        settings['keep_intermediate'] = (settings_file_args['keep_intermediate'].lower() == 'true')
        settings_file_args.pop('keep_intermediate')


    settings['alignment_extension'] = cmd_args.alignment_extension
    if 'alignment_extension' in settings_file_args:
        settings['alignment_extension'] = int(settings_file_args['alignment_extension'])
        settings_file_args.pop('alignment_extension')

    settings['cuts'] = cmd_args.cuts
    if 'cuts' in settings_file_args:
        settings['cuts'] = settings_file_args['cuts'].split(" ")
        settings_file_args.pop('cuts')
    if 'cut_sites' in settings_file_args:
        settings['cuts'].extend(settings_file_args['cut_sites'].split(" "))
        settings_file_args.pop('cut_sites')
    for cut in settings['cuts']:
        if ":" not in cut:
            parser.print_usage()
            raise Exception('Error: cut specification %s is in the incorrect format (must be given as chr1:234',cut)

    settings['arm_min_seen_bases'] = cmd_args.arm_min_seen_bases
    if 'arm_min_seen_bases' in settings_file_args:
        settings['arm_min_seen_bases'] = int(settings_file_args['arm_min_seen_bases'])
        settings_file_args.pop('arm_min_seen_bases')

    settings['arm_min_matched_start_bases'] = cmd_args.arm_min_matched_start_bases
    if 'arm_min_matched_start_bases' in settings_file_args:
        settings['arm_min_matched_start_bases'] = int(settings_file_args['arm_min_matched_start_bases'])
        settings_file_args.pop('arm_min_matched_start_bases')

    settings['ignore_n'] = cmd_args.ignore_n
    if 'ignore_n' in settings_file_args:
        settings['ignore_n'] = (settings_file_args['ignore_n'].lower() == 'true')
        settings_file_args.pop('ignore_n')

    settings['crispresso_min_count'] = cmd_args.crispresso_min_count
    if 'crispresso_min_count' in settings_file_args:
        settings['crispresso_min_count'] = int(settings_file_args['crispresso_min_count'])
        settings_file_args.pop('crispresso_min_count')

    settings['crispresso_min_aln_score'] = cmd_args.crispresso_min_aln_score
    if 'crispresso_min_aln_score' in settings_file_args:
        settings['crispresso_min_aln_score'] = int(settings_file_args['crispresso_min_aln_score'])
        settings_file_args.pop('crispresso_min_aln_score')

    settings['crispresso_quant_window_size'] = cmd_args.crispresso_quant_window_size
    if 'crispresso_quant_window_size' in settings_file_args:
        settings['crispresso_quant_window_size'] = int(settings_file_args['crispresso_quant_window_size'])
        settings_file_args.pop('crispresso_quant_window_size')


    settings['samtools_command'] = cmd_args.samtools_command
    if 'samtools_command' in settings_file_args:
        settings['samtools_command'] = settings_file_args['samtools_command']
        settings_file_args.pop('samtools_command')

    settings['bowtie2_command'] = cmd_args.bowtie2_command
    if 'bowtie2_command' in settings_file_args:
        settings['bowtie2_command'] = settings_file_args['bowtie2_command']
        settings_file_args.pop('bowtie2_command')

    settings['crispresso_command'] = cmd_args.crispresso_command
    if 'crispresso_command' in settings_file_args:
        settings['crispresso_command'] = settings_file_args['crispresso_command']
        settings_file_args.pop('crispresso_command')

    settings['n_processes'] = cmd_args.n_processes
    if 'n_processes' in settings_file_args:
        settings['n_processes'] = settings_file_args['n_processes']
        settings_file_args.pop('n_processes')

    if settings['n_processes'] == 'max':
        settings['n_processes'] = mp.cpu_count()
    else:
        settings['n_processes'] = int(settings['n_processes'])

    settings['fastq_r1'] = cmd_args.fastq_r1
    if 'fastq_r1' in settings_file_args:
        settings['fastq_r1'] = settings_file_args['fastq_r1']
        settings_file_args.pop('fastq_r1')
    if not settings['fastq_r1']:
        parser.print_usage()
        raise Exception('Error: fastq_r1 file must be provided (--fastq_r1)')
    if not os.path.isfile(settings['fastq_r1']):
        parser.print_usage()
        raise Exception('Error: fastq_r1 file %s does not exist',settings['fastq_r1'])

    settings['fastq_r2'] = cmd_args.fastq_r2
    if 'fastq_r2' in settings_file_args:
        settings['fastq_r2'] = settings_file_args['fastq_r2']
        settings_file_args.pop('fastq_r2')
    if settings['fastq_r2'] and not os.path.isfile(settings['fastq_r2']):
        raise Exception('Error: fastq_r2 file %s does not exist',settings['fastq_r2'])

    settings['genome'] = cmd_args.genome
    if 'genome' in settings_file_args:
        settings['genome'] = settings_file_args['genome']
        settings_file_args.pop('genome')

    if not settings['genome']:
        parser.print_usage()
        raise Exception('Error: the genome reference file must be provided (--genome)')
    if not os.path.isfile(settings['genome']):
        if os.path.isfile(settings['genome']+".fa"):
            settings['genome'] = settings['genome']+".fa"
        else:
            parser.print_usage()
            raise Exception('Error: The genome file %s does not exist'%settings['genome'])
    genome_len_file=settings['genome']+'.fai'
    if not os.path.isfile(genome_len_file):
        raise Exception('Error: The genome length file %s does not exist'%genome_len_file)

    settings['bowtie2_genome'] = cmd_args.bowtie2_genome
    if 'bowtie2_genome' in settings_file_args:
        settings['bowtie2_genome'] = settings_file_args['bowtie2_genome']
        settings_file_args.pop('bowtie2_genome')

    if settings['bowtie2_genome'] is None:
        potential_bowtie2_path = re.sub('.fa$','',settings['genome'])
        if os.path.isfile(potential_bowtie2_path+'.1.bt2'):
            settings['bowtie2_genome']= potential_bowtie2_path
        else:
            raise Exception('Error: bowtie2_genome is required in settings file, pointing to a bowtie2 genome (minus trailing .X.bt2)\nAlternatively, set genome to the .fa file in a bowtie2 directory.')

    unused_keys = settings_file_args.keys()
    if len(unused_keys) > 0:
        unused_key_str = '"'+'", "'.join(unused_keys)+'"'
        raise Exception('Unused keys in settings file: ' + unused_key_str)

    # read previous settings (we'll compare to these later to see if we can cache command outputs)
    settings_used_output_file = settings['root']+".settingsUsed.txt"
    previous_settings = {}
    if os.path.exists(settings_used_output_file):
        with open(settings_used_output_file,'r') as fin:
            for line in fin.readlines():
                settings_els = line.strip().split("\t")
                previous_settings[settings_els[0]] = settings_els[1]


    can_use_previous_analysis = False
    if len(previous_settings) > 0:
        can_use_previous_analysis = True
        for setting in settings:
            if setting in ['debug','n_processes']:
                continue
            if setting not in previous_settings:
                can_use_previous_analysis = False
                logger.info(('Not using previous analyses - got new setting %s (%s)')%(setting,settings[setting]))
                break
            elif str(settings[setting]) != str(previous_settings[setting]):
                can_use_previous_analysis = False
                logger.info(('Not using previous analyses - setting for %s has changed (%s (new) vs %s (old))')%(setting,settings[setting],previous_settings[setting]))
                break

    if can_use_previous_analysis:
        logger.info('Repeated settings detected. Using previous analyses if completed.')

    with open (settings_used_output_file,'w') as fout:
        for setting in settings:
            fout.write("%s\t%s\n"%(str(setting),str(settings[setting])))

    settings['can_use_previous_analysis'] = can_use_previous_analysis

    return settings, logger


def assert_dependencies(samtools_command='samtools',bowtie2_command='bowtie2',crispresso_command='CRISPResso'):
    """
    Asserts the presence of required software (faidx, bowtie2, CRISPResso)

    params:
        samtools_command: location of samtools to run
        bowtie2_command: location of bowtie2 to run
        crispresso_command: location of crispresso to run

    Raises exception if any command is not found
    """
    logger = logging.getLogger('CRISPRatac')
    logger.info('Checking dependencies')

    # check faidx
    try:
        faidx_result = subprocess.check_output('%s faidx'%samtools_command, stderr=subprocess.STDOUT,shell=True)
    except Exception:
        raise Exception('Error: samtools faidx is required')
    if not 'Usage: samtools faidx' in str(faidx_result):
        raise Exception('Error: samtools faidx is required')

    #check bowtie2
    try:
        bowtie_result = subprocess.check_output('%s --version'%bowtie2_command, stderr=subprocess.STDOUT,shell=True)
    except Exception:
        raise Exception('Error: bowtie2 is required')

    #check crispresso
    try:
        crispresso_result = subprocess.check_output('%s --version'%crispresso_command, stderr=subprocess.STDOUT,shell=True)
    except Exception as e:
        raise Exception('Error: CRISPResso2 is required ' + str(e))

nt_complement=dict({'A':'T','C':'G','G':'C','T':'A','N':'N','a':'t','c':'g','g':'c','t':'a','n':'n','_':'_','-':'-'})
def reverse(seq):
        return "".join([c for c in seq[-1::-1]])

def reverse_complement(seq):
        return "".join([nt_complement[c] for c in seq[-1::-1]])

def complement(seq):
    return "".join([nt_complement[c] for c in seq[:]])

def get_av_read_len(fastq,number_reads_to_check=50):
    """
    Reads the first few reads of a file to determine read length

    param:
        fastq: read1 file
        number_reads_to_check: the number of reads to read in

    returns:
        av_read_len: average read length
    """
    sum_len = 0

    cmd=('z' if fastq.endswith('.gz') else '' ) +('cat < \"%s\"' % fastq)+\
               r''' | awk 'BN {n=0;s=0;} NR%4 == 2 {s+=length($0);n++;} END { printf("%d\n",s/n)}' '''
    av_read_len = float(subprocess.check_output(cmd,shell=True).decode(sys.stdout.encoding).strip())

    return(int(av_read_len))

def get_num_reads_fastq(fastq):
    """
    Counts the number of reads in the specified fastq file

    param:
        fastq: fastq file

    returns:
        num_reads: number of reads in the fastq file
    """

    if fastq.endswith('.gz'):
        cmd = 'gunzip -c %s | wc -l'%fastq
    else:
        cmd = 'wc -l %s'%fastq

    res = int(subprocess.check_output(cmd,shell=True).decode(sys.stdout.encoding).split(" ")[0])

    return int(res/4)

def make_artificial_targets(root,cuts,cut_annotations,genome,target_length,target_padding,samtools_command='samtools',bowtie2_command='bowtie2',bowtie2_threads=1,can_use_previous_analysis=False):
    """
    Generates fasta sequences surrounding cuts for alignment and bowtie2 index
    At each cut point, sequence of length target_length is generated
    Combinations of sequences at each cut site are produced

    params:
        cuts: array of cut locations
        cut_annotations: dict of cut_site->annotation for description of cut (e.g. either On-target, Off-target, Known, Casoffinder, etc)
        genome: location of fasta genome
        target_length: how long the query fragment should be (on one side of the cut) (not including padding). 
        target_padding: sequence (bp) padding around target
        samtools_command: location of samtools to run
        bowtie2_command: location of bowtie2 to run
        bowtie2_threads: number of threads to run bowtie2 with
        can_use_previous_analysis: boolean for whether we can use previous analysis or whether the params have changed and we have to rerun from scratch

    returns:
        custom_index_fasta: fasta of artificial targets
        target_names: array of target names (corresponding to targets)
        target_info: hash of information for each target_name
            target_info[target_name]['sequence']: fasta sequence of artificial target
            target_info[target_name]['class']: class of targets (corresponding to targets)
            target_info[target_name]['cut1_chr']: cut information for cut 1
            target_info[target_name]['cut1_site']
            target_info[target_name]['cut1_anno']
            target_info[target_name]['cut2_chr']: cut information for cut 2
            target_info[target_name]['cut2_site']
            target_info[target_name]['cut2_anno']
            target_info[target_name]['query_pos']: genomic start of query (bp)
            target_info[target_name]['target_cut_idx']: bp of cut in target
            target_info[target_name]['target_cut_str']: string designating the cut site, as well as the direction each side of the read comes off of the cut site e.g. w-chr1:50_chr2:60+c means that the left part of the read started on the left side (designated by the -) watson-strand of chr1:50, and the right part of the read started at chr2:60 and extended right on the crick-strand (complement of reference)
    """
    logger = logging.getLogger('CRISPRatac')
    logger.info('Making artificial targets')
    target_names = []
    target_info = {}
    info_file = root + '.info'
    if os.path.isfile(info_file) and can_use_previous_analysis:
        target_count = -1
        with open (info_file,'r') as fin:
            head_line = fin.readline()
            line_els = fin.readline().rstrip('\n').split("\t")
            if len(line_els) == 3:
                (custom_index_fasta,target_count_str,target_list_file) = line_els
                target_count = int(target_count_str)
                read_target_count = 0
                if os.path.isfile(target_list_file):
                    target_cut_str_index = 0
                    with open(target_list_file,'r') as fin:
                        head_line = fin.readline().rstrip('\n')
                        head_line_els = head_line.split("\t")
                        for line in fin:
                            read_target_count += 1
                            (target_name,target_cut_str,target_class,cut1_chr,cut1_site_str,cut1_anno,cut2_chr,cut2_site_str,cut2_anno,query_pos,target_cut_idx,sequence) = line.rstrip('\n').split("\t")

                            cut1_site = int(cut1_site_str)
                            cut2_site = int(cut2_site_str)

                            target_names.append(target_name)
                            target_info[target_name] = {
                                'sequence':sequence,
                                'class':target_class,
                                'cut1_chr':cut1_chr,
                                'cut1_site':int(cut1_site),
                                'cut1_anno':cut1_anno,
                                'cut2_chr':cut2_chr,
                                'cut2_site':int(cut2_site),
                                'cut2_anno':cut2_anno,
                                'query_pos':int(query_pos),
                                'target_cut_idx':int(target_cut_idx),
                                'target_cut_str':target_cut_str
                                }
                    if read_target_count == target_count:
                        logger.info('Using ' + str(target_count) + ' previously-created targets')
                        return custom_index_fasta,target_names,target_info
                    else:
                        logger.info('In attempting to recover previuosly-created custom targets, expecting ' + str(target_count) + ' targets, but only read ' + str(read_target_count) + ' from ' + target_list_file + '. Recreating.')
                else:
                    logger.info('Could not recover previously-created targets. Recreating.')

    fasta_cache = {}

    for i,cut in enumerate(cuts):
        chr_els = cut.split(":")
        chr_A = chr_els[0]
        site_A = int(chr_els[1])
        cut_start_A = site_A - (target_length + target_padding)
        cut_start_A_stop = site_A - 1
        cut_end_A = site_A + (target_length + target_padding)

        left_bit_A_key = '%s %s %d %d'%(genome,chr_A,cut_start_A,cut_start_A_stop)

        if left_bit_A_key not in fasta_cache:
            fasta_cache[left_bit_A_key] = subprocess.check_output(
                '%s faidx -n 10000 %s %s:%d-%d | tail -n 1'%(samtools_command,genome,chr_A,cut_start_A,cut_start_A_stop),shell=True).decode(sys.stdout.encoding).strip()
        left_bit_A = fasta_cache[left_bit_A_key]

        right_bit_A_key = '%s %s %d %d'%(genome,chr_A,site_A,cut_end_A)
        if right_bit_A_key not in fasta_cache:
            fasta_cache[right_bit_A_key] = subprocess.check_output(
                '%s faidx -n 10000 %s %s:%d-%d | tail -n 1'%(samtools_command,genome,chr_A,site_A,cut_end_A),shell=True).decode(sys.stdout.encoding).strip()
        right_bit_A = fasta_cache[right_bit_A_key]

        #only add both sides if primer_chr is not given -- otherwise, every read has to originate on one end from the primer site.
        wt_seq = left_bit_A + right_bit_A
        target_name = 'CRISPRatac_WT'+str(i)
        target_names.append(target_name)
        target_info[target_name] = {
                'sequence': wt_seq,
                'class': 'Linear',
                'cut1_chr':chr_A,
                'cut1_site':site_A,
                'cut1_anno':cut_annotations[cut],
                'cut2_chr':chr_A,
                'cut2_site':site_A,
                'cut2_anno':cut_annotations[cut],
                'query_pos':cut_start_A,
                'target_cut_idx':target_padding + target_length,
                'target_cut_str':"w-%s:%s~%s:%s+w"%(chr_A,site_A,chr_A,site_A)
                }

        LALA = left_bit_A + reverse(left_bit_A)
        target_name = 'CRISPRatac_L' + str(i) + 'L' + str(i)
        target_names.append(target_name)
        target_info[target_name] = {
                'sequence': LALA,
                'class': 'Chimera',
                'cut1_chr':chr_A,
                'cut1_site':site_A,
                'cut1_anno':cut_annotations[cut],
                'cut2_chr':chr_A,
                'cut2_site':site_A,
                'cut2_anno':cut_annotations[cut],
                'query_pos':cut_start_A,
                'target_cut_idx':target_padding + target_length,
                'target_cut_str':"w-%s:%s~%s:%s-w"%(chr_A,site_A,chr_A,site_A)
                }


        LALAc = left_bit_A + reverse_complement(left_bit_A)
        target_name = 'CRISPRatac_L' + str(i) + 'Lc' + str(i)
        target_names.append(target_name)
        target_info[target_name] = {
                'sequence': LALAc,
                'class': 'Chimera',
                'cut1_chr':chr_A,
                'cut1_site':site_A,
                'cut1_anno':cut_annotations[cut],
                'cut2_chr':chr_A,
                'cut2_site':site_A,
                'cut2_anno':cut_annotations[cut],
                'query_pos':cut_start_A,
                'target_cut_idx':target_padding + target_length,
                'target_cut_str':"w-%s:%s~%s:%s-c"%(chr_A,site_A,chr_A,site_A)
                }

        RARA = right_bit_A + reverse(right_bit_A)
        target_name = 'CRISPRatac_R' + str(i) + 'R' + str(i)
        target_names.append(target_name)
        target_info[target_name] = {
                'sequence': RARA,
                'class': 'Chimera',
                'cut1_chr':chr_A,
                'cut1_site':site_A,
                'cut1_anno':cut_annotations[cut],
                'cut2_chr':chr_A,
                'cut2_site':site_A,
                'cut2_anno':cut_annotations[cut],
                'query_pos':cut_start_A,
                'target_cut_idx':target_padding + target_length,
                'target_cut_str':"w+%s:%s~%s:%s+w"%(chr_A,site_A,chr_A,site_A)
                }

        RARAc = right_bit_A + reverse_complement(right_bit_A)
        target_name = 'CRISPRatac_R' + str(i) + 'Rc' + str(i)
        target_names.append(target_name)
        target_info[target_name] = {
                'sequence': RARAc,
                'class': 'Chimera',
                'cut1_chr':chr_A,
                'cut1_site':site_A,
                'cut1_anno':cut_annotations[cut],
                'cut2_chr':chr_A,
                'cut2_site':site_A,
                'cut2_anno':cut_annotations[cut],
                'query_pos':cut_start_A,
                'target_cut_idx':target_padding + target_length,
                'target_cut_str':"w+%s:%s~%s:%s+c"%(chr_A,site_A,chr_A,site_A)
                }


        for j in range(i+1,len(cuts)):
            chr_els = cuts[j].split(":")
            chr_B = chr_els[0]
            site_B = int(chr_els[1])
            cut_start_B = site_B - (target_length + target_padding)
            cut_start_B_stop = site_B - 1
            cut_end_B = site_B + (target_length + target_padding)


            left_bit_B_key = '%s %s %d %d'%(genome,chr_B,cut_start_B,cut_start_B_stop)
            if left_bit_B_key not in fasta_cache:
                fasta_cache[left_bit_B_key] = subprocess.check_output(
                    '%s faidx -n 10000 %s %s:%d-%d | tail -n 1'%(samtools_command,genome,chr_B,cut_start_B,cut_start_B_stop),shell=True).decode(sys.stdout.encoding).strip()
            left_bit_B = fasta_cache[left_bit_B_key]

            right_bit_B_key = '%s %s %d %d'%(genome,chr_B,site_B,cut_end_B)
            if right_bit_B_key not in fasta_cache:
                fasta_cache[right_bit_B_key] = subprocess.check_output(
                    '%s faidx -n 10000 %s %s:%d-%d | tail -n 1'%(samtools_command,genome,chr_B,site_B,cut_end_B),shell=True).decode(sys.stdout.encoding).strip()
            right_bit_B = fasta_cache[right_bit_B_key]


            LARB = left_bit_A + right_bit_B
            target_name = 'CRISPRatac_L' + str(i) + 'R' + str(j)
            target_names.append(target_name)
            target_info[target_name] = {
                    'sequence': LARB,
                    'cut1_chr':chr_A,
                    'cut1_site':site_A,
                    'cut1_anno':cut_annotations[cuts[i]],
                    'cut2_chr':chr_B,
                    'cut2_site':site_B,
                    'cut2_anno':cut_annotations[cuts[j]],
                    'query_pos':cut_start_A,
                    'target_cut_idx':target_padding + target_length,
                    'target_cut_str':"w-%s:%s~%s:%s+w"%(chr_A,site_A,chr_B,site_B)
                    }
            if chr_A == chr_B:
                target_info[target_name]['class'] = 'Large deletion'
            else:
                target_info[target_name]['class'] = 'Translocation'

            LARBc = left_bit_A + complement(right_bit_B)
            target_name = 'CRISPRatac_L' + str(i) + 'Rc' + str(j)
            target_names.append(target_name)
            target_info[target_name] = {
                    'sequence': LARBc,
                    'cut1_chr':chr_A,
                    'cut1_site':site_A,
                    'cut1_anno':cut_annotations[cuts[i]],
                    'cut2_chr':chr_B,
                    'cut2_site':site_B,
                    'cut2_anno':cut_annotations[cuts[j]],
                    'query_pos':cut_start_A,
                    'target_cut_idx':target_padding + target_length,
                    'target_cut_str':"w-%s:%s~%s:%s+c"%(chr_A,site_A,chr_B,site_B)
                    }
            if chr_A == chr_B:
                target_info[target_name]['class'] = 'Large deletion'
            else:
                target_info[target_name]['class'] = 'Translocation'

            LBRA = left_bit_B + right_bit_A
            target_name = 'CRISPRatac_L' + str(j) + 'R' + str(i)
            target_names.append(target_name)
            target_info[target_name] = {
                    'sequence': LBRA,
                    'cut1_chr':chr_B,
                    'cut1_site':site_B,
                    'cut1_anno':cut_annotations[cuts[j]],
                    'cut2_chr':chr_A,
                    'cut2_site':site_A,
                    'cut2_anno':cut_annotations[cuts[i]],
                    'query_pos':cut_start_B,
                    'target_cut_idx':target_padding + target_length,
                    'target_cut_str':"w-%s:%s~%s:%s+w"%(chr_B,site_B,chr_A,site_A)
                    }
            if chr_A == chr_B:
                target_info[target_name]['class'] = 'Large deletion'
            else:
                target_info[target_name]['class'] = 'Translocation'

            LBRAc = left_bit_B + complement(right_bit_A)
            target_name = 'CRISPRatac_L' + str(j) + 'Rc' + str(i)
            target_names.append(target_name)
            target_info[target_name] = {
                    'sequence': LBRAc,
                    'cut1_chr':chr_B,
                    'cut1_site':site_B,
                    'cut1_anno':cut_annotations[cuts[j]],
                    'cut2_chr':chr_A,
                    'cut2_site':site_A,
                    'cut2_anno':cut_annotations[cuts[i]],
                    'query_pos':cut_start_B,
                    'target_cut_idx':target_padding + target_length,
                    'target_cut_str':"w-%s:%s~%s:%s+c"%(chr_B,site_B,chr_A,site_A)
                    }
            if chr_A == chr_B:
                target_info[target_name]['class'] = 'Large deletion'
            else:
                target_info[target_name]['class'] = 'Translocation'

            LALB = left_bit_A + reverse(left_bit_B)
            target_name = 'CRISPRatac_L' + str(i) + 'L' + str(j)
            target_names.append(target_name)
            target_info[target_name] = {
                    'sequence': LALB,
                    'cut1_chr':chr_A,
                    'cut1_site':site_A,
                    'cut1_anno':cut_annotations[cuts[i]],
                    'cut2_chr':chr_B,
                    'cut2_site':site_B,
                    'cut2_anno':cut_annotations[cuts[j]],
                    'query_pos':cut_start_A,
                    'target_cut_idx':target_padding + target_length,
                    'target_cut_str':"w-%s:%s~%s:%s-w"%(chr_A,site_A,chr_B,site_B)
                    }
            if chr_A == chr_B:
                target_info[target_name]['class'] = 'Large inversion'
            else:
                target_info[target_name]['class'] = 'Translocation'

            LALBc = left_bit_A + reverse_complement(left_bit_B)
            target_name = 'CRISPRatac_L' + str(i) + 'Lc' + str(j)
            target_names.append(target_name)
            target_info[target_name] = {
                    'sequence': LALBc,
                    'cut1_chr':chr_A,
                    'cut1_site':site_A,
                    'cut1_anno':cut_annotations[cuts[i]],
                    'cut2_chr':chr_B,
                    'cut2_site':site_B,
                    'cut2_anno':cut_annotations[cuts[j]],
                    'query_pos':cut_start_A,
                    'target_cut_idx':target_padding + target_length,
                    'target_cut_str':"w-%s:%s~%s:%s-c"%(chr_A,site_A,chr_B,site_B)
                    }
            if chr_A == chr_B:
                target_info[target_name]['class'] = 'Large inversion'
            else:
                target_info[target_name]['class'] = 'Translocation'

            RARB = reverse(right_bit_A) + right_bit_B
            target_name = 'CRISPRatac_R' + str(i) + 'R' + str(j)
            target_names.append(target_name)
            target_info[target_name] = {
                    'sequence': RARB,
                    'cut1_chr':chr_A,
                    'cut1_site':site_A,
                    'cut1_anno':cut_annotations[cuts[i]],
                    'cut2_chr':chr_B,
                    'cut2_site':site_B,
                    'cut2_anno':cut_annotations[cuts[j]],
                    'query_pos':cut_start_A,
                    'target_cut_idx':target_padding + target_length,
                    'target_cut_str':"w+%s:%s~%s:%s+w"%(chr_A,site_A,chr_B,site_B)
                    }
            if chr_A == chr_B:
                target_info[target_name]['class'] = 'Large inversion'
            else:
                target_info[target_name]['class'] = 'Translocation'

            RARBc = reverse_complement(right_bit_A) + right_bit_B
            target_name = 'CRISPRatac_Rc' + str(i) + 'R' + str(j)
            target_names.append(target_name)
            target_info[target_name] = {
                    'sequence': RARBc,
                    'cut1_chr':chr_A,
                    'cut1_site':site_A,
                    'cut1_anno':cut_annotations[cuts[i]],
                    'cut2_chr':chr_B,
                    'cut2_site':site_B,
                    'cut2_anno':cut_annotations[cuts[j]],
                    'query_pos':cut_start_A,
                    'target_cut_idx':target_padding + target_length,
                    'target_cut_str':"c+%s:%s~%s:%s+w"%(chr_A,site_A,chr_B,site_B)
                    }
            if chr_A == chr_B:
                target_info[target_name]['class'] = 'Large inversion'
            else:
                target_info[target_name]['class'] = 'Translocation'



    custom_index_fasta = None
    if len(target_names) > 0:
        custom_index_fasta = root + '.fa'
        logger.info('Printing ' + str(len(target_names)) + ' targets to custom index (' + custom_index_fasta + ')')
        with open(custom_index_fasta,'w') as fout:
            for i in range(len(target_names)):
                fout.write('>'+target_names[i]+'\n'+target_info[target_names[i]]['sequence']+'\n')

        logger.info('Indexing custom targets using ' + bowtie2_command + '-build (' + custom_index_fasta + ')')
        this_command = bowtie2_command + '-build --offrate 3 --threads ' + str(bowtie2_threads) + ' ' + custom_index_fasta + ' ' + custom_index_fasta
        logger.debug('Bowtie build command: ' + this_command)
        index_result = subprocess.check_output(this_command, shell=True,stderr=subprocess.STDOUT)

    target_list_file = root+".txt"
    with open (target_list_file,'w') as fout:
        fout.write('\t'.join(['target_name','target_cut_str','target_class','cut1_chr','cut1_site','cut1_anno','cut2_chr','cut2_site','cut2_anno','query_pos','target_cut_idx','sequence'])+"\n")
        for target_name in target_names:
            fout.write(target_name+'\t'+'\t'.join([str(target_info[target_name][x]) for x in ['target_cut_str','class','cut1_chr','cut1_site','cut1_anno','cut2_chr','cut2_site','cut2_anno','query_pos','target_cut_idx','sequence']])+"\n")


    #write info file for restarting
    with open(info_file,'w') as fout:
        fout.write("\t".join([str(x) for x in ["custom_index_fasta","target_count","target_list_file"]])+"\n")
        fout.write("\t".join([str(x) for x in [custom_index_fasta,str(len(target_names)),target_list_file]])+"\n")
    return(custom_index_fasta,target_names,target_info)


def reverse_complement_cut_str(target_cut_str):
    """
    Reverse-complements a cut string of the form 'w-chr1:50~chr2:90+c' to 'w+chr2:90~chr1:50-c'
    Useful for changing cut if a read aligned to the reverse complement of a target

    params:
        target_cut_str: cut string in form 'w-chr1:50~chr2:90+c'
    returns:
        rc_target_cut_str: reverse complement target_cut_str
    """
    (left_whole,right_whole) = target_cut_str.split("~")
    left_strand = left_whole[0]
    left_dir = left_whole[1]
    left_pos = left_whole[2:]

    right_strand = right_whole[-1]
    right_dir = right_whole[-2:-1]
    right_pos = right_whole[0:-2]

    new_left_strand = 'w' if right_strand == 'c' else 'c'
    new_right_strand = 'w' if left_strand == 'c' else 'c'

    return('%s%s%s~%s%s%s'%(new_left_strand,right_dir,right_pos,left_pos,left_dir,new_right_strand))

def align_reads(root,fastq_r1,fastq_r2,reads_name,bowtie2_reference,reference_name,target_info,arm_min_seen_bases=10,arm_min_matched_start_bases=5,bowtie2_command='bowtie2',bowtie2_threads=1,samtools_command='samtools',ignore_n=False,keep_intermediate=False,can_use_previous_analysis=False):
    """
    Aligns reads to the provided reference (either artificial targets or genome)

    params:
        root: root for written files
        fastq_r1: fastq_r1 to align
        fastq_r2: fastq_r2 to align
        reads_name: Name displayed to user about read origin (e.g. 'R1' or 'R2')
        bowtie2_reference: bowtie2 reference to align to (either artificial targets or reference)
        reference_name: Name displayed to user for updates (e.g. 'Genome' or 'Artificial Targets')
        target_info: hash of information for each target_name
        arm_min_seen_bases: number of bases that are required to be seen on each 'side' of artifical targets. E.g. if a artificial target represents a translocation between chr1 and chr2, arm_min_seen_bases would have to be seen on chr1 as well as on chr2 for the read to be counted.
        arm_min_matched_start_bases: number of bases that are required to be matching (no indels or mismatches) at the beginning of the read on each 'side' of artifical targets. E.g. if a artificial target represents a translocation between chr1 and chr2, the first arm_min_matched_start_bases of the read would have to match exactly to chr1 and the last arm_min_matched_start_bases of the read would have to match exactly to chr2
        bowtie2_command: location of bowtie2 to run
        bowtie2_threads: number of threads to run bowtie2 with
        samtools_command: location of samtools to run
        ignore_n: boolean whether to ignore N bases (if False, they count as mismatches)
        keep_intermediate: whether to keep intermediate files (if False, intermediate files will be deleted)
        can_use_previous_analysis: boolean for whether we can use previous analysis or whether the params have changed and we have to rerun from scratch

    returns:
        r1_assignments_file: file showing locations of assignments for r1
        r2_assignments_file: file showing locations of assignments for r2
        unmapped_fastq_r1_file: fastq_r1 file of reads not aligning
        unmapped_fastq_r2_file: fastq_r2 file of reads not aligning
        aligned_count: number of reads aligned
        mapped_bam_file: aligned reads
        chr_aln_plot_obj: plot object showing chr locations of alignments
        tlen_plot_object: plot object showing insert sizes for paired alignments
    """

    logger = logging.getLogger('CRISPRatac')
    logger.info('Aligning %s to %s'%(reads_name,reference_name.lower()))

    mapped_bam_file = root + ".bam"

    info_file = root + '.info'
    if os.path.isfile(info_file) and can_use_previous_analysis:
        read_count = -1
        with open (info_file,'r') as fin:
            head_line = fin.readline()
            line_els = fin.readline().rstrip('\n').split("\t")
            if len(line_els) == 11:
                (read_count_str,r1_count_str,r2_count_str,unmapped_r1_count_str,unmapped_r2_count_str,aligned_count_str,r1_assignments_file_str,r2_assignments_file_str,unmapped_fastq_r1_file_str,unmapped_fastq_r2_file_str,mapped_bam_file_str) = line_els
                r1_assignments_file = None if r1_assignments_file_str == "None" else r1_assignments_file_str
                r2_assignments_file = None if r2_assignments_file_str == "None" else r2_assignments_file_str
                r1_count = int(r1_count_str)
                r2_count = int(r2_count_str)
                unmapped_r1_count = int(unmapped_r1_count_str)
                unmapped_r2_count = int(unmapped_r2_count_str)
                unmapped_fastq_r1_file = None if unmapped_fastq_r1_file_str == "None" else unmapped_fastq_r1_file_str
                unmapped_fastq_r2_file = None if unmapped_fastq_r2_file_str == "None" else unmapped_fastq_r2_file_str
                read_count  = int(read_count_str)
                aligned_count  = int(aligned_count_str)
                mapped_bam_file = None if mapped_bam_file == "None" else mapped_bam_file_str

                chr_aln_plot_obj_str = fin.readline().rstrip('\n')
                chr_aln_plot_obj = None
                if chr_aln_plot_obj_str != "" and chr_aln_plot_obj_str != "None":
                    chr_aln_plot_obj = PlotObject.from_json(chr_aln_plot_obj_str)

                tlen_plot_obj_str = fin.readline().rstrip('\n')
                tlen_plot_obj = None
                if tlen_plot_obj_str != "" and tlen_plot_obj_str != "None":
                    tlen_plot_obj = PlotObject.from_json(tlen_plot_obj_str)
                if read_count > -1:
                    r1_aln_count = r1_count-unmapped_r1_count
                    r2_aln_count = r2_count-unmapped_r2_count
                    if r2_count > -1:
                        logger.info('Using previously-processed alignment of %d/%d R1 and %d/%d R2 %s (%d reads input) aligned to %s'%(r1_aln_count,r1_count,r2_aln_count,r2_count,reads_name,read_count,reference_name.lower()))
                    else:
                        logger.info('Using previously-processed alignment of %d/%d %s (%d reads input) aligned to %s'%(r1_aln_count,r1_count,reads_name,read_count,reference_name.lower()))
                    return r1_assignments_file,r2_assignments_file,unmapped_fastq_r1_file, unmapped_fastq_r2_file, aligned_count, mapped_bam_file,chr_aln_plot_obj,tlen_plot_obj
                else:
                    logger.info('Could not recover previously-analyzed alignments. Reanalyzing.')


    bowtie_log = root + '.bowtie2Log'
    if fastq_r2 is not None: #paired-end reads
        logger.info('Aligning paired reads to %s using %s'%(reference_name.lower(),bowtie2_command))
        aln_command = '{bowtie2_command} --seed 2248 --sam-no-qname-trunc --end-to-end --threads {bowtie2_threads} -x {bowtie2_reference} -1 {fastq_r1} -2 {fastq_r2} | {samtools_command} view -F 256 -Shu - | {samtools_command} sort -o {mapped_bam_file} - && {samtools_command} index {mapped_bam_file}'.format(
                bowtie2_command=bowtie2_command,
                bowtie2_threads=bowtie2_threads,
                bowtie2_reference=bowtie2_reference,
                fastq_r1=fastq_r1,
                fastq_r2=fastq_r2,
                samtools_command=samtools_command,
                mapped_bam_file=mapped_bam_file)
        logger.debug(aln_command)
        aln_result = subprocess.check_output(aln_command,shell=True,stderr=subprocess.STDOUT).decode(sys.stdout.encoding)
        if 'error' in aln_result.lower():
            logger.error('Error found while running command:\n'+aln_command+"\nOutput: "+aln_result)
            raise Exception('Alignment error: ' + aln_result)

        logger.debug(aln_result)
        with open (bowtie_log,'w') as lout:
            lout.write('\nAligning paired reads to %s\n===\nCommand used:\n===\n%s\n===\nOutput:\n===\n%s'%(reference_name.lower(),aln_command,aln_result))
    #unpaired reads
    else:
        logger.info('Aligning single-end reads to %s using %s'%(reference_name.lower(),bowtie2_command))
        aln_command = '{bowtie2_command} --seed 2248 --sam-no-qname-trunc --end-to-end --threads {bowtie2_threads} -x {bowtie2_reference} -U {fastq_r1} | {samtools_command} view -F 256 -Shu - | {samtools_command} sort -o {mapped_bam_file} - && {samtools_command} index {mapped_bam_file}'.format(
                bowtie2_command=bowtie2_command,
                bowtie2_threads=bowtie2_threads,
                bowtie2_reference=bowtie2_reference,
                fastq_r1=fastq_r1,
                samtools_command=samtools_command,
                mapped_bam_file=mapped_bam_file)
        logger.debug(aln_command)
        aln_result = subprocess.check_output(aln_command,shell=True,stderr=subprocess.STDOUT).decode(sys.stdout.encoding)
        if 'error' in aln_result.lower():
            logger.error('Error found while running command:\n'+aln_command+"\nOutput: "+aln_result)
            raise Exception('Alignment error: ' + aln_result)

        logger.debug(aln_result)
        with open (bowtie_log,'w') as lout:
            lout.write('\nAligning single-end reads to %s\n===\nCommand used:\n===\n%s\n===\nOutput:\n===\n%s'%(reference_name.lower(),aln_command,aln_result))

    logger.info('Analyzing reads aligned to %s'%(reference_name.lower()))

    #analyze alignments
    # - if read aligned, store read id and loc in dict r1_assignments
    # - if read unaligned, write read to unaligned fastq
    mapped_tlens = defaultdict(int) # observed fragment lengths
    mapped_chrs = {} #stores the counts aligning to chrs or to templates
    aligned_locs = {} # aligned reads will not be chopped, but keep track of where they aligned for CRISPResso output
    aligned_chr_counts = defaultdict(int)

    unmapped_fastq_r1_file = root + '.unmapped.fq'
    unmapped_fastq_r2_file = None
    r1_assignments_file = root + '.assignments.txt'
    r2_assignments_file = None
    if fastq_r2 is not None: #paired-end reads
        r1_assignments_file = root + '.assignments_r1.txt'
        unmapped_fastq_r1_file = root + '.unmapped_r1.fq'
        r2_assignments_file = root + '.assignments_r2.txt'
        unmapped_fastq_r2_file = root + '.unmapped_r2.fq'
        uf2 = open(unmapped_fastq_r2_file,'w')
        af2 = open(r2_assignments_file,'w')
        af2.write("read_id\treference_name\tclassification\tannotation\tposition\tcut\n")
    uf1 = open(unmapped_fastq_r1_file,'w')
    af1 = open(r1_assignments_file,'w')
    af1.write("read_id\treference_name\tclassification\tannotation\tposition\tcut\n")

    read_count = 0
    r1_count = 0
    r2_count = 0
    unmapped_r1_count = 0
    unmapped_r2_count = 0

    #-F 256 - not primary aligment (filter secondary alignments)
    for line in read_command_output('%s view -F 256 %s'%(samtools_command,mapped_bam_file)):
        if line.strip() == "": break

        line_els = line.split("\t")
        read_count += 1

        read_has_multiple_segments = int(line_els[1]) & 0x1
        read_is_paired_read_1 = int(line_els[1]) & 0x40 # only set with bowtie if read is from paired reads

        read_is_r1 = True
        if read_has_multiple_segments and not read_is_paired_read_1:
            r2_count += 1
            read_is_r1 = False
        else:
            r1_count += 1

        read_id = line_els[0]


        is_bad_alignment = False
        line_unmapped = int(line_els[1]) & 0x4

        left_matches = 0
        right_matches = 0
        mez_val = ""
        for line_el in line_els:
            if line_el.startswith('MD:Z'):
                mez_val = line_el[5:]
                left_matches,right_matches = getLeftRightMismatchesMDZ(mez_val)
                break

        seq = line_els[9]
        is_rc = int(line_els[1]) & 0x10

        if line_unmapped or \
                left_matches < arm_min_matched_start_bases or \
                right_matches < arm_min_matched_start_bases:

            qual = line_els[10]

            if is_rc:
                seq = reverse_complement(seq)

            if read_is_r1:
                uf1.write("@%s\n%s\n+\n%s\n"%(read_id,seq,qual))
                unmapped_r1_count += 1
            else:
                uf2.write("@%s\n%s\n+\n%s\n"%(read_id,seq,qual))
                unmapped_r2_count += 1
            continue

        line_chr = line_els[2]
        line_mapq = line_els[5]
        line_start = int(line_els[3])
        line_end = line_start+(len(seq)-1)
        if is_rc:
            tmp = line_end
            line_end = line_start
            line_start = tmp

        curr_classification = 'Linear'
        curr_position = "%s:%s-%s"%(line_chr,line_start,line_end)
        curr_annotation = ''
        curr_cut = 'NA'
        if line_chr in target_info: #if this aligned to a custom chromosome
            (left_read_bases_count, left_ref_bases_count, left_all_match_count, left_start_match_count, right_read_bases_count, right_ref_bases_count, right_all_match_count, right_start_match_count) = getMatchLeftRightOfCut(target_info[line_chr]['sequence'],line,target_info[line_chr]['target_cut_idx'],ignore_n=ignore_n)

            #if read doesn't sufficiently align to both parts of the artificial target, print it to unmapped and continue
            if left_read_bases_count < arm_min_seen_bases or right_read_bases_count < arm_min_seen_bases or left_start_match_count < arm_min_matched_start_bases or right_start_match_count < arm_min_matched_start_bases:

                qual = line_els[10]

                if is_rc:
                    seq = reverse_complement(seq)
                    qual = qual[::-1]

                if read_is_r1:
                    uf1.write("@%s\n%s\n+\n%s\n"%(read_id,seq,qual))
                    unmapped_r1_count += 1
                else:
                    uf2.write("@%s\n%s\n+\n%s\n"%(read_id,seq,qual))
                    unmapped_r2_count += 1
                continue

            curr_classification = target_info[line_chr]['class']
            curr_annotation = target_info[line_chr]['target_cut_str'] + ' ' + line_chr
            if is_rc:
                curr_annotation = reverse_complement_cut_str(target_info[line_chr]['target_cut_str']) + ' ' + line_chr

            cut1 =  target_info[line_chr]['cut1_chr'] + ":" +  str(target_info[line_chr]['cut1_site'])
            cut2 =  target_info[line_chr]['cut2_chr'] + ":" +  str(target_info[line_chr]['cut2_site'])

            left_dir = target_info[line_chr]['target_cut_str'][1:2]
            if left_dir == '-': #arm extends left in genome space after cut
                left_arm_end = target_info[line_chr]['cut1_site']
                left_arm_start = left_arm_end - left_ref_bases_count
            else: #arm extends right in genome space after cut
                left_arm_start = target_info[line_chr]['cut1_site']
                left_arm_end = left_arm_start + left_ref_bases_count

            right_dir = target_info[line_chr]['target_cut_str'][-2:-1]
            if right_dir == '+': #arm extends right in genome space after cut
                right_arm_start = target_info[line_chr]['cut2_site']
                right_arm_end = right_arm_start + right_ref_bases_count
            else:
                right_arm_end = target_info[line_chr]['cut2_site']
                right_arm_start = right_arm_end - right_ref_bases_count


            if curr_classification == 'Linear':
                #if read aligned to custom seqs in the forward direction
                if not is_rc:
                    curr_cut = cut1
                    curr_position = '%s:%s-%s'%(
                            target_info[line_chr]['cut1_chr'],left_arm_start,right_arm_end - 1
                            )
                else:
                    curr_cut = cut1
                    curr_position = '%s:%s-%s'%(
                            target_info[line_chr]['cut1_chr'],right_arm_end - 1,left_arm_start
                            )

            #if custom ref is not linear, set cuts and genome aln locations
            else:
                #if read aligned to custom seqs in the forward direction
                if not is_rc:
                    curr_cut = cut1 + "~" + cut2
                    curr_position = '%s:%s-%s~%s:%s-%s'%(
                            target_info[line_chr]['cut1_chr'],left_arm_start,left_arm_end,
                            target_info[line_chr]['cut2_chr'],right_arm_start,right_arm_end - 1
                            )
                #if read aligned to custom seqs in the reverse direction
                else:
                    curr_cut = cut2 + "~" + cut1
                    curr_position = '%s:%s-%s~%s:%s-%s'%(
                            target_info[line_chr]['cut2_chr'],right_arm_end - 1,right_arm_start,
                            target_info[line_chr]['cut1_chr'],left_arm_end,left_arm_start
                            )

        if read_is_r1:
            af1.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(read_id,reference_name,curr_classification,curr_annotation,curr_position,curr_cut))
        else:
            af2.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(read_id,reference_name,curr_classification,curr_annotation,curr_position,curr_cut))

        if read_is_paired_read_1: #only set if paired
            insert_size = int(line_els[8])
            mapped_tlens[insert_size] += 1

        if not line_chr in mapped_chrs:
            mapped_chrs[line_chr] = 0
        mapped_chrs[line_chr] += 1

        if line_chr not in aligned_locs:
            aligned_locs[line_chr] = {}
        if line_start not in aligned_locs[line_chr]:
            aligned_locs[line_chr][line_start] = 0

        aligned_locs[line_chr][line_start] += 1
        aligned_chr_counts[line_chr] += 1


    uf1.close()
    af1.close()
    if unmapped_fastq_r2_file is not None:
        uf2.close()
        af2.close()

    chr_aln_plot_root = root + ".chrs"
    keys = sorted(aligned_chr_counts.keys())
    vals = [aligned_chr_counts[key] for key in keys]
    with open(chr_aln_plot_root+".txt","w") as chrs:
        chrs.write('chr\tnumReads\n')
        for key in sorted(aligned_chr_counts.keys()):
            chrs.write(key + '\t' + str(aligned_chr_counts[key]) + '\n')

    fig = plt.figure(figsize=(12,12))
    ax = plt.subplot(111)
    if len(vals) > 0:
        ax.bar(range(len(keys)),vals,tick_label=keys)
        ax.set_ymargin(0.05)
    else:
        ax.bar(0,0)
    ax.set_ylabel('Number of Reads')
    ax.set_title('Location of ' + reads_name + ' aligned to the '+reference_name.lower())
    plt.savefig(chr_aln_plot_root+".pdf",pad_inches=1,bbox_inches='tight')
    plt.savefig(chr_aln_plot_root+".png",pad_inches=1,bbox_inches='tight')

    plot_label = 'Bar plot showing alignment location of ' + reads_name + ' aligned to ' + reference_name.lower()
    if unmapped_fastq_r2_file is not None:
        plot_label += ' (Note that this plot shows the number of reads not readpairs)'
    if len(keys) == 0:
        plot_label = '(No ' + reads_name + ' aligned to ' + reference_name.lower() + ')'

    chr_aln_plot_obj = PlotObject(
            plot_name = chr_aln_plot_root,
            plot_title = reference_name.title() + ' alignment summary',
            plot_label = plot_label,
            plot_datas = [(reference_name.capitalize() + ' alignment summary',chr_aln_plot_root + ".txt")]
            )

    tlen_plot_obj = None
    if fastq_r2 is not None: #paired-end reads
        tlen_plot_root = root + ".insertSizes"
        keys = sorted(mapped_tlens.keys())
        vals = [mapped_tlens[key] for key in keys]
        with open(tlen_plot_root+".txt","w") as fout:
            fout.write('insertSize\tnumReads\n')
            for key in keys:
                fout.write(str(key) + '\t' + str(mapped_tlens[key]) + '\n')

        fig = plt.figure(figsize=(12,12))
        ax = plt.subplot(111)
        if len(vals) > 0:
            ax.bar(keys,vals)
            ax.set_ymargin(0.05)
        else:
            ax.bar(0,0)
        ax.set_ylabel('Number of Reads')
        ax.set_title('Insert size')
        plt.savefig(tlen_plot_root+".pdf",pad_inches=1,bbox_inches='tight')
        plt.savefig(tlen_plot_root+".png",pad_inches=1,bbox_inches='tight')

        tlen_plot_obj = PlotObject(
                plot_name = tlen_plot_root,
                plot_title = reference_name.capitalize() + ' Alignment Insert Size Summary',
                plot_label = 'Bar plot showing insert size of reads aligned to ' + reference_name,
                plot_datas = [(reference_name.capitalize() + ' alignment insert size summary',tlen_plot_root + ".txt")]
                )

    aligned_count = read_count - (unmapped_r1_count + unmapped_r2_count)

    with open(info_file,'w') as fout:
        fout.write("\t".join([str(x) for x in ["read_count","r1_count","r2_count","unmapped_r1_count","unmapped_r2_count","aligned_count","r1_assignments_file","r2_assignments_file","unmapped_fastq_r1_file","unmapped_fastq_r2_file","mapped_bam_file"]])+"\n")
        fout.write("\t".join([str(x) for x in [read_count,r1_count,r2_count,unmapped_r1_count,unmapped_r2_count,aligned_count,r1_assignments_file,r2_assignments_file,unmapped_fastq_r1_file,unmapped_fastq_r2_file,mapped_bam_file]])+"\n")

        chr_aln_plot_obj_str = "None"
        if chr_aln_plot_obj is not None:
            chr_aln_plot_obj_str = chr_aln_plot_obj.to_json()
        fout.write(chr_aln_plot_obj_str+"\n")
        tlen_plot_obj_str = "None"
        if tlen_plot_obj is not None:
            tlen_plot_obj_str = tlen_plot_obj.to_json()
        fout.write(tlen_plot_obj_str+"\n")

    r1_aln_count = r1_count-unmapped_r1_count
    r2_aln_count = r2_count-unmapped_r2_count
    if r2_aln_count > 0:
        logger.info('Finished analysis of %d/%d R1 and %d/%d R2 %s (%d reads input) aligned to %s'%(r1_aln_count,r1_count,r2_aln_count,r2_count,reads_name,read_count,reference_name.lower()))
    else:
        logger.info('Finished analysis of %d/%d %s (%d reads input) aligned to %s'%(r1_aln_count,r1_count,reads_name,read_count,reference_name.lower()))
    return(r1_assignments_file,r2_assignments_file,unmapped_fastq_r1_file,unmapped_fastq_r2_file,aligned_count,mapped_bam_file,chr_aln_plot_obj,tlen_plot_obj)

def run_crispresso2(root,input_fastq_file,assignment_file,target_info,crispresso_name_suffix,crispresso_min_count,crispresso_min_aln_score,crispresso_quant_window_size,crispresso_command='CRISPResso',n_processes=1):
    """
    Runs CRISPResso to highlight changes in each custom amplicon

    Params:
        root: root for written files
        input_fastq_file
        assignment_file: list of assignment files with reads to analyze with CRISPResso
        target_info: hash of information for each target_name
        crispresso_name_suffix: string to append to CRISPResso runs (e.g. 'r1' or 'r2')
        crispresso_min_count: min number of reads at site for crispresso2 processing
        crispresso_min_aln_score: minimum score for reads to align to amplicons
        crispresso_quant_window_size: number of bases on each side of the cut site in which to consider editing events
        crispresso_command: location of crispresso to run
        n_processes: number of processes to run CRISPResso commands on

    Returns:
        crispresso_infos: dict containing info for each crispresso run (e.g. command used, etc.)
        crispresso_run_names: list of crispresso run names
        crispresso_sub_htmls: dict of name-> location of html file for crispresso run
    """

    logger = logging.getLogger('CRISPRatac')
    logger.info('Analysis of ' + crispresso_name_suffix + ' custom targets by CRISPResso2')

    id_ind = 0
    source_ind = 1
    classification_ind = 2
    annotation_ind = 3
    alignment_ind = 4
    cut_point_ind = 5

    aligned_target_counts = defaultdict(int)
    read_ids_for_crispresso = {}

    with open(assignment_file, 'r') as fin:
        head = fin.readline()
        for line in fin:
            line_els = line.split("\t")
            target_els = line_els[annotation_ind].split(" ")
            if len(target_els) < 2:
                print("Can't parse target from " + str(target_els))
                target_name = target_els[0]
            else:
                target_name = target_els[1]
            read_ids_for_crispresso[line_els[id_ind]] = target_name
            aligned_target_counts[target_name] += 1
            
    crispresso_infos = [] #information about each crispresso_name

    data_dir = root+'_data'
    if not os.path.isdir(data_dir):
        os.mkdir(data_dir)

    total_read_count = 0
    printed_read_count = 0
    filehandles = {}

    printed_targets = defaultdict(int)
    target_names = {}

    #iterate through fastq file
    # if read id is in the read_ids_for_crispresso, print it to that file
    if input_fastq_file.endswith('.gz'):
        f_in = gzip.open(input_fastq_file,'rt')
    else:
        f_in = open(input_fastq_file,'r')
    while (1):
        id_line   = f_in.readline().strip()
        seq_line  = f_in.readline().strip()
        plus_line = f_in.readline()
        qual_line = f_in.readline().strip()

        if not qual_line : break

        total_read_count+=1

        id_line_els = id_line.split("\t")
        read_id = id_line_els[0][1:] #chop off leading @

        if read_id not in read_ids_for_crispresso:
            continue

        target = read_ids_for_crispresso[read_id]
        if target not in target_names:
            target_names[target] = target +'_'+ crispresso_name_suffix
        target_name = target_names[target]

        if aligned_target_counts[target] < crispresso_min_count:
            continue

        printed_read_count += 1
        reads_file = os.path.join(data_dir,target_name+".fq")
        if reads_file not in filehandles:
            fh = open(reads_file,'w')
            filehandles[reads_file] = fh
        filehandles[reads_file].write("%s\n%s\n%s%s\n"%(id_line,seq_line.upper(),plus_line,qual_line))

        printed_targets[target] += 1

    crispresso_commands = []
    for i,target in enumerate(printed_targets):
        target_name = target_names[target]
        amp_seq = target_info[target]['sequence']
        amp_len = len(amp_seq)
        guide_len = 10
        half_amp_len = int(amp_len/2)
        guide_seq = amp_seq[half_amp_len-guide_len:half_amp_len+guide_len]
        reads_file = os.path.join(data_dir,target_name+".fq")
        processes_str = ""
        if n_processes > 4:
            processes_str = "--n_processes 4"
        output_folder = os.path.join(data_dir,'CRISPResso_on_'+target_name)
        crispresso_cmd = "%s -o %s -n %s --default_min_aln_score %d -a %s -g %s -gn cut_site -wc %s -w %s -r1 %s --fastq_output %s &> %s.log"%(crispresso_command,data_dir,target_name,crispresso_min_aln_score,
            amp_seq,guide_seq,-guide_len,crispresso_quant_window_size,reads_file,processes_str,reads_file)

        crispresso_commands.append(crispresso_cmd)
        crispresso_infos.append({
                "name":target_name,
                "amp_seq": amp_seq,
                "output_folder":output_folder,
                "reads_file": reads_file,
                "printed_read_count": printed_targets[target],
                "command": crispresso_cmd
                })

    CRISPRessoMultiProcessing.run_crispresso_cmds(crispresso_commands,n_processes,'run')

    crispresso_run_names = []
    crispresso_sub_htmls = {}
    crispresso_info_file = root + '.summary.txt'
    with open(crispresso_info_file,'w') as crispresso_info_fh:
        crispresso_info_fh.write('name\tprinted_read_count\tcommand\tsuccess\n')
        for crispresso_info in crispresso_infos:
            name = crispresso_info['name']
            run_success = False
            try:
                run_data = CRISPRessoShared.load_crispresso_info(crispresso_info['output_folder'])
            except:
                logger.debug('Could not load CRISPResso run information from ' + crispresso_info['name'] + ' at ' + crispresso_info['output_folder'])
            else:
                if 'running_info' in run_data:
                    run_success = True
                    report_filename = run_data['running_info']['report_filename']
                    report_file_loc = os.path.join(os.path.dirname(crispresso_info['output_folder']),report_filename)
                    if os.path.isfile(report_file_loc):
                        crispresso_run_names.append(name)
                        crispresso_sub_htmls[name] = report_file_loc
            crispresso_info_fh.write("\t".join([str(x) for x in [name,crispresso_info['printed_read_count'],crispresso_info['command'],run_success]])+"\n")

    return (crispresso_run_names,crispresso_sub_htmls)



def read_command_output(command):
    """
    Runs a shell command and returns an iter to read the output

    param:
        command: shell command to run

    returns:
        iter to read the output
    """

    p = subprocess.Popen(command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,shell=True,
#            encoding='utf-8',universal_newlines=True)
            universal_newlines=True,
            bufsize=-1) #bufsize system default
    return iter(p.stdout.readline, b'')

def getLeftRightMismatchesMDZ(mdz_str):
    """
    Gets the number of mismates on the left and right of an alignment from the MD:Z string

    From the Samtool spec:

    The MD string consists of the following items, concatenated without additional delimiter characters:
      [0-9]+, indicating a run of reference bases that are identical to the corresponding SEQ bases;
      [A-Z], identifying a single reference base that differs from the SEQ base aligned at that position;
      \^[A-Z]+, identifying a run of reference bases that have been deleted in the alignment.

    params: 
        mdz_str: MD:Z string from alignment

    returns:
        left_matches: number of bp that match on the left side of the alignment
        right_matches
   """
    num_set = set([str(x) for x in range(10)])
    left_matches = 0
    curr_num_str = ""
    for str_index in range(len(mdz_str)):
        if mdz_str[str_index] in num_set:
            curr_num_str += mdz_str[str_index]
        else:
            break

    if curr_num_str == '':
        left_matches = 0
    else:
        left_matches = int(curr_num_str)

    right_matches = 0
    curr_num_str = ""
    for str_index in range(len(mdz_str)-1,-1,-1):
        if mdz_str[str_index] in num_set:
            curr_num_str = mdz_str[str_index] + curr_num_str
        else:
            break

    if curr_num_str == '':
        right_matches = 0
    else:
        right_matches = int(curr_num_str)

    return left_matches,right_matches


def sam2aln(ref_seq,sam_line,include_ref_surrounding_read = True):
    """
    Creates an alignment of nucleotides from a sam line and a reference sequence by parsing the cigar string

    params:
        ref_seq: string, reference sequence
        sam_line: tab-sep string, sam line. Expected entries are alignment position (element 4) cigar (element 6), and sequence (element 10)
        include_ref_surrounding_read: boolean, whether to include ref sequence to the left or right of the sequenced read. If false, only the reference directly overlapping the read is included

    returns:
        ref_str: reference string with gaps added as appropriate
        aln_str: reference string with clipped bases removed and gap added as appropriate
        clipped_left_bp: number of bp clipped from left side of read - includes both hard and soft-clipped bases
        clipped_right_bp: number of bp clipped from right side of read

    tests from samtools spec:
    (ref_str,aln_str,clipped_left_bp,clipped_right_bp) = sam2aln('AGCATGTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCAT','r001	99	ref	7	30	8M2I4M1D3M	=	37	39	TTAGATAAAGGATACTG	*',False)
    print('ref: ' + ref_str + '\naln: ' + aln_str + '\nclipped: left: ' + str(clipped_left_bp) + ' right: ' + str(clipped_right_bp))

    (ref_str,aln_str,clipped_left_bp,clipped_right_bp) = sam2aln('AGCATGTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCAT','r002	0	ref	9	30	3S6M1P1I4M	*	0	0	AAAAGATAAGGATA	*')
    print('ref: ' + ref_str + '\naln: ' + aln_str + '\nclipped: left: ' + str(clipped_left_bp) + ' right: ' + str(clipped_right_bp))

    (ref_str,aln_str,clipped_left_bp,clipped_right_bp) = sam2aln('AGCATGTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCAT','r002	0	ref	9	30	3S6M1P1I4M5S	*	0	0	AAAAGATAAGGATAGGGGG	*')
    print('ref: ' + ref_str + '\naln: ' + aln_str + '\nclipped: left: ' + str(clipped_left_bp) + ' right: ' + str(clipped_right_bp))

    (ref_str,aln_str,clipped_left_bp,clipped_right_bp) = sam2aln('AGCATGTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCAT','r003	0	ref	9	30	5S6M	*	0	0	GCCTAAGCTAA	*	SA:Z:ref,29,-,6H5M,17,0;')
    print('ref: ' + ref_str + '\naln: ' + aln_str + '\nclipped: left: ' + str(clipped_left_bp) + ' right: ' + str(clipped_right_bp))

    (ref_str,aln_str,clipped_left_bp,clipped_right_bp) = sam2aln('AGCATGTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCAT','r004	0	ref	16	30	6M14N5M	*	0	0	ATAGCTTCAGC	*')
    print('ref: ' + ref_str + '\naln: ' + aln_str + '\nclipped: left: ' + str(clipped_left_bp) + ' right: ' + str(clipped_right_bp))

    (ref_str,aln_str,clipped_left_bp,clipped_right_bp) = sam2aln('AGCATGTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCAT','r003	2064	ref	29	17	6H5M	*	0	0	TAGGC	*	SA:Z:ref,9,+,5S6M,30,1;')
    print('ref: ' + ref_str + '\naln: ' + aln_str + '\nclipped: left: ' + str(clipped_left_bp) + ' right: ' + str(clipped_right_bp))

    (ref_str,aln_str,clipped_left_bp,clipped_right_bp) = sam2aln('AGCATGTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCAT','r001	147	ref	37	30	9M	=	7	-39	CAGCGGCAT	*	NM:i:1')
    print('ref: ' + ref_str + '\naln: ' + aln_str + '\nclipped: left: ' + str(clipped_left_bp) + ' right: ' + str(clipped_right_bp))

    """
    sam_line_els = sam_line.split("\t")
    cigar = sam_line_els[5]

    cigar_pattern = re.compile("(\d+)(\w)")
    cigar_els = []
    for (c_count,c_char) in re.findall(cigar_pattern,cigar):
        cigar_els.append((int(c_count),c_char))

    remaining_ref = ref_seq
    remaining_aln = sam_line_els[9]

    aln_str = ""
    ref_str = ""
    curr_start = int(sam_line_els[3])-1
    if include_ref_surrounding_read:
        aln_str = "-"*curr_start
        ref_str = ref_seq[0:curr_start]

    remaining_ref = ref_seq[curr_start:]

#    print('aln_str: ' + aln_str)
#    print('ref_str: ' + ref_str)
#    print('remaining_aln: ' + remaining_aln)
#    print('remaining_ref: ' + remaining_ref)

    clipped_left_bp = 0 #keep track of hard clipping
    clipped_right_bp = 0

    for idx,(c_count,c_char) in enumerate(cigar_els):
#        print('curr operation: ' + str(c_count) + ' ' + c_char)
        if c_char == 'M':
            aln_str += remaining_aln[0:c_count]
            remaining_aln = remaining_aln[c_count:]
            ref_str += remaining_ref[0:c_count]
            remaining_ref = remaining_ref[c_count:]
        elif c_char == 'I':
            aln_str += remaining_aln[0:c_count]
            remaining_aln = remaining_aln[c_count:]
            ref_str += '-'*c_count
        elif c_char == 'D' or c_char == 'N':
            aln_str += '-'*c_count
            ref_str += remaining_ref[0:c_count]
            remaining_ref = remaining_ref[c_count:]
        elif c_char == 'S':
            remaining_aln = remaining_aln[c_count:]
            if idx == 0:
                clipped_left_bp += c_count
            if idx == len(cigar_els)-1:
                clipped_right_bp += c_count
        elif c_char == 'H' or c_char == 'P':
            if idx == 0:
                clipped_left_bp += c_count
            if idx == len(cigar_els)-1:
                clipped_right_bp += c_count


            pass
        else:
            raise Exception('Unable to parse cigar character: ' + c_char)
#        print('aln_str: ' + aln_str)
#        print('ref_str: ' + ref_str)
#        print('remaining_aln: ' + remaining_aln)
#        print('remaining_ref: ' + remaining_ref)


    if include_ref_surrounding_read:
        aln_str += '-'*len(remaining_ref)
        ref_str += remaining_ref

#    print('Final')
#    print('aln_str: ' + aln_str)
#    print('ref_str: ' + ref_str)

    return(ref_str,aln_str,clipped_left_bp,clipped_right_bp)

def getMatchLeftRightOfCut(ref_seq,sam_line,cut_pos,ignore_n=False,debug=False):
    """
    Gets the number of mismatches at the beginning and the end of the read specified in the sam line. Left and right are defined by the cut position

    params:
        ref_seq: string, reference sequence
        sam_line: tab-sep string, sam line. Expected entries are alignment position (element 4) cigar (element 6), and sequence (element 10)
        cut_pos: the position of the cut that defines left and right. Technically, this cut happens after this many bases.
        ignore_n: boolean whether to ignore N bases (if False, they count as mismatches)
        debug: boolean whether to print debug

    returns:
        left_aln_bases_count: number of bases from this read that were aligned to the left part
        left_ref_bases_count: number of bases from the ref that were aligned to the left part
        left_all_match_count: number of total bases on the left part that matched
        left_start_match_count: number of bases at the beginning of the read that matched exactly
        right_aln_bases_count
        right_ref_bases_count
        right_all_match_count
        right_start_match_count
    """

    (ref_str,aln_str,clipped_left_bp,clipped_right_bp) = sam2aln(ref_seq,sam_line)
    ref_str = ref_str.upper()
    aln_str = aln_str.upper()

    if debug:
        print('cut pos: ' + str(cut_pos))
        print('ref: ' + ref_str + '\naln: ' + aln_str + '\nclipped: left: ' + str(clipped_left_bp) + ' right: ' + str(clipped_right_bp))

    left_all_match_count = 0
    left_start_match_count = 0
    left_read_bases_count = 0
    left_ref_bases_count = 0

    seen_read = False
    seen_mismatch = False
    aln_pos = 0 # index in alignment
    ref_pos = 0 # index of ref sequence (to check against cut_pos)
    while ref_pos < cut_pos:
        if ref_str[aln_pos] != '-':
            ref_pos += 1
        if aln_str[aln_pos] != '-':
            seen_read = True
            left_read_bases_count += 1
        if seen_read:
            if ref_str[aln_pos] == aln_str[aln_pos] or ignore_n:
                left_all_match_count += 1
                if not seen_mismatch:
                    left_start_match_count += 1
            else:
                seen_mismatch = True
            if ref_str[aln_pos] != '-':
                left_ref_bases_count += 1
        aln_pos += 1

    if debug:
        print('left all match count: ' + str(left_all_match_count))
        print('left start match count: ' + str(left_start_match_count))
        print('left read bases count: ' + str(left_read_bases_count))
        print('left ref bases count: ' + str(left_ref_bases_count))

    right_all_match_count = 0
    right_start_match_count = 0
    right_read_bases_count = 0
    right_ref_bases_count = 0

    seen_read = False
    seen_mismatch = False
    aln_pos = len(aln_str)-1 # index in alignment
    ref_pos = len(ref_str.replace("-","")) # index of ref sequence (to check against cut_pos)
    while ref_pos > cut_pos:
        if ref_str[aln_pos] != '-':
            ref_pos -= 1
        if aln_str[aln_pos] != '-':
            seen_read = True
            right_read_bases_count += 1
        if seen_read:
            if ref_str[aln_pos] == aln_str[aln_pos] or ignore_n:
                right_all_match_count += 1
                if not seen_mismatch:
                    right_start_match_count += 1
            else:
                seen_mismatch = True
            if ref_str[aln_pos] != '-':
                right_ref_bases_count += 1
        aln_pos -= 1

    if debug:
        print('right all match count: ' + str(right_all_match_count))
        print('right start match count: ' + str(right_start_match_count))
        print('right read bases count: ' + str(right_read_bases_count))
        print('right ref bases count: ' + str(right_ref_bases_count))

    if clipped_left_bp > 0: #if the left side of the read was soft/hard clipped
        if left_read_bases_count > 0: #if any part of this read was on the left side
            left_read_bases_count += clipped_left_bp
            left_start_match_count = 0
        else: #if this read was actually all on the right side
            right_read_bases_count += clipped_left_bp
    if clipped_right_bp > 0: #if the right side of the read was soft/hard clipped
        if right_read_bases_count > 0: #if any part of this read was on the right side
            right_read_bases_count += clipped_right_bp
            right_start_match_count = 0
        else: #if this read was actually all on the left side
            left_read_bases_count += clipped_right_bp

    return (left_read_bases_count, left_ref_bases_count, left_all_match_count, left_start_match_count,
        right_read_bases_count, right_ref_bases_count, right_all_match_count, right_start_match_count)


class PlotObject:
    """
    Holds information for plots for future output, namely:
        the plot name: root of plot (name.pdf and name.png should exist)
        the plot title: title to be shown to user
        the plot label: label to be shown under the plot
        the plot data: array of (tuple of display name and file name)
        the plot order: int specifying the order to display on the report (lower numbers are plotted first, followed by higher numbers)
    """
    def __init__(self,plot_name,plot_title,plot_label,plot_datas,plot_order=50):
        self.name = plot_name
        self.title = plot_title
        self.label = plot_label
        self.datas = plot_datas
        self.order = plot_order

    def to_json(self):
        obj = {
                'plot_name':self.name,
                'plot_title':self.title,
                'plot_label':self.label,
                'plot_datas':self.datas,
                'plot_order':self.order
                }
        obj_str = json.dumps(obj,separators=(',',':'))
        return obj_str

    #construct from a json string
    @classmethod
    def from_json(cls, json_str):
        obj = json.loads(json_str)
        return cls(plot_name=obj['plot_name'],
                plot_title=obj['plot_title'],
                plot_label=obj['plot_label'],
                plot_datas=obj['plot_datas'],
                plot_order=obj['plot_order'])




def make_report(report_file,report_name,crispratac_folder,
            crispresso_run_names,crispresso_sub_html_files,
            summary_plot_objects=[]
        ):
    """
    Makes an HTML report for a CRISPRatac run

    Parameters:
    report_file: path to the output report
    report_name: description of report type to be shown at top of report
    crispratac_folder (string): absolute path to the crispratac output

    crispresso_run_names (arr of strings): names of crispresso runs
    crispresso_sub_html_files (dict): dict of run_name->file_loc

    summary_plot_objects (list): list of PlotObjects to plot
    """

    logger = logging.getLogger('CRISPRatac')
    ordered_plot_objects = sorted(summary_plot_objects,key=lambda x: x.order)

    html_str = """
<!doctype html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
    <title>"""+report_name+"""</title>

    <!-- Bootstrap core CSS -->
<link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootswatch/4.5.2/cosmo/bootstrap.min.css">
  </head>

  <body>
<style>
html,
body {
  height: 100%;
}

body {
  padding-top: 40px;
  padding-bottom: 40px;
  background-color: #f5f5f5;
}

</style>
<div class='container'>
<div class='row justify-content-md-center'>
<div class='col-8'>
    <div class='text-center pb-4'>
    <h1 class='display-3'>CRISPRatac</h1><hr><h2>"""+report_name+"""</h2>
    </div>
"""
    data_path = ""
    if len(crispresso_run_names) > 0:
        run_string = """<div class='card text-center mb-2'>
          <div class='card-header'>
            <h5>CRISPResso Output</h5>
          </div>
          <div class='card-body p-0'>
            <div class="list-group list-group-flush">
            """
        for crispresso_run_name in crispresso_run_names:
            crispresso_run_names,crispresso_sub_html_files,
            run_string += "<a href='"+data_path+crispresso_sub_html_files[crispresso_run_name]+"' class='list-group-item list-group-item-action'>"+crispresso_run_name+"</a>\n"
        run_string += "</div></div></div>"
        html_str += run_string

    for plot_obj in ordered_plot_objects:
        plot_str = "<div class='card text-center mb-2'>\n\t<div class='card-header'>\n"
        plot_str += "<h5>"+plot_obj.title+"</h5>\n"
        plot_str += "</div>\n"
        plot_str += "<div class='card-body'>\n"
        plot_str += "<a href='"+data_path+plot_obj.name+".pdf'><img src='"+data_path + plot_obj.name + ".png' width='80%' ></a>\n"
        plot_str += "<label>"+plot_obj.label+"</label>\n"
        for (plot_data_label,plot_data_path) in plot_obj.datas:
            plot_str += "<p class='m-0'><small>Data: <a href='"+data_path+plot_data_path+"'>" + plot_data_label + "</a></small></p>\n";
        plot_str += "</div></div>\n";
        html_str += plot_str

    html_str += """
                </div>
            </div>
        </div>
    </body>
</html>
"""
    with open(report_file,'w') as fo:
        fo.write(html_str)
    logger.info('Wrote ' + report_file)

main()
