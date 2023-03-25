import argparse
import os,sys
import shutil
import time
import glob
from pybedtools import BedTool
import subprocess

start_time = time.process_time()
print('Process started at ' + str(start_time))
def get_parser():
    parser = argparse.ArgumentParser(prog='methylcon', description='Amplicon BS-Seq Explorer')
    return parser

def get_args():
    parser = argparse.ArgumentParser(prog='methylcon', description='Amplicon BS-Seq Explorer')
    parser.add_argument('--inputdir', help='Location of unzipped Fastq file(s)', type=str)
    parser.add_argument('--outdir', help='Output directory', default=None, type=str)
    parser.add_argument('--softsheet', help='Path to softwares directory(ies)', default=None, type=str)
    parser.add_argument('--samplesheet', help='Path to sample information(s)', default=None, type=str)
    parser.add_argument('--genome', help='Path to genome folder', default=None, type=str)
    parser.add_argument('--genome_prefix', help='genome prefix', default=None, type=str)
    parser.add_argument('--Rscript', help='Path to the Rscript', default=None, type=str)
    parser.add_argument('--plotscript', help='Path to the plot script', default=None, type=str)
    return parser.parse_args()


def softsheet(args):
    softsdict= {}
    softsheet = open(args.softsheet)
    softsheet = softsheet.read().strip().split('\n')
    for softs in softsheet:
        softs = softs.strip().split(",")
        softsdict[softs[0]] = softs[1]
    return softsdict

def samplesheet(args):
    samplesheetdict={}
    samplesheet = open(args.samplesheet)
    samplesheet = samplesheet.read().strip().split('\n')
    for samples in samplesheet:
        samples = samples.strip().split(",")
        samplesheetdict[samples[0]] = [samples[1], samples[2]]
    return samplesheetdict


def run_tools(args, softs, samples):
    # create an output directory
    if os.path.exists(args.outdir):
        shutil.rmtree(args.outdir)
    else:
        os.mkdir(args.outdir)
        # Get subset of reference genome (multiple location can be supplied in the primer_site.bed file)
        print('Creating subset of reference genome \n')
        os.system(softs["bedtools"] + 'bedtools getfasta -fi ' + args.genome_prefix + '.fa ' + '-bed primer_site.bed -fo primer_site.fa')
        if os.path.exists('genome'):
            pass
        else:
            os.mkdir('genome')
        primer_site = open("primer_site.bed")
        primer_site = primer_site.read().strip().split('\n')
        for sites in primer_site:
            sites = sites.strip().split("\t")
            print(sites)
            os.system('samtools faidx ' + args.genome_prefix + '.fa ' + sites[0] + ' -n 50 > ' + sites[0] + '.fa')
            os.system('mv ' + sites[0] + '.fa* genome/')
        print('Subset of reference genome created \n')
        if os.path.exists(args.genome + 'Bisulfite_Genome'):
            print('Bisulfite Genome Indices were found\n')
            pass
        else:
            print('Genome preparation started\n')
            os.system(softs["bismark"] + 'bismark_genome_preparation --bowtie2 --verbose --path_to_aligner ' + softs["bowtie2"] + ' ' + args.genome)

        for input in samples:
            print(input, samples[input][0], samples[input][1])
            os.system('cp ' + samples[input][0] + ' ' + args.outdir + '/' + input + '_R1.fastq.gz')
            os.system('cp ' + samples[input][1] + ' ' + args.outdir + '/' + input + '_R2.fastq.gz')
            os.chdir(args.outdir)
            print('Read 1: ' + input + '_R1.fastq.gz')
            os.system('zgrep @ ' + input + '_R1.fastq.gz' + ' | wc -l')
            print('Read 2: ' + input + '_R2.fastq.gz')
            os.system('zgrep @ ' + input + '_R2.fastq.gz' + ' | wc -l')
            if os.path.exists('fastqc'):
                pass
            else:
                os.mkdir('fastqc')
            print('\nFastQC for Read1 ' + input + '_R1.fastq.gz' + '\n')
            os.system(softs["fastqc"] + 'fastqc ' + input + '_R1.fastq.gz')
            print('\nFastQC for Read2 ' + input + '_R2.fastq.gz' + '\n')
            os.system(softs["fastqc"] + 'fastqc ' + input + '_R2.fastq.gz')
            os.system('mv *fastqc* fastqc/')
            if os.path.exists('trimgalore'):
                pass
            else:
                os.mkdir('trimgalore')
            print ('Trimming data..: --phred33  (Sanger/Illumina 1.9+ encoding) for quality trimming. Default: ON, -q 20 (default),  BS-Seq' + input +'\n')
            os.system(softs["trimgalore"] + 'trim_galore --phred33 --length 36 -q 20 --paired ' + input + '_R1.fastq.gz ' + input + '_R2.fastq.gz ' + '--path_to_cutadapt ' + softs["cutadapt"] + 'cutadapt')
            os.system('mv *trimming_report.txt trimgalore/')
            if os.path.exists('bismark'):
                pass
            else:
                os.mkdir('bismark')
            print('Aligning data using Bismark: ' + input + '\n')
            os.system(softs["bismark"] + 'bismark --score_min L,0,-0.6 --genome ' + args.genome + ' -1 ' + input + '_R1_val_1.fq.gz -2 ' + input + '_R2_val_2.fq.gz --bowtie2 --path_to_bowtie2 ' + softs["bowtie2"])

            print('As the sample is amplicon duplicates are not removed')

            print('Give a ID to BAM...\n')
            os.system('cp ' + input + '_R1_val_1_bismark_bt2_pe.bam ' + input + '_bismark_bt2_pe.bam')

            print('Sorting data...\n')
            os.system(softs["samtools"] + 'samtools sort -o ' + input + '_bismark_bt2_pe.sort.bam ' + input + '_bismark_bt2_pe.bam')

            print('Index data...\n')
            os.system(softs["samtools"] + 'samtools index ' + input + '_bismark_bt2_pe.sort.bam')

            print('Sorting data by readname...\n')
            os.system(softs["samtools"] + 'samtools sort -n ' + input + '_bismark_bt2_pe.bam -o ' + input + '_bismark_bt2_pe_sortedByReadname.bam')

            print('BAM to fasta: Can be imported to BiQ Analyzer\n')
            os.system(softs["samtools"] + 'samtools fasta -1 ' + input + '_R1.fa -2 ' + input + '_R2.fa ' + input + '_bismark_bt2_pe.sort.bam')

            print('Methylation extraction...\n')
            os.system(softs["bismark"] + 'bismark_methylation_extractor --bedGraph -p --include_overlap --zero_based --cutoff 10 ' + input + '_bismark_bt2_pe_sortedByReadname.bam ' + '--samtools_path ' + softs["samtools"] + ' --cytosine_report --genome_folder ' + args.genome)

            os.system('mv *bismark*.txt bismark/')
            os.system('gunzip *.bedGraph.gz')

            os.chdir('..')



def get_graphics(args):
    os.chdir(args.outdir)
    if os.path.exists('merged_bedGraph.bdg'):
        os.remove('merged_bedGraph.bdg')
    #Merge bedgraphs
    for bdg in glob.glob('*.bedGraph'):
        print(bdg)
        bdgname = bdg.replace('.bedGraph', '')
        bdgname = bdgname.split('_')[0]
        mybdg = open(bdg, 'r')
        mybdg = mybdg.read().strip().split("\n")
        mybdg = mybdg[1:]
        for coordinates in mybdg:
            coordinates = coordinates.split("\t")
            with open('merged_bedGraph.bdg', 'a') as mybdgout:
                mybdgout.write(''.join(coordinates[0] + '\t' + coordinates[1] + '\t' + coordinates[2] + '\t' + coordinates[3] + '\t' + bdgname + '\n'))
                mybdgout.close()

    # Run bedtools
    if os.path.exists('merged_bedGraph_onsite.txt'):
        os.remove('merged_bedGraph_onsite.txt')
    merged_bedGraph = open("merged_bedGraph.bdg")
    primer_site = open("./../primer_site.bed")
    merged_bedGraph = BedTool(merged_bedGraph)
    primer_site = BedTool(primer_site)
    merged_bedGraph_onsite = merged_bedGraph.intersect(primer_site,wa=True, wb=True).moveto('merged_bedGraph_onsite.txt')
    print(merged_bedGraph_onsite)
    os.chdir('..')
    os.system('cp ' + args.outdir + '/merged_bedGraph_onsite.txt ./')

    # plot Line chart
    # set the path to the Rscript executable and the path to the R script
    RSCRIPT_PATH = args.Rscript
    R_FILE_PATH = args.plotscript

    # define the command to run the R script
    Rcommand = [RSCRIPT_PATH, R_FILE_PATH]
    Rresult = subprocess.run(Rcommand, capture_output=True, text=True)
    print(Rresult.stdout)

# Print help message if no arguments
myparser = get_parser()
if len(sys.argv)==1:
   myparser.print_help(sys.stderr)
   sys.exit(1)

#Get parameters
args1 = get_args()

#Get softwares
softsall = softsheet(args1)

#Get samples information
samplesall = samplesheet(args1)

#Run the main pipeline
run_tools(args1, softsall, samplesall)

#Plot charts
get_graphics(args1)

print('\n')
print('Analysis Done. Please check the files !')
print('- > BAM files can be uploaded to IGV for visualization.')
print('- > bedGraph files can be uploaded to UCSC or IGV for visualization.')
print('- > cov files can be imported to SeqMonk for quantitation of methylation.')
print('- > fasta files can be imported to BiQ Analyzer for Heatmap of Individual sites.')
print('\n')

#time taken
print('----------------------------------------------------')

print('Process finished at')
time_taken = float(time.process_time() - start_time)
print(time_taken)
hours = time_taken/3600
time_taken = time_taken - 3600*hours
minutes = time_taken/60
seconds = time_taken - 60*minutes

print('%d:%d:%d' %(hours,minutes,seconds))
print('----------------------------------------------------')


