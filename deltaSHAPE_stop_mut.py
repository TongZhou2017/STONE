import sys
import os
import argparse
import numpy as np
import scipy.stats as stats
from operator import itemgetter
from itertools import groupby
import warnings
import pandas as pd

def open_csv(filename, front, back):
    # Read the CSV file
    df = pd.read_csv(filename)
    
    # Check if the file contains fewer than 4 columns
    if df.shape[1] < 4:
        sys.exit(f'ERROR: Input file {filename} contains fewer than 4 data columns. Please provide a file in .csv format with columns: Nucleotide, Reactivity, Standard Error, Nucleotide Type. See the README for details.')

    # Extract columns and convert to arrays
    data = df['Reactivity'].to_numpy()
    errs = df['Standard Error'].to_numpy()
    sequence = ''.join(df['Nucleotide Type'].astype(str))
    
    # Mask 3' primer nucleotides
    data[-back:] = -999
    errs[-back:] = 0
    
    # Mask 5' nucleotides
    data[:front] = -999
    errs[:front] = 0
    
    # Convert reactivity values of -999 to 'NaN'
    mask = data == -999
    data[mask] = np.nan
    errs[mask] = np.nan
    
    return data, errs, sequence

def smooth(data, err, pad):
    new_data, new_err = [], []
    mask = []
    
    # Create a list to store positions to ignore
    for i in range(len(data)):
        if data[i] == -999 or np.isnan(data[i]):
            mask.append(i)
    
    # Handle nucleotides at the start where the window cannot be centered
    for i in range(pad):
        new_data.append(np.nan)
        new_err.append(np.nan)
    
    # Smooth the data using a window
    for i in range(pad, len(data) - pad):
        window = data[i - pad:i + pad + 1]
        new_data.append(np.mean(np.ma.MaskedArray(window, np.isnan(window))))
        
        errs = np.array(err[i - pad:i + pad + 1])
        squerrs = np.power(errs[~np.isnan(errs)], 2)
        total = np.sum(squerrs)
        sqrt = np.sqrt(total)
        new_err.append(sqrt / len(window))
    
    # Handle nucleotides at the end where the window cannot be centered
    for i in range(pad):
        new_data.append(np.nan)
        new_err.append(np.nan)
    
    # Set masked positions to NaN
    for i in mask:
        new_data[i] = np.nan
        new_err[i] = np.nan
    
    return np.array(new_data), np.array(new_err)

def z_factor(data1, data2, err1, err2, factor=1.96):
    z_factors = []
    for i in range(len(data1)):
        if np.isnan(data1[i]) or np.isnan(data2[i]):
            z_factors.append(np.nan)
        else:
            top = factor * (err2[i] + err1[i])
            bot = abs(data2[i] - data1[i])
            if bot == 0:
                z_factors.append(np.nan)
            else:
                z = 1 - (top / bot)
                z_factors.append(z)
    return z_factors

def calc_zScores(diffs):
    mean = np.nanmean(diffs)
    sigma = np.nanstd(diffs)
    z_scores = (diffs - mean) / sigma
    return np.array(z_scores)

if __name__ == '__main__':
    
    # Pre-parse arguments to avoid argparse interpreting negatives as flags
    for i, arg in enumerate(sys.argv):
        if (arg[0] == '-') and arg[1].isdigit():
            sys.argv[i] = ' ' + arg
    
    # Parse arguments
    parse = argparse.ArgumentParser(
        description="deltaSHAPE computes statistically significant changes in SHAPE-MaP reactivity between two conditions. See README file for further details and file descriptions.", 
        epilog="deltaSHAPE v0.91 by Matt Smola (matt.smola@gmail.com)",
        add_help=False)
    
    required = parse.add_argument_group('Required files', 'These files are required in order to run deltaSHAPE analysis.')
    required.add_argument('csvFile1', type=str, help='SHAPE-MaP .csv file for the first comparison')
    required.add_argument('csvFile2', type=str, help='SHAPE-MaP .csv file, values in this file will be subtracted from those in csvFile1 for the first comparison')
    required.add_argument('csvFile3', type=str, help='SHAPE-MaP .csv file for the second comparison')
    required.add_argument('csvFile4', type=str, help='SHAPE-MaP .csv file, values in this file will be subtracted from those in csvFile3 for the second comparison')
    
    data_opt = parse.add_argument_group('Data manipulation', 'Options to specify how SHAPE-MaP data are manipulated and analyzed.')
    data_opt.add_argument('--mask5', type=int, default=0, help="Specify the number of nucleotides at the 5' end to ignore. Default: 0")
    data_opt.add_argument('--mask3', type=int, default=0, help="Specify the number of nucleotides at the 3' end to ignore. Default: 0")
    data_opt.add_argument('-p', '--pad', type=int, default=1, help='Indicate the smoothing window size. Window = 2*pad+1. To turn off smoothing, set PAD = 0. Default: 1')
    data_opt.add_argument('-z', '--Zcoeff', type=float, default=1.96, help='Adjust the Z-factor stringency by changing the equation coefficient. See the README for details. Default: 1.96')
    data_opt.add_argument('-t', '--Zthresh', type=float, default=0, help='Adjust the Z-factor stringency by changing the cutoff threshold. See the README for details. Default: 0')
    data_opt.add_argument('-s', '--SSthresh', type=float, default=1, help='Set the cutoff threshold of standard score filtering. Default: 1.0')
    data_opt.add_argument('-f', '--FindSite', type=str, default='2,3', help='Comma-separated pair of numbers indicating the window pad size and number of required hits when finding binding sites. Default settings look for 3+ nucleotides within a 5-nucleotide window. See the README for details. Default: 2,3')
    
    out_opt = parse.add_argument_group('Output', 'Options specifying plotting and output file details.')
    out_opt.add_argument('-o', '--out', type=str, default="differences.txt", help='Name and location of output file to be written. Default: ./differences.txt')
    out_opt.add_argument('--magrank', action='store_true', help='Sort output file by decreasing deltaSHAPE magnitude. Default: OFF')
    out_opt.add_argument('--all', action='store_true', help='Output data for all nucleotides. Insignificant changes are listed as zero. Default: OFF')
    out_opt.add_argument('--pdf', action='store_true', help='Save plot as PDF. If output file is given, PDF will have same prefix. Default: OFF')
    out_opt.add_argument('--noshow', action='store_true', help='Generate the plot but do not show it. Typically used with --pdf. Default: display plot')
    out_opt.add_argument('--noplot', action='store_true', help='Skip plotting completely. Default: OFF')
    out_opt.add_argument('--dots', action='store_true', help='Plot markers indicating nucleotides that pass Z-factor and standard score filtering. This can get unwieldy for large RNAs (>1000). Standard score (open) dots are plotted above Z-factor (filled) dots. Default: OFF')
    out_opt.add_argument('--Zdots', action='store_true', help='Plot markers indicating only nucleotides that pass Z-factor filtering. Default: OFF')
    out_opt.add_argument('--SSdots', action='store_true', help='Plot markers indicating only nucleotides that pass standard score filtering. Default: OFF')
    out_opt.add_argument('--colorfill', action='store_true', help='Highlight deltaSHAPE sites with coloration beneath the plot line for "prettier" figures. Default: OFF')
    out_opt.add_argument('--ymin', type=float, default=-999, help='Set plot y-axis minimum. Default: Determined automatically')
    out_opt.add_argument('--ymax', type=float, default=-999, help='Set plot y-axis maximum. Default: Determined automatically')
    out_opt.add_argument('--xmin', type=float, default=-999, help='Set plot x-axis minimum. Default: Determined automatically')
    out_opt.add_argument('--xmax', type=float, default=-999, help='Set plot x-axis maximum. Default: Determined automatically')
    
    help_opt = parse.add_argument_group('Help')
    help_opt.add_argument('-h', '--help', action="help", help="show this help message and exit")
    args = parse.parse_args()
    
    
    # Set and check analysis parameters
    front = args.mask5
    back = args.mask3
    pad = args.pad
    z_coeff = args.Zcoeff
    z_thresh = args.Zthresh
    if z_thresh > 1:
        sys.exit('ERROR: Z-factor can never exceed 1. Change -t/--Zthresh value accordingly.')
    ss_thresh = args.SSthresh
    site_pad = int(args.FindSite.split(',')[0])
    site_min = int(args.FindSite.split(',')[1])
    if site_pad * 2 + 1 < site_min:
        sys.exit('ERROR: Binding site window size and hit minimum are incompatible.\nDouble-check the -f --FindSite flag or consult the README file.')
    
    color = 'bar' if not args.colorfill else 'fill'
    outfile = os.path.normpath(args.out)
    if args.pdf:
        pdf_file = str(os.path.normpath(args.out)).split('.')[0]+".pdf"
    if args.ymin > args.ymax:
        sys.exit('ERROR: --ymin must be less than --ymax.')
    if args.xmin > args.xmax:
        sys.exit('ERROR: --xmin must be less than --xmax.')
    
    
    # Run analysis steps
    # Step 1: Open the .csv files
    data1, err1, seq1 = open_csv(args.csvFile1, front, back)
    data2, err2, seq2 = open_csv(args.csvFile2, front, back)
    data3, err3, seq3 = open_csv(args.csvFile3, front, back)
    data4, err4, seq4 = open_csv(args.csvFile4, front, back)
    
    # Step 2: Smooth data and errors
    s_data1, s_err1 = smooth(data1, err1, pad)
    s_data2, s_err2 = smooth(data2, err2, pad)
    s_data3, s_err3 = smooth(data3, err3, pad)
    s_data4, s_err4 = smooth(data4, err4, pad)
    
    # Step 3: Calculate data differences
    diff1 = data1 - data2
    diff2 = data3 - data4
    s_diff1 = smooth(diff1, err1, pad)[0]
    s_diff2 = smooth(diff2, err3, pad)[0]

    # Step 4: Calculate Z-factors
    z_factors1 = z_factor(s_data1, s_data2, s_err1, s_err2, z_coeff)
    z_factors2 = z_factor(s_data3, s_data4, s_err3, s_err4, z_coeff)

    # Step 5: Calculate Z-scores
    z_scores1 = calc_zScores(s_diff1)
    z_scores2 = calc_zScores(s_diff2)

    # Step 6: Identify significant differences (5-nt windows)
    sigdiff1, sigdiff2 = [], []
    for i in range(site_pad, len(diff1) - site_pad):
        win1 = range(i - site_pad, i + site_pad + 1)
        win2 = range(i - site_pad, i + site_pad + 1)
        count1 = count2 = 0
        maybes1 = []
        maybes2 = []
        for j in win1:
            if z_factors1[j] > z_thresh and np.abs(z_scores1[j]) >= ss_thresh:
                count1 += 1
                maybes1.append(j)
        for j in win2:
            if z_factors2[j] > z_thresh and np.abs(z_scores2[j]) >= ss_thresh:
                count2 += 1
                maybes2.append(j)
        if count1 >= site_min:
            for k in maybes1:
                if k not in sigdiff1:
                    sigdiff1.append(k)
        if count2 >= site_min:
            for k in maybes2:
                if k not in sigdiff2:
                    sigdiff2.append(k)

    # Step 7: Prepare for plotting and output data
    pos_consec1, neg_consec1 = [], []
    pos_consec2, neg_consec2 = [], []
    for k, g in groupby(enumerate([i for i in sigdiff1 if s_diff1[i] >= 0]), lambda ix: ix[0] - ix[1]):
        pos_consec1.append(list(map(itemgetter(1), g)))
    for k, g in groupby(enumerate([i for i in sigdiff1 if s_diff1[i] < 0]), lambda ix: ix[0] - ix[1]):
        neg_consec1.append(list(map(itemgetter(1), g)))
    for k, g in groupby(enumerate([i for i in sigdiff2 if s_diff2[i] >= 0]), lambda ix: ix[0] - ix[1]):
        pos_consec2.append(list(map(itemgetter(1), g)))
    for k, g in groupby(enumerate([i for i in sigdiff2 if s_diff2[i] < 0]), lambda ix: ix[0] - ix[1]):
        neg_consec2.append(list(map(itemgetter(1), g)))

    pos_shade_bits1, pos_x_bits1 = [], []
    pos_shade_bits2, pos_x_bits2 = [], []
    pos_span1, pos_span2 = [], []
    data_out1, data_out2 = [], []
    for region in pos_consec1:
        pos_shade, pos_x = [], []
        for i in region:
            pos_shade.extend((s_diff1[i], s_diff1[i]))
            pos_x.extend((i + 0.5, i + 1.5))
            data_out1.append([i + 1, seq1[i], s_diff1[i], z_factors1[i], z_scores1[i], s_data1[i], s_data2[i], diff1[i], data1[i], data2[i]])
        pos_shade_bits1.append(pos_shade)
        pos_x_bits1.append(pos_x)
        pos_span1.append([region[0] + 0.5, region[-1] + 1.5])
    
    for region in pos_consec2:
        pos_shade, pos_x = [], []
        for i in region:
            pos_shade.extend((s_diff2[i], s_diff2[i]))
            pos_x.extend((i + 0.5, i + 1.5))
            data_out2.append([i + 1, seq3[i], s_diff2[i], z_factors2[i], z_scores2[i], s_data3[i], s_data4[i], diff2[i], data3[i], data4[i]])
        pos_shade_bits2.append(pos_shade)
        pos_x_bits2.append(pos_x)
        pos_span2.append([region[0] + 0.5, region[-1] + 1.5])
    
    neg_shade_bits1, neg_x_bits1 = [], []
    neg_shade_bits2, neg_x_bits2 = [], []
    neg_span1, neg_span2 = [], []
    for region in neg_consec1:
        neg_shade, neg_x = [], []
        for i in region:
            neg_shade.extend((s_diff1[i], s_diff1[i]))
            neg_x.extend((i + 0.5, i + 1.5))
            data_out1.append([i + 1, seq1[i], s_diff1[i], z_factors1[i], z_scores1[i], s_data1[i], s_data2[i], diff1[i], data1[i], data2[i]])
        neg_shade_bits1.append(neg_shade)
        neg_x_bits1.append(neg_x)
        neg_span1.append([region[0] + 0.5, region[-1] + 1.5])
    
    for region in neg_consec2:
        neg_shade, neg_x = [], []
        for i in region:
            neg_shade.extend((s_diff2[i], s_diff2[i]))
            neg_x.extend((i + 0.5, i + 1.5))
            data_out2.append([i + 1, seq3[i], s_diff2[i], z_factors2[i], z_scores2[i], s_data3[i], s_data4[i], diff2[i], data3[i], data4[i]])
        neg_shade_bits2.append(neg_shade)
        neg_x_bits2.append(neg_x)
        neg_span2.append([region[0] + 0.5, region[-1] + 1.5])
    
# Step 8: Plotting
if not args.noplot:
    
    if args.noshow:
        import matplotlib
        matplotlib.use('Agg')
    
    import matplotlib.pyplot as plt

    plt.figure(figsize=(15, 6))
    x1 = range(1, len(s_diff1) + 1)
    x2 = range(1, len(s_diff2) + 1)
    
    plt.plot(x1, s_diff1, drawstyle='steps-mid', color='#A9AD9B', label='STOP-signal')
    plt.plot(x2, s_diff2, drawstyle='steps-mid', color='#BD8185', label='MUT-signal')
    plt.axhline(0, color='black')
    # Adding legend in the lower right corner
    plt.legend(loc='lower right')
    # Mask primer binding regions
    plt.axvspan(0, front + 0.5, color="grey", alpha=0.25)
    plt.axvspan(len(diff1) - back + 0.5, len(diff1) + 0.5, color="grey", alpha=0.25)

    # Shade deltaSHAPE sites
    if color == 'fill':
        for i in range(len(pos_shade_bits1)):
            plt.fill_between(pos_x_bits1[i], 0, np.maximum(pos_shade_bits1[i], 0), color='#7C9C29', step='mid', alpha=1)
        for i in range(len(neg_shade_bits1)):
            plt.fill_between(neg_x_bits1[i], 0, np.minimum(neg_shade_bits1[i], 0), color='#B5C583', step='mid', alpha=1)
        for i in range(len(pos_shade_bits2)):
            plt.fill_between(pos_x_bits2[i], 0, np.maximum(pos_shade_bits2[i], 0), color='#F19465', step='mid', alpha=1)
        for i in range(len(neg_shade_bits2)):
            plt.fill_between(neg_x_bits2[i], 0, np.minimum(neg_shade_bits2[i], 0), color='#FFCA99', step='mid', alpha=1)
    elif color == 'bar':
        for i in range(len(pos_shade_bits1)):
            plt.fill_between(pos_x_bits1[i], 0, np.maximum(pos_shade_bits1[i], 0), color='#7C9C29', step='mid', alpha=1)
        for i in range(len(neg_shade_bits1)):
            plt.fill_between(neg_x_bits1[i], 0, np.minimum(neg_shade_bits1[i], 0), color='#B5C583', step='mid', alpha=1)
        for i in range(len(pos_shade_bits2)):
            plt.fill_between(pos_x_bits2[i], 0, np.maximum(pos_shade_bits2[i], 0), color='#F19465', step='mid', alpha=1)
        for i in range(len(neg_shade_bits2)):
            plt.fill_between(neg_x_bits2[i], 0, np.minimum(neg_shade_bits2[i], 0), color='#FFCA99', step='mid', alpha=1)

    # Whether to plot Z-factor and standard score dots
    dots = args.dots or args.SSdots or args.Zdots
    
    # Set axis range
    if args.xmax == -999:
        args.xmax = max(len(diff1), len(diff2))
    if args.xmin == -999:
        args.xmin = 1
            
    plt.xlim(args.xmin, args.xmax)
    
    if args.ymin == -999:
        y_min = min(filter(lambda x: not np.isnan(x), s_diff1.tolist() + s_diff2.tolist())) - 0.25
    else:
        y_min = args.ymin

    if args.ymax == -999:
        y_max = max(filter(lambda x: not np.isnan(x), s_diff1.tolist() + s_diff2.tolist())) + 0.25
        if dots:
            y_max += 0.6
    else:
        y_max = args.ymax
    
    plt.ylim(y_min, y_max)
    
    # Plot Z-factor/standard score dots
    if dots:
        for i in range(len(z_scores1)):
            if (args.dots or args.SSdots) and abs(z_scores1[i]) >= ss_thresh:
                plt.scatter(i + 1, y_max - 0.2, marker="s", s=5, color='none', edgecolor='black', zorder=3)
            if (args.dots or args.Zdots) and z_factors1[i] >= z_thresh:
                plt.scatter(i + 1, y_max - 0.4, marker="o", s=5, color='black', zorder=3)
        for i in range(len(z_scores2)):
            if (args.dots or args.SSdots) and abs(z_scores2[i]) >= ss_thresh:
                plt.scatter(i + 1, y_max - 0.2, marker="s", s=5, color='none', edgecolor='red', zorder=3)
            if (args.dots or args.Zdots) and z_factors2[i] >= z_thresh:
                plt.scatter(i + 1, y_max - 0.4, marker="o", s=5, color='red', zorder=3)

    # Labels and tick formatting
    plt.xlabel("Nucleotide")
    plt.ylabel(r'$\Delta$SHAPE')
    plt.tick_params(which='both', direction='out', top=False, right=False)
    plt.legend(loc='best')
    
    # Temporarily disable UserWarnings
    warnings.simplefilter("ignore", UserWarning)
    plt.tight_layout()
    warnings.resetwarnings()
    
    if args.pdf:
        plt.savefig(pdf_file, format='pdf')
    if not args.noshow:
        plt.show()
    
    # Step 9: Generate output files
    if args.all:
        for i in filter(lambda x: x + 1 not in [j[0] for j in data_out1], range(len(seq1))):
            data_out1.append([i + 1, seq1[i], 0, z_factors1[i], z_scores1[i], s_data1[i], s_data2[i], diff1[i], data1[i], data2[i]])
        for i in filter(lambda x: x + 1 not in [j[0] for j in data_out2], range(len(seq3))):
            data_out2.append([i + 1, seq3[i], 0, z_factors2[i], z_scores2[i], s_data3[i], s_data4[i], diff2[i], data3[i], data4[i]])
    
    with open(outfile, 'w') as o:
        o.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('Nuc', 'Seq', 'DeltaSHAPE', 'Z-factor', 'Std_Score', 'Smoothed_Data1', 'Smoothed_Data2', 'Unsmoothed_Diff', 'Data1', 'Data2'))
        
        if args.magrank:
            data_out1.sort(reverse=True, key=lambda x: abs(x[2]))
            data_out2.sort(reverse=True, key=lambda x: abs(x[2]))
        else:
            data_out1.sort(key=lambda x: x[0])
            data_out2.sort(key=lambda x: x[0])
        
        for i in data_out1:
            for j in range(2, len(i)):
                if np.isnan(i[j]):
                    i[j] = -999
            o.write('\t'.join(map(str, i)) + "\n")
        o.write('\n')  # Add a newline between two results
        for i in data_out2:
            for j in range(2, len(i)):
                if np.isnan(i[j]):
                    i[j] = -999
            o.write('\t'.join(map(str, i)) + "\n")
