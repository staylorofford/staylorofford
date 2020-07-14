import argparse
import matplotlib.pyplot as plt
import pandas as pd

# Parse arguments and data
parser = argparse.ArgumentParser(description='Perform analysis on sunlight data. Script can either plot distribution '
                                             'of sunlight data values for a given title or can return a variant of the'
                                             'input CSV with only those titles with data matching criteria set in the'
                                             'script input.')
parser.add_argument('--path', type=str,
                    help='Path to CSV file containing titles and sunlight data')
parser.add_argument('--title', type=str,
                    help='title_no to perform analysis on')
parser.add_argument('--type', type=str,
                    help='Data type to analyse: beam_rad, insol_time')
parser.add_argument('--relative', type=bool,
                    help='Whether to use relative data type variant')
parser.add_argument('--threshold', type=str,
                    help='Threshold value to return all titles with all (unless specified with --percentile argument)'
                         ' values above for the given data type')
parser.add_argument('--percentile', type=str,
                    help='Percentile of values above given threshold at which threshold is considered met.')
args = parser.parse_args()
if not args.path:
    print('Path to CSV file must be set with --path argument. Use -h flag for help.')
    exit()
else:
    data = pd.read_csv(args.path, index_col='title_no')
if not args.type:
    print('Data type for analysis must be set with --type argument. Use -h flag for help.')
    exit()
else:
    data_type = args.type
    if args.relative:
        data = data.filter(regex='^relative*' + data_type + '*', axis=1)
    else:
        for column in data.columns:
            if data_type in column and 'relative' not in column:
                prefix = column[:column.index(data_type)]
                break
        data = data.filter(regex='^' + prefix + data_type + '*', axis=1)
if not args.title and not args.threshold:
    print('Either a --title or --threshold argument must be given for the script to do work. Use -h flag for help.')
    exit()
if args.title:
    title = args.title
    data = data.loc[title]
    # Reformat index to be sortable time references
    reindex = []
    for value in data.index:
        reindex.append(int(value.split('_')[-2]))
    data.index = reindex
    # Sort the data by DOY
    data.sort_index(inplace=True)
if args.threshold:
    threshold = float(args.threshold)
if args.percentile and not args.threshold:
    print('Percentile argument only valid when --threshold argument is given. Use -h flag for help.')
    exit()
elif args.percentile:
    percentile = float(args.percentile)

data.plot.hist()
plt.show()
