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
parser.add_argument('--type', type=str,
                    help='Data type to analyse: beam_rad, insol_time')
parser.add_argument('--relative', type=bool,
                    help='Whether to use relative data type variant: give any value for yes')
parser.add_argument('--title', type=str,
                    help='title_no to perform analysis on')
parser.add_argument('--statistics', type=str,
                    help='type of statistics to perform on data: minimum, median, or maximum value distribution')
parser.add_argument('--threshold', type=str,
                    help='Threshold value to return all titles with all values above for the given data type')
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
    for column in data.columns:
        if data_type in column and 'relative' not in column:
            prefix = column[:column.index(data_type)]
            break
    if args.relative:
        data = data.filter(regex='^relative_' + prefix + data_type + '*', axis=1)
    else:
        data = data.filter(regex='^' + prefix + data_type + '*', axis=1)
if not args.title and not args.statistics and not args.threshold:
    print('Either a --title, --statistics, or --threshold argument must be given for the script to do work. '
          'Use -h flag for help.')
    exit()

# If desired, do title-based analysis

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
    # Plot the data
    data.plot.hist()
    plt.show()
    exit()  # Exit the code here

# If desired, do statistical analysis of all title data

if args.statistics:
    method = args.statistics
    # Extract the desired values from each title's data
    if method == 'minimum':
        subset = data.min(axis=1)
    elif method == 'maximum':
        subset = data.max(axis=1)
    elif method == 'median':
        subset = data.median(axis=1)
    elif method == 'sum':
        subset = data.sum(axis=1)
        # Remove 0 value sums as erroneous
        subset = subset.loc[subset != 0]
    # Calculate inequality graph
    values, sums, exclusions = [], [], []
    for q in range(100):
        quantile_value = subset.quantile(q / 100)
        quantile_subset = subset.loc[subset <= quantile_value]
        quantile_sum = quantile_subset.sum() / subset.sum()
        quantile_exclusion = 1 - quantile_sum
        values.append(quantile_value)
        sums.append(100 * quantile_sum)
        exclusions.append(100 * quantile_exclusion)
    subset.plot.hist(bins=365)
    plt.show()
    plt.plot([0, 100], [0, 100], color='k')
    plt.plot(range(100), sums, color='b')
    plt.plot(range(100), exclusions, color='r')
    plt.xlabel('Distribution Percentage')
    plt.ylabel('Percentage of Distribution Values')
    plt.show()
    exit()  # Exit the code here

# Otherwise, do threshold-based analysis

if args.threshold:
    # First, find all titles for which the threshold is met
    threshold = float(args.threshold)
    for column in data.columns:
        data = data.loc[data[column] >= threshold]

    # Then, write to file those titles
    with open('titles_' + str(threshold) + '_' + data_type + '.csv', 'w') as outfile:
        outfile.write('title_no,\n')
    with open('titles_' + str(threshold) + '_' + data_type + '.csv', 'a') as outfile:
        for title in data.index:
            outfile.write(str(title) + ',THRESHOLD MET\n')
