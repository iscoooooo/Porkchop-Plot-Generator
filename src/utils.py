import csv
import numpy as np

def stateReader(filespec):
    # This function assumes that the loaded Ephmeris Data is in comma separated value (csv) text format.

    # Initialize lists to store Julian Dates and states
    julianDates = []
    states = []

    expr1 = 'SOE' # text to look for in the file
    expr2 = 'EOE' # text to look for in the file

    # Loop to find where the data range is at
    N = 0 # line count
    idx1 = None
    idx2 = None

    with open(filespec,'r') as file:
        while True:
            dataline = file.readline() # read one line at a time

            N = N + 1 # find the line number

            if isinstance(dataline, str): # if the line is a string

                # found the 1st string that I'm seeking
                if expr1 in dataline:
                    idx1 = N
                
                # found the 2nd string that I'm seeking
                if expr2 in dataline:
                    idx2 = N

            # break loop if both indices are found
            if (idx1 is not None) and (idx2 is not None):
                break
        
        # scan through data, skipping header, and collecting data
        with open(filespec, 'r') as file:
            reader = csv.reader(file,delimiter=',')

            # skip header lines to idx1
            for _ in range(idx1):
                next(reader)

            # collect data starting from line idx1 to idx2
            for row in reader:
                if reader.line_num >= idx2:
                    break

                julianDates.append(float(row[0]))
                states.append([float(value) for value in row[2:8]])
            
        # Convert lists to numpy arrays for further processing
        julianDates = np.array(julianDates)
        states = np.array(states)

    return julianDates, states