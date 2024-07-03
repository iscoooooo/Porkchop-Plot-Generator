def stateReader(filespec):
    # This function assumes that the loaded Ephmeris Data is in comma separated value (csv) text format.

    expr1 = 'SOE'
    expr2 = 'EOE'

    # Loop 1: to find where the data range is at
    N = 0 # line count

    with open(filespec,'r') as file:
        while True:
            dataline = file.readline() # read one line at a time

            N = N + 1 # find the line number

            if isinstance(dataline, str): # if the line is a string

                # found the 1st string that I'm seeking
                if expr1 in dataline:
                    idx1 = N
                    print(f'idx1: {idx1}')
                
                # found the 2nd string that I'm seeking
                if expr2 in dataline:
                    idx2 = N
                    print(f'idx2: {idx2}')

            # break loop if there is no longer any text
            if not isinstance(dataline, str):
                break

filespec = 'data/Earth_Ephemeris_Data.txt'
stateReader(filespec)