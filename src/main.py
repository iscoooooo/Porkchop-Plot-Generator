from utils import stateReader

def main():
    # Define the filespec for the data file
    filespec = '../data/Earth_Ephemeris_Data.txt'
    
    # Call the stateReader function
    julianDates, states = stateReader(filespec)
    
    # Print results (or handle them as needed)
    print(f'\n Julian Dates: \n {julianDates}')
    print(f'\n\n States: \n {states}\n')

if __name__ == "__main__":
    main()