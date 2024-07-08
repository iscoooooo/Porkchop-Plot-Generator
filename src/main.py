from utils import stateReader
from planetary_data import earth
from ephemeris_query import *


# def main():
#     # Define the filespec for the data file
#     filespec = '../data/Earth_Ephemeris_Data.txt'
    
#     # Call the stateReader function
#     julianDates, states = stateReader(filespec)
    
#     # Print results (or handle them as needed)
#     print(f'\n Julian Dates: \n {julianDates}')
#     print(f'\n\n States: \n {states}\n')

def main():
    # Define inputs
    ID = earth['ID']
    start_time = '2020-07-01'
    stop_time = '2020-09-19'
    step_size = '1'

    # Call function
    url = generate_url(ID, start_time, stop_time, step_size)

    # Define relative path for the output file
    output_dir = '../data/departure_data/'
    output_filename = f"ephemeris_data_{start_time}_{stop_time}.txt"
    output_path = output_dir + output_filename

    save_query_to_file(url, output_path)

if __name__ == "__main__":
    main()