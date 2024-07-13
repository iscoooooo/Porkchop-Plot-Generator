'''
Module containing functions related to querying the SSD Horizons API
'''

# Python Standard Libraries
import sys
import requests
import csv

# Third-Party Libraries
import numpy as np


def generate_url(ID, start_time, stop_time, step_size):
    '''
    Generates a URL for querying the Horizons API.

    Parameters:
    object_id (str): The ID of the object to query.
    start_time (str): The start time for the ephemeris data.
    stop_time (str): The stop time for the ephemeris data.
    step_size (str): The step_size for the ephemeris data.

    Returns:
        str: The generated URL. 
    '''
    base_url = "https://ssd.jpl.nasa.gov/api/horizons.api?format=text"
    params = {
        "COMMAND"       : f"'{ID}'",
        "OBJ_DATA"      : "'NO'",
        "MAKE_EPHEM"    : "'YES'",
        "EPHEM_TYPE"    : "'VECTORS'",
        "REF_PLANE"     : "'ECLIPTIC'",
        "REF_SYSTEM"    : "'J2000'",
        "VEC_TABLE"     : "'2'",
        "CSV_FORMAT"    : "'YES'",
        "CENTER"        : "'500@0'",
        "START_TIME"    : f"'{start_time}'",
        "STOP_TIME"     : f"'{stop_time}'",  
        "STEP_SIZE"     : f"'{step_size}d'", # days
        "QUANTITIES"    : "'2'"
    }

    query_string = "&".join([f"{key}={value}" for key, value in params.items()])

    return f"{base_url}&{query_string}"


def save_query_to_file(url, output_filename):
    '''
    Submits API request and saves the response text to a file.

    Parameters:
    url (str): The URL to query
    output_filename (str): The name of the output text file.
    '''

    response = requests.get(url)

    if response.status_code == 200:
        try:
            with open(output_filename, "w") as file:
                file.write(response.text)
                print(f"Ephemeris data saved to {output_filename}")
        except OSError as err:
            print(f"Unable to open file '{output_filename}'")
    else:
        print("Error: Request failed")
        print(f"Response code: {response.status_code}")
        print(response.text)
        sys.exit(1)

def stateReader(filespec):
    '''
    Reads ephemeris data from a text file in CSV format and extracts Julian dates and state vectors.

    The function expects the file to contain ephemeris data in a specific format.
    It looks for the start ("SOE") and end ("EOE") markers to determine the data range
    and extracts the data between these markers. The extracted data is then converted
    to numpy arrays for further processing.

    Parameters:
    filespec : str
        The file path to the txt file containing the ephemeris data.

    Returns:
    julianDates : ndarray
        An array of Julian dates.
    states : ndarray
        An array of state vectors, where each state vector contains position and velocity
        components (x, y, z, vx, vy, vz).
    '''

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

def encode_value(value):
    '''
    Manually encodes a string for use in a URL query parameter. This function ensures the special character `'` in query parameter values are correctly encoded.
    
    Parameters:
    value (str): The string to encode.
    
    Returns:
    str: The encoded string.
    '''
    replacements = {
        ' ': '%20',
        '\n': '%0A',
        '#': '%23',
        '$': '%24',
        '&': '%26',
        '+': '%2B',
        ',': '%2C',
        '/': '%2F',
        ':': '%3A',
        ';': '%3B',
        '=': '%3D',
        '?': '%3F',
        '@': '%40',
        '[': '%5B',
        ']': '%5D',
        "'": '%27'
    }

    for char, code in replacements.items():
        value = value.replace(char, code)
    return value