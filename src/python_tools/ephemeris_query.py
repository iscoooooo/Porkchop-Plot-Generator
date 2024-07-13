'''
Module containing functions related to querying the SSD Horizons API
'''

# Libraries
import sys
import requests


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