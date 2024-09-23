import os

def get_test_data_path(relative_path):
    """
    Get the absolute path for test data files.
    Works in both local and Docker environments.
    
    :param relative_path: Relative path to the test data file
    :return: Absolute path to the test data file
    """
    # Check if we're running in Docker
    if os.path.exists('/app'):
        base_dir = '/app'
    else:
        # If not in Docker, use the directory of this file as the base
        base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    return os.path.abspath(os.path.join(base_dir, relative_path))

def get_opt_data_path(filename):
    """
    Get the path for files in the /opt/data directory.
    Works in both local and Docker environments.
    
    :param filename: Name of the file in /opt/data
    :return: Absolute path to the file
    """
    if os.path.exists('/opt/data'):
        return os.path.join('/opt/data', filename)
    else:
        # If not in Docker, assume /opt/data contents are in the project's data directory
        return get_test_data_path(os.path.join('data', filename))