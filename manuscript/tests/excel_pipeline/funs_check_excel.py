import os

def check_file_exists(file_path: str) -> bool:
    """Check if the file exists."""
    return os.path.isfile(file_path)