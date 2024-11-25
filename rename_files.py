import os
import re

def rename_files(directory):
    for filename in os.listdir(directory):
        # Create the new filename by replacing spaces and tabs with underscores
        new_filename = re.sub(r'[ \t]+', '_', filename)
        
        # Check if the filename needs to be changed
        if new_filename != filename:
            src = os.path.join(directory, filename)
            dst = os.path.join(directory, new_filename)
            
            # Rename the file
            os.rename(src, dst)
            print(f'Renamed: "{filename}" to "{new_filename}"')

if __name__ == "__main__":
    # Specify your directory path
    directory = r'D:\X\ND\Files'
    
    # Call the rename function
    rename_files(directory)
